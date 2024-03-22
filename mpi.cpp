#define EPS 1e-14
// 19:40. Осталось добавить:
// 1). Подсчёт обоих невязок.
// 2). Обработку ошибок (если матрица не обращается типо).
// 3). Посчитать норму матрицы и нормально использовать её в коде.
// 4). Подсчёт времени работы одного процесса.
// 5). Написать makefile, чтобы программа собиралась так, как требуется.
#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <string>
#include <cmath>

int l2g(int n, int m, int p, int k, int i_loc) {
    int i_loc_m = i_loc / m;
    int i_glob_m = i_loc_m * p + k;
    return i_glob_m * m + i_loc % m;
}

// global to local, если разделение по строкам
int g2l(int n, int m, int p, int i_glob) {
    int i_glob_m = i_glob / m;
    int i_loc_m = i_glob_m / p;
    return i_loc_m * m + i_glob % m;
}

// максимальное число блочных строк на процесс
int get_max_rows(int n, int m, int p) {
    int b = (n + m - 1) / m;
    return (b + p - 1) / p;
}

// число строк на процесс k
int get_rows(int n, int m, int p, int k) {
    // число блочных строк
    int b = (n + m - 1) / m;
    return (b % p) >= k + 1 ? (b / p) + 1 : b / p; 
}

// в каком процессе лежит строка i_glob
int get_k(int n, int m, int p, int i_glob) {
    int i_glob_m = i_glob / m;
    return i_glob_m % p;
}

double f(int s, int n, int i, int j) { 
    switch(s) {
        case 1:
            return n - std::max(i, j);
        case 2:
            return std::max(i, j) + 1;
        case 3:
            return std::fabs(i-j);
        case 4:
            return 1. / (i + j + 1);
        default:
            return 1;
    }
}

void init_matrix(double* a, /* своя часть */ int n, int m, int p, int k, int s, double (*f)(int, int, int, int)) {
    int rows = get_rows(n, m, p, k); // число блочных строк в процессе k

    for (int i_loc_m = 0; i_loc_m < rows; ++i_loc_m) {
        for (int i_loc = i_loc_m * m; i_loc < (i_loc_m + 1) * m; ++i_loc) {
            int i_glob = l2g(n, m, p, k, i_loc);
            if (i_glob >= n) {
                break;
            }
            for (int j = 0; j < n; ++j) {
                a[n*i_loc + j] = f(s, n, i_glob, j);
            }
        }
    }
}

void init_b(double* a, double* b, int n, int m, int p, int k) {
    int rows = get_rows(n, m, p, k); // число блочных строк в процессе k

    for (int i_loc_m = 0; i_loc_m < rows; ++i_loc_m) {
        for (int i_loc = i_loc_m * m; i_loc < (i_loc_m + 1) * m; ++i_loc) {
            int i_glob = l2g(n, m, p, k, i_loc);
            if (i_glob >= n) {
                break;
            }

            double sum = 0;
            for (int k = 0; k <= (n-1) / 2; ++k) {
                sum += a[n * i_loc + 2*k];
            }
            b[i_loc] = sum; 
        }
    }
}

int read_array(FILE* fp, double* buf, int len) {
    for (int i = 0; i < len; ++i) {
        if (fscanf(fp, "%lf", buf + i) != 1) {
            return -1;
        }
    }

    return 0;
}

// читает 1 процесс и рассылает остальным, читает в буфер
// в остальных процессах есть буфер для получения
int read_matrix(double* a, int n, int m, int p, int k, 
    const char* name, double* buf /* буфер m x n на блочную строку */, MPI_Comm com) {

    int main_k = 0;
    FILE* fp = nullptr;
    int err = 0;
    if (k == main_k) {
        fp = fopen(name, "r");
        if (fp == nullptr) { err = 1; }
    }
    MPI_Bcast(&err, 1, MPI_INT, main_k, com);
    if (err) { return err; }
    // в процессе main_k файл открыт
    memset(buf, 0, n*m*sizeof(double));

    // число блочных строк 
    int max_b = (n + m - 1) / m;
    for (int b = 0; b < max_b; ++b) {
        // цикл по блочным строкам в каждом процессе идет
        // owner of row
        int owner = b % p;
        int num_rows = ((b + 1) * m  <= n ? m : n - b * m);
        // локальный номер блочной строки
        int b_loc = b / p;
        if (k == main_k) {
            err += read_array(fp, buf, n*num_rows);
            if (owner == k) {
                // копируем на место 
                memcpy(a + b_loc*m*n, buf, n*num_rows*sizeof(double)); // ТУТ ВРОДЕ ДИЧЬ НАПИСАНА ПРО B_LOC
            } else {
                // отправляем получателю
                MPI_Send(buf, n*num_rows, MPI_DOUBLE, owner, 0 /*tag*/, com);
            }
        } else if (k == owner) {
            // получает строку buf
            MPI_Status st;
            MPI_Recv(a + b_loc*m*n, n*num_rows, MPI_DOUBLE, main_k, 0 /*tag*/, com, &st);
        }
    }

    if (k == main_k) {
        fclose(fp); 
        fp = nullptr;
    }

    MPI_Bcast(&err, 1, MPI_INT, main_k, com);
    if (err) { return err; }

    return 0;
}

int print_array(double* a, int n, int m, int printed_rows, int max_print);
void print_matrix(double* a, int n, int m, int p, int k,
    double* buf, int max_print, MPI_Comm com) {

    int main_k = 0; /*только 0*/
    int max_b = (n + m - 1) / m;
    int printed_rows = 0;
    for (int b = 0; b < max_b; ++b) {
        // owner of row
        int owner = b % p;
        // num of rows
        int num_rows = (b + 1) * m <= n ? m : n - b * m; 
        // loc num of row
        int b_loc = b / p;
        if (k == main_k) {
            if (owner == main_k) {
                // строка у главного процесса
                printed_rows += print_array(a + b_loc*m*n, n, num_rows, printed_rows, max_print);
            } else {
                MPI_Status st;
                MPI_Recv(buf, n*num_rows, MPI_DOUBLE, owner, 0/*tag*/, com, &st);
                // вывести из буфера
                printed_rows += print_array(buf, n, num_rows, printed_rows, max_print);
            }
        } else if (owner == k) {
            // остальные 
            MPI_Send(a + b_loc*m*n, n*num_rows, MPI_DOUBLE, main_k, 0/*tag*/, com);
        }
    }
}

void print_b(double* a, int n, int m, int p, int k,
    double* buf, int max_print, MPI_Comm com) {

    int main_k = 0; /*только 0*/
    int max_b = (n + m - 1) / m;
    int printed_rows = 0;
    for (int b = 0; b < max_b; ++b) {
        // owner of row
        int owner = b % p;
        // num of rows
        int num_rows = (b + 1) * m <= n ? m : n - b * m; 
        // loc num of row
        int b_loc = b / p;
        if (k == main_k) {
            if (owner == main_k) {
                // строка у главного процесса
                printed_rows += print_array(a + b_loc*m, 1, num_rows, printed_rows, max_print);
            } else {
                MPI_Status st;
                MPI_Recv(buf, num_rows, MPI_DOUBLE, owner, 0/*tag*/, com, &st);
                // вывести из буфера
                printed_rows += print_array(buf, 1, num_rows, printed_rows, max_print);
            }
        } else if (owner == k) {
            // остальные 
            MPI_Send(a + b_loc*m, num_rows, MPI_DOUBLE, main_k, 0/*tag*/, com);
        }
    }
}

// печать прямоугольной матрицы m x n не более чем max_print
// при этом printed_rows уже выведены
int print_array(double* a, int n, int m, int printed_rows, int max_print) {
    if (printed_rows >= max_print) {
        return 0;
    }

    int p_n = (n > max_print ? max_print : n);
    int p_m = printed_rows + m < max_print ? m : max_print - printed_rows;
    for (int i = 0; i < p_m; ++i) {
        for (int j = 0; j < p_n; ++j) {
            printf(" %10.3e", a[i*n + j]);
        }
        printf("\n");
    }

    return p_m;
}

void solution(int argc, char* argv[], MPI_Comm com, int p, int k);
void get_block(int i, int j, int n, int m, int f, int l, double* matrix, double* block1);
void put_block(int i, int j, int n, int m, int k, int l, double* block, double* matrix);
void put_vector(int i, int m, int k, int l, double* b_i, double* b);
void swap_rows(double* matrix, int n, int i, int j);
bool inverse_matrix(int m, double* matrix, double* identity, double a_norm);
void matrix_product(int n, int m, int k, double* a, double* b, double* c);
void subtract_matrix_inplace(int n, int m, double* a, double* b);

int main(int argc, char* argv[]) {

    int p, k;
    MPI_Comm com = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(com, &p);
    MPI_Comm_rank(com, &k);

    solution(argc, argv, com, p, k); 
    MPI_Finalize();

    return 0;
}

void solution(int argc, char* argv[], MPI_Comm com, int p, int k) {
    int n = std::stoi(argv[1]);
    int m = std::stoi(argv[2]);
    int r = std::stoi(argv[3]);
    int s = std::stoi(argv[4]);

    int rows = get_rows(n, m, p, k);
    double* a = new double[rows*m*n];
    double* b = new double[rows*m];
    double* x = new double[n];
    double* buf = new double[m*n];
    double* buf_b = new double[m];
    double* block1 = new double[m*m];
    double* block2 = new double[m*m];
    double* block3 = new double[m*m];

    if (s) {
        init_matrix(a, /* своя часть */ n, m, p, k, s, &f);
    } else {
        const char* name_file = argv[argc - 1];
        read_matrix(a, n, m, p, k, name_file, buf /* буфер m x n на блочную строку */, com);
    }

    init_b(a, b, n, m, p, k);
    //print_matrix(a, n, m, p, k, buf, r, com);
    //print_b(b, n, m, p, k, buf, r, com);

    int f = n / m;
    int l = n - f * m;
    int h = l ? f + 1 : f;

    for (int i_glob_m = 0; i_glob_m < f; ++i_glob_m) {
        int i_loc_m = i_glob_m / p;
        int main_k = i_glob_m % p; 
        if (k == main_k) {
            get_block(i_loc_m, i_glob_m, n, m, f, l, a, block1);
            inverse_matrix(m, block1, block2, 1e-10); // НОРМУ МАТРИЦЫ ПОТОМ ПОСЧИТАТЬ И ЗАМЕНИТЬ АРГ.            
            for (int j_glob_m = i_glob_m; j_glob_m < f; ++j_glob_m) {
                get_block(i_loc_m, j_glob_m, n, m, f, l, a, block1); 
                matrix_product(m, m, m, block2, block1, block3); 
                put_block(i_loc_m, j_glob_m, n, m, f, l, block3, a);
            }
            
            if (l) {
                get_block(i_loc_m, f, n, m, f, l, a, block1);
                matrix_product(m, m, l, block2, block1, block3);
                put_block(i_loc_m, f, n, m, f, l, block3, a);
            }

            matrix_product(m, m, 1, block2, b + m * i_loc_m, block3); 
            put_vector(i_loc_m, m, f, l, block3, b);

            memcpy(buf, a + m*n*i_loc_m, n*m*sizeof(double)); 
            memcpy(buf_b, b + m*i_loc_m, m*sizeof(double));
        } 

        // насколько тут всё плохо и можно ли сделать лучше, чтобы не было двух обменов?
        MPI_Bcast(buf, n*m, MPI_DOUBLE, main_k, com);
        MPI_Bcast(buf_b, m, MPI_DOUBLE, main_k, com); 
        
        int offset = main_k <= k ? k - main_k : p + k - main_k; 
        int loc_start_m = (i_glob_m + offset) / p;
        if (k == main_k) { loc_start_m++; }

        for (int i = loc_start_m; i < rows; ++i) {
            get_block(i, i_glob_m, n, m, f, l, a, block1);
            int multiplier_rows = (m + l2g(n, m, p, k, i * m) <= n ? m : l);
            for (int j = 0; j < h; ++j) {
                int block_cols = j < f ? m : l;
                get_block(0, j, n, m, f, l, buf, block2); 
                matrix_product(multiplier_rows, m, block_cols, block1, block2, block3);
                get_block(i, j, n, m, f, l, a, block2);
                subtract_matrix_inplace(m, m, block2, block3);
                put_block(i, j, n, m, f, l, block2, a);
            }

            matrix_product(multiplier_rows, m, 1, block1, buf_b, block3);
            subtract_matrix_inplace(1, multiplier_rows, b + m*i, block3);
        }
    }

    // if l and обратный ход гаусса летс го делать.
    if (l && f % p == k) {
        get_block(f/p, f, n, m, f, l, a, block1);  
        inverse_matrix(l, block1, block2, 1e-10); // НОРМУ МАТРИЦЫ ПОТОМ ПОСЧИТАТЬ И ЗАМЕНИТЬ АРГ.
        matrix_product(l, l, 1, block2, b + m*(f/p), block3);
    }

    if (l) { // сюда зайдут либо все процессы, либо никакой, поэтому всё норм
        MPI_Bcast(block3, l, MPI_DOUBLE, f % p, com);
        put_vector(f, m, f, l, block3, x);
    }

    for (int i_glob_m = f - 1; i_glob_m >= 0; --i_glob_m) {
        memset(block2, 0, m*sizeof(double));
        int i_loc_m = i_glob_m / p;
        int owner = i_glob_m % p;
        if (k == owner) { 
            for (int j = i_glob_m + 1; j < h; ++j) {
                get_block(i_loc_m, j, n, m, f, l, a, block1);
                int cols = j < f ? m : l;
                matrix_product(m, cols, 1, block1, x + m*j, block3);

                for (int p = 0; p < m; ++p) {
                    block2[p] += block3[p];
                }                             
            }

            for (int vv = 0; vv < m; ++vv) {
                buf_b[vv] = b[m*i_loc_m + vv] - block2[vv]; 
            }
        }

        MPI_Bcast(buf_b, m, MPI_DOUBLE, owner, com);
        put_vector(i_glob_m, m, f, l, buf_b, x);
    }

    //print_matrix(a, n, m, p, k, buf, r, com);
    //print_b(b, n, m, p, k, buf, r, com);

    delete[] a;
    delete[] b;
    delete[] x;
    delete[] buf;
    delete[] buf_b;
    delete[] block1;
    delete[] block2;
    delete[] block3;
}

void get_block(int i, int j, int n, int m, int f, int l, double* matrix, double* block1) {
    int h = i < f ? m : l;
    int w = j < f ? m : l;
    
    int ind = 0;

    for (int p = 0; p < h; ++p) {
        for (int q = 0; q < w; ++q) {
            block1[ind] = matrix[n * (m * i + p) + m * j + q];
            ind++;
        }
    }
}

void put_block(int i, int j, int n, int m, int k, int l, double* block, double* matrix) {
    int h = i < k ? m : l;
    int w = j < k ? m : l;

    int ind = 0;
    for (int p = 0; p < h; ++p) {
        for (int q = 0; q < w; ++q) {
            matrix[n * (m * i + p) + m * j + q] = block[ind];
            ind++;
        }
    }
}

void put_vector(int i, int m, int k, int l, double* b_i, double* b) {
    int length = i < k ? m : l;
    for (int p = 0; p < length; ++p) {
        b[m*i + p] = b_i[p];
    }
}

void swap_rows(double* matrix, int n, int i, int j) {
    for (int k = 0; k < n; ++k) {
        std::swap(matrix[n*i + k], matrix[n*j + k]);
    }
}

void subtract_matrix_inplace(int n, int m, double* a, double* b) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            a[m*i + j] -= b[m*i + j];           
        }
    }
}

bool inverse_matrix(int m, double* matrix, double* identity, double a_norm) {

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            identity[i * m + j] = (i != j) ? 0 : 1;
        }
    }

    for (int i = 0; i < m; ++i) {
        double max_elem = matrix[i * m + i];
        int row_max_elem = i;
        for (int j = i + 1; j < m; ++j) {
            if (std::fabs(matrix[j * m + i]) > std::fabs(max_elem)) {
                max_elem = matrix[j * m + i];
                row_max_elem = j;
            }
        }

        swap_rows(matrix, m, i, row_max_elem);
        swap_rows(identity, m, i, row_max_elem);

        if (std::fabs(max_elem) < EPS * a_norm) {
            return false;    
        }

        double factor = 1 / max_elem;
        for (int s = 0; s < i; ++s) {
            identity[i * m + s] *= factor;
        }

        for (int s = i; s < m; ++s) {
            matrix[i * m + s] *= factor;
            identity[i * m + s] *= factor;
        }

        for (int k = i + 1; k < m; ++k) {
            double multiplier = -matrix[k * m + i];
            for (int p = 0; p < i + 1; ++p) {
                identity[k * m + p] += identity[i * m + p] * multiplier;
            }

            for (int p = i + 1; p < m; ++p) { 
                matrix[k * m + p] += matrix[i * m + p] * multiplier;
                identity[k * m + p] += identity[i * m + p] * multiplier;
            }
        }
    }

    for (int i = m - 1; i > 0; --i) {
        for (int k = i - 1; k >= 0; --k) {
            double multiplier = -matrix[k * m + i];
            for (int p = 0; p < m; ++p) { 
                identity[k * m + p] += identity[i * m + p] * multiplier;
            }
        }
    }

    return true;
}

void matrix_product(int n, int m, int k, double* a, double* b, double* c) {
    
    for (int i = 0; i < n*k; ++i) {
        c[i] = 0;
    }
    
    double sum00, sum01, sum02, sum10, sum11, sum12, sum20, sum21, sum22;
    
    int v3 = n%3;
    int h3 = k%3;
    
    for (int i = 0; i < v3; ++i) {
        for (int j = 0; j < h3; ++j) {
            sum00 = 0;
            for (int p = 0; p < m; ++p) {
                sum00 += a[m*i + p] * b[k*p + j];
            }
            
            c[k*i + j] = sum00;
        }
        for (int j = h3; j < k; j+=3) {
            sum00 = 0; sum01 = 0; sum02 = 0;
            for (int p = 0; p < m; ++p) {
                double factor = a[m*i + p];
                sum00 += factor * b[k*p + j];
                sum01 += factor * b[k*p + j + 1];
                sum02 += factor * b[k*p + j + 2];
            }
            c[k*i + j] = sum00;
            c[k*i + j + 1] = sum01;
            c[k*i + j + 2] = sum02;
        }
    }
    
    for (int i = v3; i < n; i+=3) {   
        for (int j = 0; j < h3; ++j) {
            sum00 = 0; sum01 = 0; sum02 = 0;
            for (int p = 0; p < m; ++p) {
                double factor = b[k*p + j];
                sum00 += a[m*i + p] * factor;
                sum01 += a[m*(i + 1) + p] * factor;
                sum02 += a[m*(i + 2) + p] * factor;
            }
            c[k*i + j] = sum00;
            c[k*(i + 1) + j] = sum01;
            c[k*(i + 2) + j] = sum02;
        }
        
        for (int j = h3; j < k; j+=3) {
            sum00 = 0; sum01 = 0; sum02 = 0; sum10 = 0; sum11 = 0; sum12 = 0;
            sum20 = 0; sum21 = 0; sum22 = 0;
            for (int p = 0; p < m; ++p) {
                sum00 += a[m*i + p] * b[k*p + j];
                sum01 += a[m*i + p] * b[k*p + j + 1];
                sum02 += a[m*i + p] * b[k*p + j + 2];
                sum10 += a[m*(i + 1) + p] * b[k*p + j];
                sum11 += a[m*(i + 1) + p] * b[k*p + j + 1];
                sum12 += a[m*(i + 1) + p] * b[k*p + j + 2];
                sum20 += a[m*(i + 2) + p] * b[k*p + j];
                sum21 += a[m*(i + 2) + p] * b[k*p + j + 1];
                sum22 += a[m*(i + 2) + p] * b[k*p + j + 2];
            }
            c[k*i + j] = sum00;
            c[k*i + j + 1] = sum01; 
            c[k*i + j + 2] = sum02;
            c[k*(i + 1) + j] = sum10;
            c[k*(i + 1)  + j + 1] = sum11;
            c[k*(i + 1) + j + 2] = sum12;
            c[k*(i + 2) + j] = sum20;
            c[k*(i + 2) + j + 1] = sum21;
            c[k*(i + 2) + j + 2] = sum22;
        }
    }    
}
