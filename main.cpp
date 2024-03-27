#define EPS 1e-14

#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>

#include <string>
#include <cmath>
#include "main.h"

double get_cpu_time() {
    struct rusage buf;
    getrusage(RUSAGE_THREAD, &buf);
    return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec * 1e-6;
}

int l2g(int m, int p, int k, int i_loc) {
    int i_loc_m = i_loc / m;
    int i_glob_m = i_loc_m * p + k;
    return i_glob_m * m + i_loc % m;
}

// global to local, если разделение по строкам
int g2l(int m, int p, int i_glob) {
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
int get_k(int m, int p, int i_glob) {
    int i_glob_m = i_glob / m;
    return i_glob_m % p;
}

void init_matrix(double* a, int n, int m, int p, int k, int s, double (*func)(int, int, int, int)) {
    int rows = get_rows(n, m, p, k); // число блочных строк в процессе k

    for (int i_loc_m = 0; i_loc_m < rows; ++i_loc_m) {
        for (int i_loc = i_loc_m * m; i_loc < (i_loc_m + 1) * m; ++i_loc) {
            int i_glob = l2g(m, p, k, i_loc);
            if (i_glob >= n) {
                break;
            }
            for (int j = 0; j < n; ++j) {
                a[n*i_loc + j] = func(s, n, i_glob, j);
            }
        }
    }
}

void init_b(double* a, double* b, int n, int m, int p, int k) {
    int rows = get_rows(n, m, p, k); // число блочных строк в процессе k

    for (int i_loc_m = 0; i_loc_m < rows; ++i_loc_m) {
        for (int i_loc = i_loc_m * m; i_loc < (i_loc_m + 1) * m; ++i_loc) {
            int i_glob = l2g(m, p, k, i_loc);
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

void subtract_matrix_inplace(int n, int m, double* a, double* b) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            a[m*i + j] -= b[m*i + j];           
        }
    }
}

void swap_rows(double* matrix, int n, int i, int j) {
    for (int k = 0; k < n; ++k) {
        std::swap(matrix[n*i + k], matrix[n*j + k]);
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

double matrix_norm(double* a, /* своя часть */ int n, int m, int p, int k, MPI_Comm com) {
    int rows = get_rows(n, m, p, k);
    
    double norm = -1;
    double global_sum;
    double sum;
    for (int j = 0; j < n; ++j) {
        sum = 0;
        for (int i_loc_m = 0; i_loc_m < rows; ++i_loc_m) {
            for (int i_loc = i_loc_m * m; i_loc < (i_loc_m + 1) * m; ++i_loc) {
                int i_glob = l2g(m, p, k, i_loc);
                if (i_glob >= n) {
                    break;
                }

                sum += std::fabs(a[n*i_loc + j]);
            }
        }

        MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, com);
        norm = std::max(norm, global_sum);
    }

    return norm;
}

double func(int s, int n, int i, int j) { 
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
                memcpy(a + b_loc*m*n, buf, n*num_rows*sizeof(double)); 
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

double r2_eval(int n, double* x) {   
    double r2 = 0;
    for (int i = 1; i <= n; ++i) {
        r2 += std::fabs(x[i - 1] - (i % 2));
    }

    return r2;
}

void free_memory(double* a, double* b, double* x, double* buf, double* buf_b,
    double* block1, double* block2, double* block3) {

    delete[] a;
    delete[] b;
    delete[] x;
    delete[] buf;
    delete[] buf_b;
    delete[] block1;
    delete[] block2;
    delete[] block3;
}

void solution(int argc, char* argv[], MPI_Comm com, int p, int k);
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
    const int task = 9;

    int n = std::stoi(argv[1]);
    int m = std::stoi(argv[2]);
    int r = std::stoi(argv[3]);
    int s = std::stoi(argv[4]);

    int rows = get_rows(n, m, p, k);
    double* a = new double[rows*m*n];
    double* buf = new double[m*n + m];

    if (s) {
        init_matrix(a, /* своя часть */ n, m, p, k, s, &func);
    } else {
        const char* name_file = argv[argc - 1];
        if (read_matrix(a, n, m, p, k, name_file, buf /* буфер m x n на блочную строку */, com)) {
            
            if (k == 0) {
                printf (
                "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
                argv[0], task, -1., -1., 0., 0., s, n, m, p);
            }

            delete[] a;
            delete[] buf;
            return;
        }
    }

    double* b = new double[rows*m];
    double* x = new double[n];
    double* buf_b = new double[m];
    double* block1 = new double[m*m];
    double* block2 = new double[m*m];
    double* block3 = new double[m*m];

    init_b(a, b, n, m, p, k);

    double a_norm = matrix_norm(a, n, m, p, k, com);
    double t1 = get_cpu_time();

    int f = n / m;
    int l = n - f * m;
    int h = l ? f + 1 : f;

    int err = 0;
    for (int i_glob_m = 0; i_glob_m < f; ++i_glob_m) {
        int i_loc_m = i_glob_m / p;
        int main_k = i_glob_m % p; 
        if (k == main_k) {
            get_block(i_loc_m, i_glob_m, n, m, f, l, a, block1);
            if (!inverse_matrix(m, block1, block2, a_norm)) { err = -1; }             
            for (int j_glob_m = i_glob_m + 1; j_glob_m < f; ++j_glob_m) {
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
            memcpy(buf + n*m, b + m*i_loc_m, m*sizeof(double)); 
        } 

        MPI_Bcast(buf, n*m + m, MPI_DOUBLE, main_k, com); 
        MPI_Bcast(&err, 1, MPI_INT, main_k, com);

        if (err) {
            if (k == 0) {
                printf (
                "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
                argv[0], task, -1., -1., t1, 0., s, n, m, p);
            }

            free_memory(a, b, x, buf, buf_b, block1, block2, block3);
            return;
        }
        
        int offset = main_k <= k ? k - main_k : p + k - main_k; 
        int loc_start_m = (i_glob_m + offset) / p;
        if (k == main_k) { loc_start_m++; }

        for (int i = loc_start_m; i < rows; ++i) {
            get_block(i, i_glob_m, n, m, f, l, a, block1);
            int multiplier_rows = (m + l2g(m, p, k, i * m) <= n ? m : l);
            for (int j = i_glob_m + 1; j < h; ++j) {
                int block_cols = j < f ? m : l;
                get_block(0, j, n, m, f, l, buf, block2); 
                matrix_product(multiplier_rows, m, block_cols, block1, block2, block3);
                get_block(i, j, n, m, f, l, a, block2);
                subtract_matrix_inplace(multiplier_rows, block_cols, block2, block3);
                put_block(i, j, n, m, f, l, block2, a);
            }

            matrix_product(multiplier_rows, m, 1, block1, buf + n*m, block3);
            subtract_matrix_inplace(1, multiplier_rows, b + m*i, block3);
        }
    }

    int var = f/p; 
    if (l && f == k + var*p) {
        get_block(var, f, n, m, f, l, a, block1);  
        if (!inverse_matrix(l, block1, block2, a_norm)) { err = -1; } 
        matrix_product(l, l, 1, block2, b + m*var, block3);
    }

    MPI_Bcast(&err, 1, MPI_INT, f - var*p, com);

    if (err) {
        if (k == 0) {
            printf (
            "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
            argv[0], task, -1., -1., t1, 0., s, n, m, p);
        }

        free_memory(a, b, x, buf, buf_b, block1, block2, block3);
        return;
    }    

    if (l) { // сюда зайдут либо все процессы, либо никакой, поэтому всё норм
        MPI_Bcast(block3, l, MPI_DOUBLE, f - var*p, com);
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

    t1 = get_cpu_time() - t1;

    // снова инициализируем a и b, чтобы посчитать r1 невязку.
    if (s) {
        init_matrix(a, /* своя часть */ n, m, p, k, s, &func);
    } else {
        const char* name_file = argv[argc - 1];
        if (read_matrix(a, n, m, p, k, name_file, buf /* буфер m x n на блочную строку */, com)) {
            
            if (k == 0) {
                printf (
                "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
                argv[0], task, -1., -1., t1, 0., s, n, m, p);
            }

            free_memory(a, b, x, buf, buf_b, block1, block2, block3);
            return;
        }
    }

    init_b(a, b, n, m, p, k);
    
    double t2 = get_cpu_time();

    double sum = 0;
    double global_sum;
    
    for (int i_loc_m = 0; i_loc_m < rows; ++i_loc_m) {
        int max_i_loc = ((i_loc_m + 1) * m) - 1;
        int multiplier_rows = l2g(m, p, k, max_i_loc) < n ? m : l;
        memset(block3, 0, multiplier_rows*sizeof(double));
        for (int j = 0; j < h; ++j) { 
            get_block(i_loc_m, j, n, m, f, l, a, block1);
            int cols = j < f ? m : l;
            matrix_product(multiplier_rows, cols, 1, block1, x + m*j, block2);
                    
            for (int u = 0; u < multiplier_rows; ++u) {
                block3[u] += block2[u];
            }          
        }

        for (int u = 0; u < multiplier_rows; ++u) {
            sum += std::fabs(block3[u] - b[m*i_loc_m + u]);
        }  
    }

    MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, com);
    double b_norm = matrix_norm(b, 1, m, p, k, com);
    double r1 = global_sum / b_norm;

    if (k == 0) { printf("A:\n"); }
    print_matrix(a, n, m, p, k, buf, r, com);
    
    if (k == 0) { 
        double r2 = r2_eval(n, x);
        t2 = get_cpu_time() - t2;
        printf("\nx:\n");
        for (int i = 0; i < std::min(r, n); ++i) {
            printf(" %10.3e", x[i]);
        }
        printf("\n\n");
        printf (
        "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
        argv[0], task, r1, r2, t1, t2, s, n, m, p);
    }

    free_memory(a, b, x, buf, buf_b, block1, block2, block3);
}
