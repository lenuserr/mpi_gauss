#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <string>
#include <cmath>
// а за ноут сел так-то вообще в 11 00.
// 11:45. Переписываю, чтобы разделение стало по столбцам.
// Весь день есть кодить, максимум Димон наберет и усё, ваще не проблема.
// А самые оптимальные звонки теперь это в поездках в автобусе, реально, супер оптимизация.
// Ну и по пути в ашан после пар бгч например.
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
int get_max_cols(int n, int m, int p) {
    int b = (n + m - 1) / m;
    return (b + p - 1) / p;
}

// число столбцов на процесс k
int get_cols(int n, int m, int p, int k) {
    // число блочных столбцов
    int b = (n + m - 1) / m;
    return (b % p) >= k + 1 ? (b / p) + 1 : b / p; 
}

// в каком процессе лежит столбец i_glob
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
    int cols = get_cols(n, m, p, k); // число блочных столбцов на процесс k
    for (int i_loc_m = 0; i_loc_m < cols; ++i_loc_m) {
        for (int j = 0; j < n; ++j) {
            for (int i_loc = i_loc_m * m; i_loc < (i_loc_m + 1) * m; ++i_loc) {
                int i_glob = l2g(n, m, p, k, i_loc);
                if (i_glob >= n) {
                    break; 
                }
                
                a[cols*m*j + i_loc] = f(s, n, j, i_glob);
            }
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
    const char* name, double* buf /* буфер n x m на блочный столбец */, MPI_Comm com) {
    
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

    // число блочных столбцов
    int max_b = (n + m - 1) / m;
    for (int b = 0; b < max_b; ++b) {
        // цикл по блочным столбцам в каждом процессе идет
        // owner of col
        int owner = b % p;
        int num_rows = ((b + 1) * m  <= n ? m : n - b * m);
        // локальный номер блочного столбца
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

int print_array(double* a, int n, int m, int printed_rows, int max_print);
void print_matrix(double* a, int n, int m, int p, int k,
    double* buf, int max_print, MPI_Comm com) {
    
    int main_k = 0; /*только 0*/
    int max_b = (n + m - 1) / m;
    int printed_rows = 0;
    for (int b = 0; b < max_b; ++b) {
        // owner of col
        int owner = b % p;
        // num of cols
        int num_rows = (b + 1) * m <= n ? m : n - b * m; 
        // loc num of col
        int b_loc = b / p;
        if (k == main_k) {
            if (owner == main_k) {
                // столбец у главного процесса
                printed_rows += print_array(a + b_loc*m*n, n, num_rows, printed_rows, max_print);
            } else {
                MPI_Status st;
                MPI_Recv(buf, n*num_rows, MPI_DOUBLE, owner, 0/*tag*/, com, &st);
                // вывести из буфера
                printed_rows += print_array(buf, n, num_rows, printed_rows, max_print);
            }
        } else if (owner != main_k) {
            // остальные 
            MPI_Send(a + b_loc*m*n, n*num_rows, MPI_DOUBLE, main_k, 0/*tag*/, com);
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

void initialization(int argc, char* argv[], MPI_Comm com, int p, int k);

int main(int argc, char* argv[]) {
    
    int p, k;
    MPI_Comm com = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(com, &p);
    MPI_Comm_rank(com, &k);
    
    initialization(argc, argv, com, p, k); 
    MPI_Finalize();
    return 0;
}

void initialization(int argc, char* argv[], MPI_Comm com, int p, int k) {
    int n = std::stoi(argv[1]);
    int m = std::stoi(argv[2]);
    int r = std::stoi(argv[3]);
    int s = std::stoi(argv[4]);
    
    int cols = get_cols(n, m, p, k);
    int l = n % m;
    double* a = new double[n*cols*m];
    double* buf = new double[n*m];
    
    if (s) {
        init_matrix(a, /* своя часть */ n, m, p, k, s, &f);
    } else {
        const char* name_file = argv[argc - 1];
        read_matrix(a, n, m, p, k, name_file, buf /* буфер m x n на блочную строку */, com);
    }
    

    delete[] a;
    delete[] buf;
}
