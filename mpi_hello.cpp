#include "mpi.h"
#include <stdio.h>
#include <string.h>

void my_main(int argc, char* argv[], MPI_Comm com, int p, int k);

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
    
    int p, k;
    MPI_Comm com = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(com, &p);
    MPI_Comm_rank(com, &k);
    
    my_main(argc, argv, com, p, k); 
    MPI_Finalize();
    
    return 0;
}

void my_main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[], MPI_Comm com, int p, int k) {
    
    const int len = 1234;
    char buf[len];
    snprintf(buf, len, "Hello from process %d\n", k);
    
    int tag = 0;
    if (k == 0) {
        printf("%s", buf);
        [[maybe_unused]] int i;
        MPI_Status st;
        for (k = 1; k < p; ++k) {
            MPI_Recv(buf, len, MPI_CHAR, MPI_ANY_SOURCE, tag, com, &st);
            printf("%s", buf);
        }
    } else {
        MPI_Send(buf, strlen(buf) + 1, MPI_CHAR, 0, tag, com);
    }
}
