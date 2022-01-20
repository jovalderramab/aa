#include "mpi.h"
#include <iostream>
#include <stdlib.h>
#include <vector>

void print( const std::vector<int>& matrix, const int& columns, const int& rows);

int main(int argc, char **argv)
{
    int pId;    // rank of process
    int nP;     // number of processes
    int N = atoi(argv[1]);  // size of matrix

    MPI_Init(&argc, &argv);     /* Mandatory */
    MPI_Comm_size(MPI_COMM_WORLD, &nP);
    MPI_Comm_rank(MPI_COMM_WORLD, &pId);

    int rows = N/nP;    // rows of local matrix
    int Lower = pId*rows;
    int Upper = (pId+1)*rows;

    std::vector<int> Local(N*rows,0);   // local matrix (1D)
    for (int ii=0; ii<rows; ii++){      // fill local matrix
        Local[N*ii+(Lower+ii)] = 1;
    }

    int tag = 0;

    if (pId == 0){
        std::vector<double> bandwidth(nP-1,0);
        print(Local, N, rows);
        for (int src=1; src < nP; src++){
            double tstart = MPI_Wtime();
            MPI_Recv(&Local[0], N*rows, MPI_INT, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            double tend = MPI_Wtime();
            bandwidth[src-1]= tend - tstart;     // time period spent by Recv
            print(Local, N, rows);      // print local matrix
        }

        std::cout << "\nBandwidth" << std::endl;
        double average{0};
        for (int ii=0; ii<nP-1; ii++){
            std::cout << "Process " << ii+1 << " :\t" << N*rows*sizeof(int)/(bandwidth[ii])/1.0e6 << "\tMB/s."<< std::endl;
            average += 1/bandwidth[ii];
        }
        std::cout << "Average : \t" << N*rows*sizeof(int)*average/(1.0e6*(nP-1)) << "\tMB/s."<< std::endl;

    } else{        
        MPI_Send(&Local[0], N*rows, MPI_INT, 0, tag, MPI_COMM_WORLD);
    }

    MPI_Finalize();     /* Mandatory */

    return 0;
}

void print( const std::vector<int>& matrix, const int& columns, const int& rows){
    int total=columns*rows;
    for (int ii=0; ii<rows; ii++){
        for (int jj=0; jj<columns; jj++){
            std::cout << matrix[columns*ii+jj] << "\t";
        }
        std::cout << std::endl;
    }
}