#include "mpi.h"
#include <iostream>
#include <stdio.h>
#include <vector>

void matrixId(int N, int pid, int np);

int main(int argc, char **argv) {

  /* MPI Variables */
  int np, pid;

  /* MPI setup */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  // get the size of matrix
  int N = atoi(argv[1]);

  /*matrixId return the complete identity matrix
    and the average bandwith of the comumunication with each node*/

  matrixId(N, pid, np);

  /* finish */
  MPI_Finalize();

  return 0;
}

void matrixId(int N, int pid, int np) {

  /* MPI Variables */

  int dest, src, tag;
  MPI_Status status;

  /* Problem variables */

  int Nrows = N / np;
  std::vector<double>arr(Nrows * N);
  std::fill(arr.begin(),arr.end(),0.0);
  //double *arr = new double[Nrows * N]{0};
  double lower, upper;
  double time_node = 0;

  lower = pid * Nrows;
  upper = (pid + 1) * Nrows;

  /* Fill matrix for subproblem */

  for (int ilocal = 0; ilocal < Nrows; ++ilocal) {
    arr[(ilocal * N) + Nrows * pid + ilocal] = 1.0;
  }
  // for (int i = 0; i < Nrows; ++i) {
  // for (int j = 0; j < N; j++) {

  // if (j == pid * Nrows + i) {
  // arr[(i * N) + j] = 1;
  //}
  //}

  /* Collect info and print results */

  tag = 0;

  if (0 == pid) { /* Master */

    // print its part

    for (int i = 0; i < Nrows; ++i) {
      for (int j = 0; j < N; ++j) {
        std::cout << arr[(i * N) + j] << " ";
      }
      std::cout << "\n";
    }

    // Print the other parts

    for (src = 1; src < np; ++src) {

      double tstart = MPI_Wtime();

      MPI_Recv(&arr[0], Nrows * N, MPI_DOUBLE, src, tag, MPI_COMM_WORLD,
               &status);

      double tend = MPI_Wtime();

      // accumulates communication times with each node
      time_node = time_node + (tend - tstart);

      for (int i = 0; i < Nrows; ++i) {
        for (int j = 0; j < N; ++j) {
          std::cout << arr[(i * N) + j] << " ";
        }
        std::cout << "\n";
      }
    }

    std::cout << "The average bandwith is:"
              << "\t" << Nrows * N * sizeof(double) / (time_node / (np - 1))
              << std::endl;

  }

  else { /* slaves only send */
    dest = 0;
    MPI_Send(&arr[0], Nrows * N, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
  }
}
