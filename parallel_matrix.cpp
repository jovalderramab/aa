#include "mpi.h"
#include <iostream>
#include <stdio.h>

void printdata(double d[], int N, int Nlocal);

int main(int argc, char **argv) {
  int dest, noProcesses, pid, tag, Nlocal;
  MPI_Status status;
  double sumtime = 0;
  const int N = std::atoi(argv[1]);
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &noProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  Nlocal = N / noProcesses;
  double *data = new double[N * Nlocal]{0.0};
  for (int ilocal = 0; ilocal < Nlocal; ++ilocal) {
    int jlocal = Nlocal * pid + ilocal;
    data[ilocal * N + jlocal] = 1.0;
  }
  tag = 0;
  if (0 == pid) {
    double tstart = MPI_Wtime();

    printdata(data, N, Nlocal);

    for (int ipid = 1; ipid < noProcesses; ++ipid) {
      MPI_Recv(&data[0], N * Nlocal, MPI_DOUBLE, ipid, tag, MPI_COMM_WORLD,
               &status);

      double tend = MPI_Wtime();
      double total_time = tend - tstart;
      sumtime = total_time + sumtime;

      printdata(data, N, Nlocal);
    }
    std::cout << "Ancho de banda:"
              << "\t"
              << N * Nlocal * sizeof(double) / sumtime / 1.0e6 / noProcesses
              << std::endl;
  } else {
    dest = 0;
    MPI_Send(&data[0], N * Nlocal, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
  }
  MPI_Finalize();
  delete[] data;
  return 0;
}

void printdata(double d[], int N, int Nlocal) {
  for (int ii = 0; ii < Nlocal; ii++) {
    for (int jj = 0; jj < N; jj++) {
      std::cout << d[jj + N * ii] << "\t";
    }
    std::cout << std::endl;
  }
}
