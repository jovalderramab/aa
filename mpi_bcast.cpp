#include "mpi.h"
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <numeric>

void send_data_collective(int size, int pid, int np);
void broadcast(int size, int pid , int np);

int main(int argc, char **argv)
{
  int np, pid;

  /* MPI setup */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  const int SIZE = std::atoi(argv[1]);
  send_data_collective(SIZE, pid, np);

  MPI_Finalize();

  return 0;
}

void send_data_collective(int size, int pid, int np)
{
  
  double * data = new double [size];
 
  
  double start = MPI_Wtime();
  MPI_Bcast(&data[0], size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  double mid = MPI_Wtime();
  broadcast(size, pid ,  np);
  double end = MPI_Wtime();

  if (0 == pid) {
    int datasize = sizeof(double)*size;
    std::cout << size << "\t" << size*datasize/(mid-start)/1.0e6<<"\t"<< size*datasize/(end-mid)/1.0e6 << "\n";
  }
  delete [] data;
}

void broadcast(int size, int pid , int np)
{
  double * data= new double [size];
  if (0 ==pid ) {
    
    for (int i = 1; i < np; i++) {
        MPI_Send(&data[0], size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    }
  } else {
    	MPI_Recv(&data[0], size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  delete [] data;
}
