#include "mpi.h"
#include <iostream>
#include <cstdlib>
#include <vector>

void fill(std::vector<int>& Teil, int N, int Nrows, int pid);
void print(std::vector<int>& Teil, int N, int Nrows);
double mean(std::vector<double>& Vec, int size);

int main(int argc, char **argv)
{
  /* MPI Variables */
  int np, pid;

  /* MPI setup */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  MPI_Status status;
  int tag = 0;

  /*create matrix*/
  const int N = std::atoi(argv[1]);
  int Nrows = N/np;
  std::vector<int> Teil(N*Nrows,0);

  fill(Teil, N, Nrows, pid);
  
  /*print matrix*/
  if (pid == 0){
    std::vector<double> Time(np-1,0.0);
    double tstart = 0;
    double tend = 0;
    double total_time = 0;
    print(Teil, N, Nrows);
    for (int kk = 1; kk < np; kk++){
      
      tstart = MPI_Wtime();
      MPI_Recv(&Teil[0], N*Nrows, MPI_INT, kk, tag, MPI_COMM_WORLD, &status);
      tend = MPI_Wtime();

      Time[kk-1] = tend - tstart;
     
      print(Teil, N, Nrows);
    }
    total_time = mean(Time,np-1);
    std::cout<<"Local size [bytes]: "<< N*Nrows*sizeof(int)<<"\t Ancho de banda [MB/s]: "<< N*Nrows*sizeof(int)/total_time/1.0e6<<"\n";
  }
  else{
    MPI_Send(&Teil[0], N*Nrows, MPI_INT, 0, tag, MPI_COMM_WORLD);
  }
  /* finish */
  MPI_Finalize();

  return 0;
}

void fill(std::vector<int>& Teil, int N, int Nrows, int pid){
  int wo = 0;
  for(int ii = 0; ii < Nrows; ii++){
    wo = (Nrows*pid)+(N+1)*ii;
    Teil[wo] = 1;  
  }
}

void print(std::vector<int>& Teil, int N, int Nrows){
  for (int ii = 0; ii < Nrows; ii++){
    for (int jj = 0; jj < N; jj++){
      std::cout<<Teil[N*ii+jj]<<" ";
     }
     std::cout<<"\n";
  }
}
double mean(std::vector<double>& Vec, int size){
  double Suma = 0;
  for (int ii = 0; ii < size; ii++){
    Suma+=Vec[ii];
  }
  return Suma/size;
}
