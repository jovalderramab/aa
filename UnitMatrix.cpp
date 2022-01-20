#include "mpi.h"
#include <cstdio>
#include <cstdlib>




int main(int argc, char **argv){
  
  
  const int N= std::atoi(argv[1]); //Size matrix

  /* MPI variables */

  int pid = 0, np = 0;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  
  int i,k;
  int Mat[N*N/np]={0}; //size array for each process


  int Nlocal = N/np;

  /*fill matrix  */
  for (i = 0; i < N*Nlocal; ++i) {
    k=pid*Nlocal;
    Mat[k]=1;
    for(int j=1;j<=Nlocal;++j){
      Mat[k+j*(N+1)]=1;}
  }

 
  int tag = 0;
  if (0 == pid) {
    for(int ll=0;ll<Nlocal;++ll){
    for(int ii=0;ii<N;++ii){
	printf(" %d  ",Mat[ii+N*ll]);}
        printf(" \n");}
    
    //double tstart=MPI_Wtime();
    for (int src = 1; src < np; ++src) {
      double tstart = MPI_Wtime();
      MPI_Recv(&Mat, N*Nlocal, MPI_INT, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      double tend = MPI_Wtime();
      double time = tend -tstart; //Time send-recv
      //Print matrix
      for(int ll=0;ll<Nlocal;++ll){
          for(int ii=0;ii<N;++ii){
	    printf(" %d  ",Mat[ii+N*ll]);}
	  printf(" \n");}
      printf("Para %d el ancho de banda %lf ",src,sizeof(double)/(time*1000000));printf(" \n");
      }
     } else {
    //send matrix to pid 0;
    int dest = 0;
    MPI_Send(&Mat,N*Nlocal, MPI_INT, dest, tag, MPI_COMM_WORLD);
    }
  
  MPI_Finalize();
  return 0;
}
