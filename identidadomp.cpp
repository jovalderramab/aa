#include "mpi.h"
#include <iostream>
#include <cstdlib>

void identityomp(int N, int pid, int np);
void printmat(int N, int cols, int *a);
void filldiag(int N, int pid, int cols, int *a);

int main(int argc, char **argv)
{
  /* MPI Variables */
  int np, pid;

  /* MPI setup */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  int N= std::atoi(argv[1]);
   
  identityomp( N, pid, np);

  /* finish */
  MPI_Finalize();

  return 0;
}

void identityomp(int N, int pid, int np)
{
  MPI_Status status;
  int tag = 0, cols=N/np;
  int *a= new int[N*cols]{0};
  double start, end=0;
if(np>1){  
    if (0 == pid) { /* Master */
      int *buff= new int[N*cols]{0}; 
      double *anchos= new double[np-1]{0};
      filldiag( N, pid,cols,a);
      printmat( N,cols,a);
      double w=0;
      for (int kk=1;kk<np;++kk){
	     start=MPI_Wtime();
	       
	      MPI_Recv(buff, N*cols, MPI_INT, kk, tag, MPI_COMM_WORLD, &status);
		
		end=MPI_Wtime();
		
		w=(N*cols*sizeof(double))/((end-start)*1.0e6);
		anchos[kk-1]=w;
 	     printmat( N,cols,buff);
	}
     	std::cout.precision(7); 
      for(int ss=0;ss<np-1;++ss){
	      std::cout<<"El ancho de banda del proceso "<<ss+1<<" al proceso 0, en megabytes por segundo, es: "<<std::fixed<<anchos[ss]<<"\n";
      }
      double prom=0;
      for(int rr=0;rr<np-1;rr++)prom+=anchos[rr];

      std::cout<<"El ancho de banda promedio, en megabytes por segundo, es: "<<std::fixed<<prom/(np-1)<<"\n";
      delete []a;
    }
    else { /* slaves only send */
      
      filldiag( N,pid,cols,a);
      
      
      MPI_Send(a, N*cols, MPI_INT, 0, tag, MPI_COMM_WORLD);

      delete []a;
    }
  }
	else{
		filldiag(N,0,N,a);
		printmat(N,N,a);
		std::cout<<"Solo hay un proceso, por ende no hay anchos de banda de comunicaciones"<<"\n";	
	}
}

void printmat(int N, int cols, int *a){
  for(int ii=0;ii<cols;++ii){
	for(int jj=0;jj<N;++jj){
	  std::cout<<a[ii*N+jj]<<" ";
	}
	std::cout<<"\n";
      }


}

void filldiag(int N, int pid, int cols, int *a){

  for(int ii=0; ii<cols;++ii) a[pid*cols+ii*(N+1)]=1;


}

