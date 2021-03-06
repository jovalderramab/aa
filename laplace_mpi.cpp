#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "mpi.h"

// constants
// const double DELTA = 0.05;
// const double L = 1.479;
// const int N = int(L/DELTA)+1;
// const int STEPS = 200;

typedef std::vector<double> Matrix; // alias

void initial_conditions(Matrix & m, int nrows, int ncols);
void boundary_conditions(Matrix & m, int nrows, int ncols);
void evolve(Matrix & m, int nrows, int ncols);
void print(const Matrix & m, double delta, int nrows, int ncols);
void init_gnuplot(void);
void plot_gnuplot(const Matrix & m, double delta, int nrows, int ncols);

// parallel versions
void initial_conditions(Matrix & m, int nrows, int ncols, int pid, int np);
void boundary_conditions(Matrix & m, int nrows, int ncols, int pid, int np);
void print_screen(const Matrix & m, int nrows, int ncols, int pid, int np);
void print_slab(const Matrix & m, int nrows, int ncols);

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  int pid, np;
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  const int N = std::atoi(argv[1]);
  const double L = std::atof(argv[2]);
  const int STEPS = std::atoi(argv[3]);
  const double DELTA = L/N;

  // problem partition
  int NCOLS = N, NROWS = N/np +2; // include ghosts
  Matrix data(NROWS*NCOLS); // include ghosts cells
  initial_conditions(data, NROWS, NCOLS, pid, np);
  if (0 == pid) {std::cout << " After initial conditions ...\n";}
  print_screen(data, NROWS, NCOLS, pid, np); // todo
  boundary_conditions(data, NROWS, NCOLS, pid, np); // todo
  if (0 == pid) {std::cout << " After boundary conditions ...\n";}
  print_screen(data, NROWS, NCOLS, pid, np); // todo

  /*
  // Serial version
  Matrix data(N*N);
  initial_conditions(data, N, N, ...);
  print_screen(...);
  boundary_conditions(data, N, N, ...);
  init_gnuplot();
  for (int istep = 0; istep < STEPS; ++istep) {
    evolve(data, N, N);
    plot_gnuplot(data, DELTA, N, N);
  }
  */

  MPI_Finalize();
  return 0;
}
/////////////////////////////////////////////////////
// Parallel implementations
void initial_conditions(Matrix & m, int nrows, int ncols, int pid, int np)
{
  // same task for all pids, but fill with the pids to distinguish among thems
  for(int ii=0; ii<nrows; ++ii) {
    for(int jj=0; jj<ncols; ++jj) {
      m[ii*ncols + jj] = pid;
    }
  }
}

void boundary_conditions(Matrix & m, int nrows, int ncols, int pid, int np)
{
  if(pid == 0)
  {
    int ii = 0, jj = 0;

    ii = 0;
    for (jj = 0; jj < ncols; ++jj)
      m[ii*ncols + jj] = 100;

    jj = 0;
    for (ii = 1; ii < nrows-1; ++ii)
      m[ii*ncols + jj] = 0;

    jj = ncols-1;
    for (ii = 1; ii < nrows-1; ++ii)
      m[ii*ncols + jj] = 0;
  }
  else if(pid == np-1)
  {
    int ii = 0, jj = 0;
    ii = nrows-1;
    for (jj = 0; jj < ncols; ++jj)
      m[ii*ncols + jj] = 0;

    jj = 0;
    for (ii = 1; ii < nrows-1; ++ii)
      m[ii*ncols + jj] = 0;

    jj = ncols-1;
    for (ii = 1; ii < nrows-1; ++ii)
      m[ii*ncols + jj] = 0;
  }
  else
  {
    int ii = 0, jj = 0;
    jj = 0;
    for (ii = 1; ii < nrows-1; ++ii)
      m[ii*ncols + jj] = 0;

    jj = ncols-1;
    for (ii = 1; ii < nrows-1; ++ii)
      m[ii*ncols + jj] = 0;
  }
}

void print_screen(const Matrix & m, int nrows, int ncols, int pid, int np)
{
  // TODO
  MPI_Status status;
  int tag = 0;
  // Master pid prints
  if (0 == pid) {
    // print master data
    print_slab(m, nrows, ncols);
    // now receive in buffer and print other pids data
    Matrix buffer(nrows*ncols);
    // TODO
    for(int ii = 1; ii<np; ii++)
    {
      MPI_Recv(& buffer, nrows*ncols, MPI_DOUBLE, ii, tag, MPI_COMM_WORLD, &status);
      print_slab(buffer, nrows, ncols);
    }
  }
  else { // workers send
    // TODO
    MPI_Send(&m, ncols*nrows, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
  }
}

void print_slab(const Matrix & m, int nrows, int ncols)
{
  // ignore ghosts
  for(int ii = 1; ii < nrows-1; ++ii) {
    for(int jj = 0; jj < ncols; ++jj) {
      std::cout << std::setw(3) <<  m[ii*ncols + jj] << " ";
    }
    std::cout << "\n";
  }
}

/////////////////////////////////////////////////////
// SERIAL VERSIONS

void initial_conditions(Matrix & m, int nrows, int ncols)
{
  for(int ii=0; ii<nrows; ++ii) {
    for(int jj=0; jj<ncols; ++jj) {
      m[ii*ncols + jj] = 1.0;
    }
  }
}
void boundary_conditions(Matrix & m, int nrows, int ncols)
{
  int ii = 0, jj = 0;

  ii = 0;
  for (jj = 0; jj < ncols; ++jj)
    m[ii*ncols + jj] = 100;

  ii = nrows-1;
  for (jj = 0; jj < ncols; ++jj)
    m[ii*ncols + jj] = 0;

  jj = 0;
  for (ii = 1; ii < nrows-1; ++ii)
    m[ii*ncols + jj] = 0;

  jj = ncols-1;
  for (ii = 1; ii < nrows-1; ++ii)
    m[ii*ncols + jj] = 0;
}

void evolve(Matrix & m, int nrows, int ncols)
{
  for(int ii=0; ii<nrows; ++ii) {
    for(int jj=0; jj<ncols; ++jj) {
      // check if boundary
      if(ii == 0) continue;
      if(ii == nrows-1) continue;
      if(jj == 0) continue;
      if(jj == ncols-1) continue;
      // evolve non boundary
      m[ii*ncols+jj] = (m[(ii+1)*ncols + jj] +
                    m[(ii-1)*ncols + jj] +
                    m[ii*ncols + jj + 1] +
                    m[ii*ncols + jj - 1] )/4.0;
    }
  }
}

void print(const Matrix & m, double delta, int nrows, int ncols)
{
  for(int ii=0; ii<nrows; ++ii) {
    for(int jj=0; jj<ncols; ++jj) {
      std::cout << ii*delta << " " << jj*delta << " " <<  m[ii*ncols + jj] << "\n";
    }
    std::cout << "\n";
  }
}

void init_gnuplot(void)
{
  std::cout << "set contour " << std::endl;
  //std::cout << "set terminal gif animate " << std::endl;
  //std::cout << "set out 'anim.gif' " << std::endl;
}

void plot_gnuplot(const Matrix & m, double delta, int nrows, int ncols)
{
  std::cout << "splot '-' w pm3d " << std::endl;
  print(m, delta, nrows, ncols);
  std::cout << "e" << std::endl;
}