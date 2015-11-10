#include <iostream>
#include <armadillo>
#include <functional>

using namespace std;
using namespace arma;


double dx(vec t, double x, double y){
  //return -x + (2*y) + pow((x-y),2);
   return -x + 3*y;
}

double dy(vec t, double x, double y){
  //return pow((x-y),2) + y;
   return 3 * x + y;
}

void fbeuler_jb( double x0, double T, long N, double tol, long maxiter){
  double h = .0001;//(T+1)/N;
  vec t = linspace<vec>(0,h,T);
  vec x = x0 * ones<vec>(N);
  vec y = zeros<vec>(N);
  vec xlast = zeros<vec>(N);
  vec ylast = zeros<vec>(N);
  vec convInd = zeros<vec>(maxiter);
  long m = 0;

  do {
    xlast = x;
    for (long i = 0; i < N-1; i++){
      x(i+1) = xlast(i) + h * dx(t,xlast(i),ylast(i));
    }

    ylast = y;
    for (long i = N-1; i > 0; i--){
      y(i-1) = ylast(i) - h * dy(t,xlast(i),ylast(i));
    }

    convInd[m] = max(norm(y-ylast), norm(x-xlast));
    // printf("diff: %6.13f  y: %6.13f\n", convInd[m], y(0));

    m++;
  } while ((m < maxiter) && convInd[m-1] > tol);
  printf("diff: %6.13f  y: %6.13f\n", convInd[m-1], y(0));
}

int main(int argc, char** argv)
  {
  double T = 4.1;
  long N = 5000;
  double x0 = -0.01;
  double tol = pow(10, -7);
  long maxiter = 10000;

  fbeuler_jb(x0, T, N, tol, maxiter);

  return 0;
  }
