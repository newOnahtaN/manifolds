#include <iostream>
#include <armadillo>
#include <functional>

using namespace std;
using namespace arma;


double dx(vec t, double x, double y){
  return -x + (2*y) + pow((x-y),2);
}

double dy(vec t, double x, double y){
  return pow((x-y),2) + y;
}

void fbeuler_jb( double x0, double T, long N, double tol, long maxiter){
  double h = (T+1)/N;
  vec t = linspace<vec>(0,h,T);
  vec x = x0 * ones<vec>(t.n_rows);
  vec y = zeros<vec>(t.n_rows);
  vec xlast = zeros<vec>(t.n_rows);;
  vec ylast = zeros<vec>(t.n_rows);;
  vec convInd = zeros<vec>(maxiter);
  long m = 0;

  do {
    xlast = x;
    for (long i = 0; i < t.n_rows-1; i++){
      x(i+1) = xlast(i) + h * dx(t,xlast(i),ylast(i));
    }

    ylast = y;
    for (long i = t.n_rows-1; i > 0; i--){
      y(i-1) = ylast(i) - h * dy(t,xlast(i),ylast(i));
    }

    convInd[m] = max(norm(y-ylast), norm(x-xlast));
    printf("%6.13f\n", convInd[m]);

    m++;
  } while ((m < maxiter) && convInd[m-1] > tol);
}

int main(int argc, char** argv)
  {
  // mat A = randu<mat>(4,5);
  // mat B = randu<mat>(4,5);
  // vec x = 5 * ones<vec>(6);
  // vec y = zeros<vec>(6);
  //
  // cout << A*B.t() << endl;
  // cout << y << endl;
  double T = 4.1;
  long N = 310000;
  double x0 = -0.01;
  double tol = pow(10, -7);
  long maxiter = 10000;

  fbeuler_jb(x0, T, N, tol, maxiter);

  return 0;
  }
