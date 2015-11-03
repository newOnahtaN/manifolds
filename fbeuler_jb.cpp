#include <iostream>
#include <armadillo>
#include <functional>

using namespace std;
using namespace arma;


double dx(vec t, double x, double y){
  return (-x+2)*y + pow((x-y),2);
}

double dy(vec t, double x, double y){
  return pow((x-y),2) + y;
}

void fbeuler_jb( double x0, double T, double N, double tol, long maxiter){
  double h = (T+1)/N;
  //vector t = 0:h:T use span
  vec x = x0 * ones<vec>(5/*t.n_rows*/);
  vec y = zeros<vec>(5/*t.n_rows*/);
  vec xlast = x * 0;
  vec ylast = y * 0;

  vec convInd = zeros<vec>(maxiter);

  long m = 1;
  while ((m < maxiter) && max(norm(y-ylast), norm(x-xlast)) > tol){
    xlast = x;
    for (long i = 0; i < 5/*t.n_rows*/-1; i++ ){
      x[i-1] = xlast[i] - h * dx(/*t*/,xlast[i],ylast[i]);
    }

    ylast = y;
    for (long i = 5/*t.n_rows*/-1; i > 0; i-- ){
      y[i-1] = ylast[i] - h * dy(/*t*/,xlast[i],ylast[i]);
    }

    convInd[m] = max(norm(y-ylast), norm(x-xlast));
    m++;
  }
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
