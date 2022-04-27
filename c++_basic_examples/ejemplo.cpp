#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main()
{
  int Seed=1111;
  arma_rng::set_seed(Seed);
  mat A = randu<mat>(4,5);
  mat B = randu<mat>(4,5);

  cout << A*B.t() << endl;

  return 0;
}
