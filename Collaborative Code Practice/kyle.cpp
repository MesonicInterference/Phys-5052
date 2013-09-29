#include <cmath>
double func_kyle(double x) // my actual function
{
  x = sin(x);   // have to include cmath here as well if I want to use a sine
  x = pow(x,2);
  x = x - 1;
  return x;
}
