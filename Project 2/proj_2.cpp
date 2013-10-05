#include <iostream>
#include <cmath>
#include <vector>       // use vectors instead of arrays due to avoidance of
                        // issues with dynamic memory allocation
#define SIZE 11
// future: include header files for:
//      > LU decomposition
//      > solving sysem of equations using L and U matrices
//

using namespace std;

  vector<double> position;                      // the positions from Fig 3.5
  vector<double> t;                             // the times from Fig 3.5
  
  double sum_t   = 0;                           // sum of time
  double sum_t_2 = 0;                           // sum of time^2
  double sum_t_3 = 0;                           // sum of time^3
  double sum_t_4 = 0;                           // sum of time^4
  double sum_t_5 = 0;                           // sum of time^5
  double sum_t_6 = 0;                           // sum of time^6

  double sum_p     = 0;                         // sum of position
  double sum_p_t   = 0;                         // sum of position
  double sum_p_t_2 = 0;                         // sum of position
  double sum_p_t_3 = 0;                         // sum of position

  vector<vector<double> > A;                     // matrix of coefficients
  vector<double> b;                             // resulting vector
  
// initialize data
void initialize()
{
  position.resize(11);
  t.resize(11);
// initializing the positions
  position[0]  = 1.67203;
  position[1]  = 1.79792;
  position[2]  = 2.37791;
  position[3]  = 2.66408;
  position[4]  = 2.11245;
  position[5]  = 2.43969;
  position[6]  = 1.88843;
  position[7]  = 1.59447;
  position[8]  = 1.79634;
  position[9]  = 1.07810;
  position[10] = 0.21066;

// initializing the times
  t[0]  = 0.00;
  t[1]  = 0.10;
  t[2]  = 0.20;
  t[3]  = 0.30;
  t[4]  = 0.40;
  t[5]  = 0.50;
  t[6]  = 0.60;
  t[7]  = 0.70;
  t[8]  = 0.80;
  t[9]  = 0.90;
  t[10] = 1.00;
  
}

// calculating the necessary sums
void summation()
{
  int counter;                                  // a counter variable
  for (counter = 0; counter < SIZE; counter++)
  {
    sum_t     += t[counter];                    // summing t[i]
    sum_t_2   += pow(t[counter],2);             // summing t[i]^2
    sum_t_3   += pow(t[counter],3);             // summing t[i]^3
    sum_t_4   += pow(t[counter],4);             // summing t[i]^4
    sum_t_5   += pow(t[counter],5);             // summing t[i]^5
    sum_t_6   += pow(t[counter],6);             // summing t[i]^6
    
    sum_p     += position[counter];             // summing p[i]
    sum_p_t   += position[counter]
                 * t[counter];                  // summing p[i] * t[i]
    sum_p_t_2 += position[counter]
                 * pow(t[counter],2);           // summing p[i] * t[i]^2
    sum_p_t_3 += position[counter]
                 * pow(t[counter],3);           // summing p[i] * t[i]^3
  }
}

// find the linear fit
void lin_fit()
{
  A.resize(2);
  b.resize(2);
  
  A[0,0] = SIZE;
  A[0,1] = sum_t;
  A[1,0] = sum_t;
  A[1,1] = sum_t_2;
  
  b[0,0] = sum_p;
  b[1,0] = sum_p_t;
}

// find the quadratic fit
void quad_fit()
{
  
}

// find the cubic fit
void cubic_fit()
{
  
}

// the main algorithm
int main()
{
  initialize();
  summation();
  
  return 0;
}