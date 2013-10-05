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

  vector<vector<double> > A;                    // matrix of coefficients
  vector<double> b;                             // resulting vector
  
  int order;                                    // the order of the polynomial
                                                // used in the least-squares fit
  
// initialize data
void initialize()
{
  int row, column;                              // counter variables
  
  A.resize(order+1);                            // need A to be a (order+1)
                                                // square matrix
  b.resize(order+1);                            // need b to be a (order+1)-term
                                                // column vector
  
  position.resize(SIZE);                        // give position the appropriate
                                                // size
  t.resize(SIZE);                               // give t the appropriate size
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
  
// fill A and b with zeroes
  for (row = 0; row < SIZE; row++)
    {b[row] = 0;
     for (column = 0; column < SIZE; column++)
      {A[row,column] = 0;}}
  
}

// calculating the necessary sums
void summation()
{
  int row, column, counter;                     // counter variables
  
  // treat elements on a per-row basis
  for (row = 0; row < SIZE; row++)
  {
    // SETTING A
    // treat elements on a per-column basis
    for (column = 0; column < SIZE; column++)
    {
      // sum over the appropriate terms
      for (counter = 0; counter < order+1; counter++)
      {
        A[row,column] += pow(t[counter],column+row);
/*
 Algorithm Explanation
 ---------------------
 Since the elements in our coefficient matrix follow a pattern qualitatively
 described by multiplying each successive row by another factor of t, the
 above algorithm handles this dynamically by noting that the index of t
 at each particular location is simply the sum of the row and column of that
 position (with the indexing necessarily starting at zero).  Thus, we may then
 simply vary the type of least-squares fit by altering the global parameter
 order, which affects both this subroutine and intiailize().
 */
      }
    }
    
    for (counter = 0; counter < order+1; counter++)
    {
      b[row] += position[counter] * pow(t[counter],row);
/*
 Algorithm Explanation
 ---------------------
 By the same method explained above, the sums comprising the resultant vector
 are determined by sums of the position data multiplied by increasing powers of
 t, with those powers being equivalent to, as before, the row+column indices,
 again beginning indexing at zero.  Since this is a column matrix, the column
 index is always zero, leaving the power of t being equal to the row index.
 */
    }
  }
      
      
      
//     sum_t     += t[counter];                    // summing t[i]
//     sum_t_2   += pow(t[counter],2);             // summing t[i]^2
//     sum_t_3   += pow(t[counter],3);             // summing t[i]^3
//     sum_t_4   += pow(t[counter],4);             // summing t[i]^4
//     sum_t_5   += pow(t[counter],5);             // summing t[i]^5
//     sum_t_6   += pow(t[counter],6);             // summing t[i]^6
//     
//     sum_p     += position[counter];             // summing p[i]
//     sum_p_t   += position[counter]
//                  * t[counter];                  // summing p[i] * t[i]
//     sum_p_t_2 += position[counter]
//                  * pow(t[counter],2);           // summing p[i] * t[i]^2
//     sum_p_t_3 += position[counter]
//                  * pow(t[counter],3);           // summing p[i] * t[i]^3
}

// find the linear fit
void lin_fit()
{
  order = 1;                                    // use a first order polynomial
  initialize();                                 // set everything back to normal
                                                // and appropriately resize the
                                                // matrices
  summation();                                  // perform the necessary summations
  
//   A.resize(2);
//   b.resize(2);
//   
//   A[0,0] = SIZE;
//   A[0,1] = sum_t;
//   A[1,0] = sum_t;
//   A[1,1] = sum_t_2;
//   
//   b[0,0] = sum_p;
//   b[1,0] = sum_p_t;
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