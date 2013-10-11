#include <iostream>
#include <cmath>
#include <vector>       // use vectors instead of arrays due to avoidance of
                        // issues with dynamic memory allocation
#include "LU.h"
#define SIZE 11
// future: include header files for:
//      > LU decomposition of A
//      > solving sysem of equations using L and U matrices

using namespace std;

  vector<double> position;                      // the positions from Fig 3.5
  vector<double> t;                             // the times from Fig 3.5

  vector<vector<double> > A;                    // matrix of coefficients
//   vector<vector<double> > L;                    // decomposed lower matrix
//   vector<vector<double> > U;                    // decomposed upper matrix
  vector<double> b;                             // resulting vector

  int order_;                                    // the order_ of the polynomial
                                                // used in the least-squares fit

// output the contents of A to the screen; mostly for debugging
void print_A()
{
  int row, column;

  for (row = 0; row < order_+1; row++)
  {for (column = 0; column < order_+1; column++)
     {cout << "A(" << row <<","<< column<< ") = "<<A[row][column] << endl;}}
}

// output the contents of L to the screen; mostly for debugging
void print_L()
{
  int row, column;

  for (row = 0; row < order_+1; row++)
  {for (column = 0; column < order_+1; column++)
     {cout << "L(" << row <<","<< column<< ") = "<<L[row][column] << endl;}}
}

// output the contents of U to the screen; mostly for debugging
void print_U()
{
  int row, column;

  for (row = 0; row < order_+1; row++)
  {for (column = 0; column < order_+1; column++)
     {cout << "U(" << row <<","<< column<< ") = "<<U[row][column] << endl;}}
}

// output the contents of b to the screen; mostly for debugging
void print_b()
{
  int row;

  for (row = 0; row < order_+1; row++)
  {  cout << "b(" << row <<") = "<<b[row] << endl;}
}

// initialize data
void initialize()
{
  int row, column, counter;                     // counter variables

  // resizing and initializing A, L, U, and b
  A.resize(order_+1);                            // need A, L, and U to be
//   L.resize(order_+1);                            // (order_+1) square matrices, so
//   U.resize(order_+1);                            // fill with (order_+1) vectors

  for (counter=0; counter < order_+1; counter++) // since A, L, and U are vectors
  {A[counter].resize(order_+1,0.0);              // of vectors, need to resize
//    L[counter].resize(order_+1,0.0);              // each element individually
//    U[counter].resize(order_+1,0.0);}
  }
    
  b.resize(order_+1,0.0);                        // need b to be a (order_+1)-term
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
}

// calculating the necessary sums
void summation()
{
  int row, column, counter;                     // counter variables

  // treat elements on a per-row basis
  for (row = 0; row < order_+1; row++)
  {
    // SETTING A
    // treat elements on a per-column basis
    for (column = 0; column < order_+1; column++)
    {
      // sum over the appropriate terms
      for (counter = 0; counter < SIZE; counter++)
      {
        A[row][column] += pow(t[counter],column+row);
/*
 Algorithm Explanation
 ---------------------
 Since the elements in our coefficient matrix follow a pattern qualitatively
 described by multiplying each successive row by another factor of t, the
 above algorithm handles this dynamically by noting that the index of t
 at each particular location is simply the sum of the row and column of that
 position (with the indexing necessarily starting at zero).  Thus, we may then
 simply vary the type of least-squares fit by altering the global parameter
 order_, which affects both this subroutine and intiailize().
 */
      }
    }

    // SETTING B
    // sum over the appropriate terms
    for (counter = 0; counter < order_+1; counter++)
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
}

// perform the least squares fits
void least_squares_fit(int x)
{
  order_ = x;                                    // set the order_ being used in the fit
  initialize();                                 // set everything back to normal
                                                // and appropriately resize the
                                                // matrices
  summation();                                  // perform the necessary summations

  print_A();
//   print_L();
//   print_U();
  print_b();
}

// the main algorithm
int main()
{
  int counter;                                  // a counter variable

  // do the linear (order_ = 1), quadratic (order_ = 2), and cubic (order_ = 3) fit
// //   for (counter=1; counter<4; counter++)
// //   {least_squares_fit(counter);}

  least_squares_fit(1);
//   LU_init(A, order_+1)

  LU_init(A, order_+1);
  LU_Solve(L, U, b, order_+1);

  for (int i = 0; i < (order_)+1; i++)
  {  cout << "x" << x[i] << endl;}
  
  return 0;
}
