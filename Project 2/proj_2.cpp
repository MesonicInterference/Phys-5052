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

  vector<vector<double> > A;                    // matrix of coefficients
  vector<double> b;                             // resulting vector
  
  int order;                                    // the order of the polynomial
                                                // used in the least-squares fit

// output the contents of A to the screen; mostly for debugging
void print_A()
{
  int row, column;

  for (row = 0; row < order+1; row++)
  {for (column = 0; column < order+1; column++)
     {cout << "A(" << row <<","<< column<< ") = "<<A[row][column] << endl;}}
}

// output the contents of b to the screen; mostly for debugging
void print_b()
{
  int row;

  for (row = 0; row < order+1; row++)
  {  cout << "b(" << row <<") = "<<b[row] << endl;}
}


// initialize data
void initialize()
{
  int row, column, counter;                     // counter variables
  
  // resizing and initializing A and b
  A.resize(order+1);                            // need A to be a (order+1)
                                                // square matrix, so fill A with
                                                // (order+1) vectors

  for (counter=0; counter < order+1; counter++) // since A is a vector of
  {A[counter].resize(order+1,0.0);}             // vectors, need to resize each
                                                // element individually

  b.resize(order+1,0.0);                        // need b to be a (order+1)-term
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
  for (row = 0; row < order+1; row++)
  {
    // SETTING A
    // treat elements on a per-column basis
    for (column = 0; column < order+1; column++)
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
 order, which affects both this subroutine and intiailize().
 */
      }
    }
    
    // SETTING B
    // sum over the appropriate terms
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
}

// use the least squares fits
void least_squares_fit(int x)
{
  order = x;                                    // use a first order polynomial
  initialize();                                 // set everything back to normal
                                                // and appropriately resize the
                                                // matrices
  summation();                                  // perform the necessary summations

  print_A();
  print_b();
}

// // find the quadratic fit
// void quad_fit()
// {
//   order = 2;                                    // use a first order polynomial
//   initialize();                                 // set everything back to normal
//                                                 // and appropriately resize the
//                                                 // matrices
//   summation();                                  // perform the necessary summations
// 
//   print_A();
//   print_b();
// }
// 
// // find the cubic fit
// void cubic_fit()
// {
//   order = 3;                                    // use a first order polynomial
//   initialize();                                 // set everything back to normal
//                                                 // and appropriately resize the
//                                                 // matrices
//   summation();                                  // perform the necessary summations
// 
//   print_A();
//   print_b();
// }

// the main algorithm
int main()
{
  int counter;                                  // a counter variable
  
  // do the linear (order = 1), quadratic (order = 2), and cubic (order = 3) fit
  for (counter=1; counter<4; counter++)
  {least_squares_fit(counter);}
  
  return 0;
}