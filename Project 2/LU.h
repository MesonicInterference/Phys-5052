#ifndef LU_run
#define LU_run

/***************************************************************************
LU Decomposition
by Jeremy Primus
10/9/13
This function solves the linear set of equations Ax = b
by the method of LU decomposition
***************************************************************************/

#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

//global variable declaration
vector< vector<double> > L;//(2, vector<double>(2,1));
vector< vector<double> > U;//(2, vector<double>(2,1));
// vector< vector<double> > A(2, vector<double>(2,1));
// vector<double> b(2,1);
vector<double> x;
vector<double> order;



void LU_init(vector< vector<double> > a, int N)
//function to compute lower triangular and upper tiangular
//matrices for LU Decomposition
{
//declarations
vector<double> scale;
double max, imax, det;


  // resizing and initializing A, L, U, and b
                                                 // need A, L, and U to be
  L.resize(N);                            // (N) square matrices, so
  U.resize(N);                            // fill with (N) vectors
  order.resize(N);
  scale.resize(N);
  
  for (int counter=0; counter < N+1; counter++) // since A, L, and U are vectors
  {                                              // of vectors, need to resize
   L[counter].resize(N,0.0);              // each element individually
   U[counter].resize(N,0.0);}

det = 1;

for (int i = 0; i < N; i++)  //determine a scaling factor for each row
{
  order[i] = i+1;
  double max = 1;
  for (int j = 0; j < N; j++)
  { if (fabs(a[i][j]) > max)
    {max = fabs(a[i][j]);}}
  scale[i] = 1/max;
}


for (int k = 0; k < N; k++)  //calculate a column of L
{
  for (int i = k; i < N; i++)
  {
    double sumL;
    sumL = a[i][k];
    cout << "i = " << i << endl;
    cout << "k = " << k << endl;
    for (int j = 0; j < k+1; j++)
    { sumL = sumL - a[i][j]*a[j][k];
  cout << "finished the j loop" << " " << i << " " << j << " " << k << endl;}
    L[i][k] = sumL;
  cout << "it's done something" << " " << i << endl;
  }

  cout << "finished that loop" << endl;

  max = 0;

  for (int i = k; i < N; i++)  //check for largest element of computed column of L
  {
    if (scale[i]*L[i][k] >= max)
    {
      max = scale[i]*L[i][k];
      imax = i;
    }
  }

  if (imax != k){}  //The following is not necessary if max is already on diagonal
    else
    {
      det = -det;
      for (int j = 0; j < N; j++)
      {
        double temp;
        temp = L[imax][j];
        L[imax][j] = L[k][j];
        L[k][j] = temp;
      }

//scale row
      double temporary;
      temporary = scale[imax];
      scale[imax] = scale[k];
      scale[k] = temporary;

//reorder
      double itemp;
      itemp = order[imax];
      order[imax] = order[k];
      order[k] = itemp;
  }

  det = det*L[k][k];

//compute U matrix
  double sumU;
  if (k == 1)
  {
    for (int j = 1; j < N; j++)
    {  a[0][j] = a[0][j]/a[0][0]; }
  }
    else
    {
      for (int j = k+1; j < N; j++)
      {
        sumU = a[k][j];
        for (int i = 0; i < k; i++)
        {  sumU = sumU - a[k][i]*a[i][j];}
      U[k][j] = sumU/a[k][k];
      }
    }
}

//Last element of L
  double sumLL;
  sumLL = a[N][N];
  for (int j = 0; j < N; j++)
  {  sumLL = sumLL - a[N][j]*a[j][N];}
  L[N][N] = sumLL;
  det = det* a[N][N];
return;
}

void LU_Solve(vector< vector<double> > L, vector< vector<double> > U, vector<double> b, int N)
//Function to solve for matrix 'x' in matrix equations of the form Ax = b
//given a lower triangular matrix L, and upper triangular matrix U
//have been initialized by the function LU_init()
{
//variable declarations
  double sumL, sumU;

  for (int i = 0; i < N; i ++)
  { x[i] = b[order[i]];}

  x[0] = x[0]/L[0][0];
  for (int i =1; i < N; i++)
  {
    sumL = x[i];
    for (int k = 1; k < i; k++)
    {  sumL = sumL - L[i][k]*x[k]; }
    x[i] = sumL/L[i][i];
  }

  for (int i = 0; i < N; i++)
  {
    sumU = x[N-i];
    for (int k = N-i; k < N; k++)
    {  sumU = sumU - U[N-i][k]*x[k];}
    x[N-i] = sumU;
  }
return;
}

// int main()
// {
// LU_init(A, N);
// LU_Solve(L, U, b, N);
// 
// for (int i = 0; i < N+1; i++)
// {  cout << "x" << x[i] << endl;}

// return 0;
// }

#endif