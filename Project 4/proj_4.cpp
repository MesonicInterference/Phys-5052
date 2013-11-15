/*
  Project 3
  Phys 5052
  A. Kafle, J. Primus, K. Thomsen
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;

typedef double F(double,double,double); // defining a type allowing the program
                                        // to more closely parallel the book's
                                        // notation

vector<vector<double> > v_results;      // vector of vectors containing the
                                        // results of the calculations
                                        // row 0: t values
                                        // row 1: RK4 results

vector<vector<double> > x_results;      // vector of vectors containing the
                                        // results of the integration
                                        // row 0: t values
                                        // row 1: RK4 results

  double h;                             // step size
  int case_num;                         // case number returning different
                                        // derivatives via prime()
  int num_iterations;                   // number of iterations

  double epsilon;                       // the van der Pol oscillator constant

  double x0;                            // x(0)
  double v0;                            // x'(0)
  double low;                           // lower limit
  double high;                          // upper limit

// uses the fourth-order RK method to populate x_results and v_results
//
// input:
//      f - a renaming of the function input, here prime(), in order
//          to mirror the book's notation
// output:
//      none
void rk4(F f)
{
    double counter = 0;                 // a counter variable
    double f0a, f1a, f2a, f3a;          // intermediate function values for x
    double f0b, f1b, f2b, f3b;          // intermediate function values for v

// setting the initial conditions
    double x = x0;
    double v = v0;
    x_results[0][0] = 0;
    x_results[1][0] = x;
    v_results[0][0] = 0;
    v_results[1][0] = v;

// looping over time
    for (double t = 0; t < (high-low); t += h)
    {
// setting case_num so that f will return the derivative of x
      case_num = 0;
// intermediate function evaluations: x
      f0a = h * f(t, v, x);
      f1a = h * f(t + h / 2, v, x + f0a / 2);
      f2a = h * f(t + h / 2, v, x + f1a / 2);
      f3a = h * f(t + h, v, x + f2a);
// combining the intermediate function evaluations: x
      x += (f0a + 2 * f1a + 2 * f2a + f3a) / 6;
// setting the x results
      x_results[0][counter+1] = t;
      x_results[1][counter+1] = x;

// setting case_num so that f will return the derivative of v
      case_num = 1;
// intermediate function evaluations: v
      f0b = h * f(t, v, x);
      f1b = h * f(t + h / 2, v + f0b / 2, x);
      f2b = h * f(t + h / 2, v + f1b / 2, x);
      f3b = h * f(t + h, v + f2b, x);
// combining the intermediate function evaluations: v
      v += (f0b + 2 * f1b + 2 * f2b + f3b) / 6;
// setting the v results
      v_results[0][counter+1] = t;
      v_results[1][counter+1] = v;
// incrementing the counter
      counter++;}
}

// the two first-order DEs derived from the problem statement
//
// input:
//      t - time
//      x - x
//      v - x'
// output:
//      derivative of either x or v, depending on case_num
double prime(double t, double v, double x)
{
// first DE
  if(case_num == 0) {return v;}

// second DE
  if(case_num == 1) {return -x - epsilon * (pow(x,2) - 1) * v;}
}

// create output files for making graphs with gnuplot
//
// input:
//      none
// output:
//      none
void print_gnuplot()
{
// setting up x output files
  fstream f_rk4_pos;
  f_rk4_pos.open("rk4_pos.dat", ios::out);

// setting up v output files
  fstream f_rk4_vel;
  f_rk4_vel.open("rk4_vel.dat", ios::out);

// setting up 2D phase space output files
  fstream f_rk4_phase_2d;
  f_rk4_phase_2d.open("rk4_phase_2d.dat", ios::out);

// setting up 3D phase space output files
  fstream f_rk4_phase_3d;
  f_rk4_phase_3d.open("rk4_phase_3d.dat", ios::out);

// outputting RK4 results to files
  for (int counter = low; counter < num_iterations; counter++)
    {
// for x vs t plot: (t, x)
      f_rk4_pos << x_results[0][counter] << "" <<  ","
      <<" " << x_results[1][counter] << endl;
// for x' vs t plot: (t, x')
      f_rk4_vel << v_results[0][counter] << "" <<  ","
      <<" " << v_results[1][counter] << endl;
// for x vs x' plot: (x, x')
      f_rk4_phase_2d << x_results[1][counter]<< "" <<  ","
      <<" " << v_results[1][counter] << endl;
// for t vs x vs x' plot: (t, x, x')
      f_rk4_phase_3d << x_results[0][counter]<< "" <<  ","
      <<" " << x_results[1][counter]<< "" <<  ","
      <<" " << v_results[1][counter] << endl;
    }

// closing the created files
  f_rk4_vel.close();
  f_rk4_pos.close();
  f_rk4_phase_2d.close();
  f_rk4_phase_3d.close();
}

// calling the functions necessary to do the calculations
//
// input:
//      none
// output:
//      none
void computations()
{
// resize the results vectors based on the step size
    for (int counter = 0; counter < 2; counter++)
      {v_results[counter].resize((num_iterations+1),0.0);
       x_results[counter].resize((num_iterations+1),0.0);}
// RK4 calculation
    rk4(prime);
// output gnuplot files for streamlined plotting of results
    print_gnuplot();
}

int main()
{
// setting the lower and upper limits
      low = 0;                          // lower limit
      high = 8 * (2.0e0 * acos(0));     // upper limit

// user input for the step size
      cout << "What step size?\n";
      cout << "Recall that the interval of interest is t:[0,8π].\n";
      cin >> h;
      cout << endl;

// calculate the number of iterations to be done
      num_iterations = static_cast<int>(abs((high-low)/h));

// user input for ε
      cout << "What value for ε?\n";
      cin >> epsilon;
      cout << endl;

// user input for x(0)
      cout << "What value for x(0)?\n";
      cin >> x0;
      cout << endl;

// user input for x'(0)
      cout << "What value for x'(0)?\n";
      cin >> v0;
      cout << endl;

// resize vectors to allow for independent variable values and two different
// methods of solving the given DE
  v_results.resize(2);
  x_results.resize(2);
// use a function to perform the calculations
  computations();
}