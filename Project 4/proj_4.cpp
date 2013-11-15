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

typedef double F(double,double,double);        // defining a type allowing the program
                                        // to more closely parallel the book's
                                        // notation, namely renaming
                                        // v'(t,v) = f(t,v)

vector<vector<double> > vel_results;    // vector of vectors containing the
                                        // results of the calculations
                                        // row 0: t values
                                        // row 1: RK4 results

vector<vector<double> > pos_results;    // vector of vectors containing the
                                        // results of the integration
                                        // row 0: t values
                                        // row 1: uses RK4 results

  double h;                             // step size
  int case_num;                         // case number for easily running
                                        // test cases
  int num_iterations;                   // number of iterations

// constant declarations

  double epsilon;                       // the van der Pol oscillator constant

  double v0, x0, low, high;

// calculating v(t) via RK4
void rk4(F f, double h)
{
    double counter = 0;
    double f0a, f1a, f2a, f3a;
    double f0b, f1b, f2b, f3b;

// setting the initial conditions
    double x = x0;
    double v = v0;
    pos_results[0][0] = 0;
    pos_results[1][0] = x;
    vel_results[0][0] = 0;
    vel_results[1][0] = v;


    for (double t = 0; t < (high-low); t += h)
    {
// intermediate function evaluations: position
      case_num = 0;
      f0a = h * f(t, v, x);
      f1a = h * f(t + h / 2, v, x + f0a / 2);
      f2a = h * f(t + h / 2, v, x + f1a / 2);
      f3a = h * f(t + h, v, x + f2a);
      x += (f0a + 2 * f1a + 2 * f2a + f3a) / 6;
      pos_results[0][counter+1] = t;
      pos_results[1][counter+1] = x;

// intermediate function evaluations: velocity
      case_num = 1;
      f0b = h * f(t, v, x);
      f1b = h * f(t + h / 2, v + f0b / 2, x);
      f2b = h * f(t + h / 2, v + f1b / 2, x);
      f3b = h * f(t + h, v + f2b, x);
      v += (f0b + 2 * f1b + 2 * f2b + f3b) / 6;
      vel_results[0][counter+1] = t;
      vel_results[1][counter+1] = v;
// setting result first so that the term indexed 0 will be equal to v0
// combining the intermediate function evaluations
      counter++;}
}

// the DE expressed as an equation for the first derivative
// inputs:
//      t - time
//      x - position
//      v - x'
double v_prime(double t, double v, double x)
{
// first DE
  if(case_num == 0) {return v;}

// second DE
  if(case_num == 1) {return -x - epsilon * (pow(x,2) -1) * v;}
}

// create output files for making graphs with gnuplot
void print_gnuplot()
{
// setting up velocity output files
  fstream f_rk4_vel;
  f_rk4_vel.open("rk4_vel.dat", ios::out);

// setting up position output files
  fstream f_rk4_pos;
  f_rk4_pos.open("rk4_pos.dat", ios::out);
  
// setting up phase space output files
  fstream f_rk4_phase_2d;
  f_rk4_phase_2d.open("rk4_phase_2d.dat", ios::out);
  
// setting up phase space output files
  fstream f_rk4_phase_3d;
  f_rk4_phase_3d.open("rk4_phase_3d.dat", ios::out);

// outputting RK4 results to files
  for (int counter = low; counter < num_iterations; counter++)
    {
      f_rk4_vel << vel_results[0][counter] << "" <<  ","
      <<" " << vel_results[1][counter] << endl;
      f_rk4_pos << pos_results[0][counter] << "" <<  ","
      <<" " << pos_results[1][counter] << endl;
      f_rk4_phase_2d << pos_results[1][counter]<< "" <<  ","
      <<" " << vel_results[1][counter] << endl;
      f_rk4_phase_3d << pos_results[0][counter]<< "" <<  ","
      <<" " << pos_results[1][counter]<< "" <<  ","
      <<" " << vel_results[1][counter] << endl;
    }

// close files
  f_rk4_vel.close();
  f_rk4_pos.close();
}

// calling the functions necessary to do the calculations
void computations()
{
// resize the results vectors based on the step size
    for (int counter = 0; counter < 2; counter++)
      {vel_results[counter].resize((num_iterations+1),0.0);
       pos_results[counter].resize((num_iterations+1),0.0);}
// RK4 calculation
    rk4(v_prime, h);
// output gnuplot files for streamlined plotting of results
    print_gnuplot();
}

int main()
{
// setting the step size for the RK method
  h = .01;
  cout << "h = " << h << endl;

      low = 0;                          // lower limit in s
      high = 8 * (2.0e0 * acos(0));     // upper limit in s
      num_iterations = static_cast<int>(abs((high-low)/h));

      cout << "What value for ε?";
      cin >> epsilon;
      cout << endl;

      cout << "What value for x(0)?";
      cin >> x0;
      cout << endl;

      cout << "What value for x'(0)?";
      cin >> v0;
      cout << endl;


// resize vectors to allow for independent variable values and two different
// methods of solving the given DE
  vel_results.resize(2);
  pos_results.resize(2);
  case_num = 0;
  computations();
  case_num = 1;
  computations();
}