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

typedef double F(double,double);        // defining a type allowing the program
                                        // to more closely parallel the book's
                                        // notation, namely renaming
                                        // v'(t,v) = f(t,v)

vector<vector<double> > vel_results;    // vector of vectors containing the
                                        // results of the calculations
                                        // row 0: t values
                                        // row 1: modified Euler results
                                        // row 2: RK4 results

vector<vector<double> > pos_results;    // vector of vectors containing the
                                        // results of the integration
                                        // row 0: t values
                                        // row 1: modified Euler results
                                        // row 2: RK4 results

vector< vector<double> > vel_err;       // vector of vectors containing the
                                        // actual errors in the calculations
                                        // row 0: t values
                                        // row 1: modified Euler errors
                                        // row 2: RK4 errors

vector< vector<double> > pos_err;       // vector of vectors containing the
                                        // actual errors in the integration
                                        // row 0: t values
                                        // row 1: modified Euler errors
                                        // row 2: RK4 errors

  double h;                             // step size
  int case_num;                         // case number for easily running
                                        // test cases

// constant declarations
  double g = 9.80665;                   // standard Earth gravity in m * s^(-2)
  double m = 1e-2;                      // mass of the projectile in kg
  double k = 1e-4;                      // drag coefficient for spherical
                                        // projectile in kg * m^(-1)

  double v0, a0, low, high;

//   double v0 = 0;                        // initial velocity in m * s^(-1)
//   double a0 = g - k / m * pow(v0,2);    // initial acceleration in m * s^(-2)
//   double low = 0, high = 10;            // lower and upper limits in s

// test case 1
//   double v0 = 0;
//   double a0 = 0;
//   double low = 0, high = 1;

// test case 2
//   double v0 = 1;
//   double a0 = 0;
//   double low = 0, high = 2.0e0 * acos(0);

// calculating v(t) via the modified Euler method
void mod_euler(F f, double h)
{
    double v = v0;
    double counter = 0;
    for (double t = 0; t < (high-low); t += h)
    {
      vel_results[0][counter] = t;
      vel_results[1][counter] = v;
      v += h * f(t+h/2, v+h/2);
      counter++;
    }
}

// calculating v(t) via RK4
void rk4(F f, double h)
{
    double v = v0;
    double counter = 0;
    double f0, f1, f2, f3;
    for (double t = low; t < (high-low); t += h)
    {
      f0 = h * f(t, v);
      f1 = h * f(t + h / 2, v + f0 / 2);
      f2 = h * f(t + h / 2, v + f1 / 2);
      f3 = h * f(t + h, v + f2);
      vel_results[2][counter] = v;
      v += (f0 + 2 * f1 + 2 * f2 + f3) / 6;
      counter++;}
}

// the DE expressed as an equation for the first derivative
double v_prime(double t, double v)
{
  if(case_num == 0) {return g - k / m * pow(v,2);}

  if(case_num == 1) {return pow(t,2) + 1;}

  if(case_num == 2) {return -sin(t);}
}

// the analytic solution to the DE
double v_actual(double t)
{
  if(case_num == 0) {return sqrt((m*g)/k) * tanh(sqrt((k*g)/m)*t);}

  if(case_num == 1) {return tan(t);}

  if(case_num == 2) {return cos(t);}
}

// the analytic solution to the integration
double x_actual(double t)
{
  if(case_num == 0) {return (sqrt((g*m)/k)*log(cosh(sqrt((g*k)/m)*t)))/sqrt((g*k)/m);}

  if(case_num == 1) {return -log(cos(t));}

  if(case_num == 2) {return sin(t);}
}

// integration using Simpson's rule to find the position
double simpson(double f0, double f1, double f2)
{
  double interval = h/2;
  return interval/3 * (f0 + 4*f1 + f2);
}

// Simpson's rule for integrating v(t) to get the position
void simpson_euler()
{
  for(int counter = 0; counter < (high-low)/h-2; counter++)
    { pos_results[0][counter] = vel_results[0][counter];
      if(counter == 0)
        {pos_results[1][counter] = simpson(vel_results[1][low+counter],
                                           vel_results[1][low+(counter+1)],
                                           vel_results[1][low+(counter+2)]);}
      else
        {pos_results[1][counter] = pos_results[1][counter-1] +
                                   simpson(vel_results[1][low+counter],
                                          vel_results[1][low+(counter+1)],
                                          vel_results[1][low+(counter+2)]);}}
}

// Simpson's rule for integrating v(t) to get the position
void simpson_rk4()
{
  for(int counter = 0; counter < (high-low)/h-2; counter++)
    {
      if(counter == 0)
        {pos_results[2][counter] = simpson(vel_results[2][low+counter],
                                           vel_results[2][low+(counter+1)],
                                           vel_results[2][low+(counter+2)]);}
      else
        {pos_results[2][counter] = pos_results[2][counter-1] +
                                   simpson(vel_results[2][low+counter],
                                          vel_results[2][low+(counter+1)],
                                          vel_results[2][low+(counter+2)]);}}
}

// print results to the terminal
void print_results()
{
  int counter, counter2;
    cout << "result tables for h = " << setprecision(2) << h << endl;
    cout << "    x    \t|\t modEuler\t|\t   RK4    \t|\n";
  for (counter = 0; counter < 10/h+1; counter++)
    {  for (counter2 = 0; counter2 <= 2; counter2++)
      {
        cout << fixed << setprecision(8) <<
        vel_results[counter2][counter] << "\t" <<  "|" <<"\t";
      }
      cout << endl;
    }
   cout << "------------------------------\n";
}

// calculate errors in v(t)
void vel_errors()
{
  for(int counter = 0; counter < 3; counter++)
    {vel_err[counter].resize(10/h+1);}
    
  for (int counter = 0; counter < 10/h+1; counter++)
    { vel_err[0][counter] = vel_results[0][counter];
      for (int counter2 = 1; counter2 <= 2; counter2++)
        {
          if (counter == 0) {vel_err[counter2][counter] = 0;} else {
          vel_err[counter2][counter] = fabs(v_actual(vel_results[0][counter])
                - vel_results[counter2][counter]);}
        }
    }
}

// calculate errors in v(t)
void pos_errors()
{
  for(int counter = 0; counter < 3; counter++)
    {pos_err[counter].resize(10/h+1);}
    
  for (int counter = 0; counter < 10/h+1; counter++)
    { pos_err[0][counter] = pos_results[0][counter];
      for (int counter2 = 1; counter2 <= 2; counter2++)
        {
          if (counter == 0) {pos_err[counter2][counter] = 0;} else {
          pos_err[counter2][counter] = fabs(x_actual(pos_results[0][counter])
                - pos_results[counter2][counter]);}
        }
    }
}

// print errors to the terminal
void print_errors()
{
  int counter, counter2;
    cout << "errors for h = " << setprecision(2) << h << endl;
    cout << "    x    \t|\t  Euler \t|\t modEuler\t|\t impEuler\t|\t   RK4    \t|\n";
  for (counter = 0; counter < 10/h+1; counter++)
    {  for (counter2 = 0; counter2 < 3; counter2++)
      {
        cout << fixed << setprecision(8) <<
        vel_err[counter2][counter] << "\t" <<  "|" <<"\t";
      }
      cout << endl;
    }
   cout << "------------------------------\n";
}

// create output files for making graphs with gnuplot
void print_gnuplot()
{
// setting up velocity output files
  fstream f_modeuler_vel;
  f_modeuler_vel.open("mod_euler_vel.dat", ios::out);
  fstream f_rk4_vel;
  f_rk4_vel.open("rk4_vel.dat", ios::out);
  fstream f_vel;
  f_vel.open("vel.dat", ios::out);

// setting up velocity error output files
  fstream f_modeuler_vel_err;
  f_modeuler_vel_err.open("mod_euler_vel_err.dat", ios::out);
  fstream f_rk4_vel_err;
  f_rk4_vel_err.open("rk4_vel_err.dat", ios::out);

// setting up position output files
  fstream f_modeuler_pos;
  f_modeuler_pos.open("mod_euler_pos.dat", ios::out);
  fstream f_rk4_pos;
  f_rk4_pos.open("rk4_pos.dat", ios::out);
  fstream f_pos;
  f_pos.open("pos.dat", ios::out);

// setting up position error output files
  fstream f_modeuler_pos_err;
  f_modeuler_pos_err.open("mod_euler_pos_err.dat", ios::out);
  fstream f_rk4_pos_err;
  f_rk4_pos_err.open("rk4_pos_err.dat", ios::out);
  
  int counter, counter2;                // counter variables

  
  for (counter = low; counter < high/h+1; counter++)
  {
      f_vel << vel_results[0][counter] << "" << 
      "," <<" " << v_actual(vel_results[0][counter]) << endl;
      f_pos << pos_results[0][counter] << "" << 
      "," <<" " << x_actual(pos_results[0][counter]) << endl;
  }
  
// outputting modified Euler results to files
  for (counter = low; counter < high/h; counter++)
    {
      f_modeuler_vel << vel_results[0][counter] << "" << 
      "," <<" " << vel_results[1][counter] << endl;
      f_modeuler_vel_err << vel_err[0][counter] << "" << 
      "," <<" " << vel_err[1][counter] << endl;
      f_modeuler_pos << pos_results[0][counter] << "" << 
      "," <<" " << pos_results[1][counter] << endl;
      f_modeuler_pos_err << pos_err[0][counter] << "" << 
      "," <<" " << pos_err[1][counter] << endl;
    }

// outputting RK4 results to files
  for (counter = low; counter < high/h; counter++)
    {
      f_rk4_vel << vel_results[0][counter] << "" <<  ","
      <<" " << vel_results[2][counter] << endl;
      f_rk4_vel_err << vel_results[0][counter] << "" << 
      "," <<" " << vel_err[2][counter] << endl;
      f_rk4_pos << pos_results[0][counter] << "" <<  ","
      <<" " << pos_results[2][counter] << endl;
      f_rk4_pos_err << pos_err[0][counter] << "" << 
      "," <<" " << pos_err[2][counter] << endl;
    }

  f_modeuler_vel.close();
  f_rk4_vel.close();
  f_modeuler_vel_err.close();
  f_rk4_vel_err.close();
  f_modeuler_pos.close();
  f_rk4_pos.close();
  f_modeuler_pos_err.close();
  f_rk4_pos_err.close();
}

// allowing for easy switching among test cases
void choice()
{
  cout << "select a case:\n";
  cout << "case 0: v'(t,v) = g - k / m * v(t)^2\n";
  cout << "case 1: v'(t,v) = v(t)^2 + 1\n";
  cout << "case 2: v'(t,v) = -sin(t)\n";
  cin >> case_num;
  while((case_num != 0) && (case_num != 1) && (case_num != 2))
    {cin >> case_num;}
    
  if(case_num == 0)
    { v0 = 0;                        // initial velocity in m * s^(-1)
      a0 = g - k / m * pow(v0,2);    // initial acceleration in m * s^(-2)
      low = 0;
      high = 10;}           // lower and upper limits in s
  if(case_num == 1)
    { v0 = 0;
      a0 = 0;
      low = 0;
      high = 1;}
  if(case_num == 2)
    { v0 = 1;
      a0 = 0;
      low = 0;
      high = 2.0e0 * acos(0);}
}

// calling the functions necessary to do the necessary calculations
void computations()
{
    for (int counter = 0; counter <= 2; counter++)
      {vel_results[counter].resize(((high-low)/h+1),0.0);
       pos_results[counter].resize(((high-low)/h+1),0.0);}
    mod_euler(v_prime, h);
    rk4(v_prime, h);
    vel_errors();
    simpson_euler();
    simpson_rk4();
    pos_errors();
    print_gnuplot();
}



int main()
{
  choice();
  vel_results.resize(3);
  pos_results.resize(3);
  vel_err.resize(3);
  pos_err.resize(3);
  h = .05;
  computations();
}