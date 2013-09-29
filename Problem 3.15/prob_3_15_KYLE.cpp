#include <iostream>
#include <cmath>

using namespace std;

double velocity[11],            // velocities
       t[11],                   // times
       m,                       // y-intercept of the linear fit
       n,                       // slope of the linear fit; also the index of interest in this problem
       sum_t = 0,               // sum of t
       sum_t_2 = 0,             // sum of t^2
       sum_velocity = 0,        // sum of velocities
       sum_velocity_t = 0;      // sum of t * velocities

// set the lists as noted in Table 3.2
void initialize()
{
  // set the velocities
  velocity[0]  = -0.10290;
  velocity[1]  =  0.37364;
  velocity[2]  =  2.43748;
  velocity[3]  =  3.93836;
  velocity[4]  =  3.31230;
  velocity[5]  =  5.49472;
  velocity[6]  =  5.43325;
  velocity[7]  =  6.39321;
  velocity[8]  =  9.06048;
  velocity[9]  =  9.36416;
  velocity[10] =  9.52066;
  
  // set the times
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

void summations()
{
  int counter;
  
  for (counter = 0; counter < 11; counter++)
  {
    sum_t          += t[counter];
    sum_t_2        += pow(t[counter],2);
    sum_velocity   += velocity[counter];
    sum_velocity_t += velocity[counter] * t[counter];
  }
}

int main()
{
  initialize();
  summations();
  
  // calculate slopes and intercepts using least-squares relationships
  double slope = (sum_velocity_t - sum_t * sum_velocity / 12.0) / (sum_t_2 - sum_t*sum_t / 12.0);
  double intercept = (sum_velocity - sum_t * slope) / 12.0;
  
  cout << "linear fit intercept is: " << intercept << endl;
  cout << "linear fit slope is: " << slope << endl;
  
  return 0;
}
