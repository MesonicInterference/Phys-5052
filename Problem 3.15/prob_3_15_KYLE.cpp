#include <iostream>
#include <cmath>

using namespace std;

double period[6],		// periods in years
       dist[6],			// relative distances in AU
       l_period[6],		// natural logs of periods
       l_dist[6],		// natural logs of distances
       m,			// y-intercept of the linear fit
       n,			// slope of the linear fit; also the index of interest in this problem
       sum_l_dist = 0,		// sum of l_dist
       sum_l_dist_2 = 0,	// sum of l_dist^2
       sum_l_period = 0,	// sum of l_period
       sum_l_period_dist = 0;	// sum of l_dist * l_period

// set the lists as noted in Table 3.2
void initialize()
{
  // set the periods
  period[0] = .242;
  period[1] = .616;
  period[2] = 1.000;
  period[3] = 1.881;
  period[4] = 11.86;
  period[5] = 29.33;
  
  // set the distances
  dist[0] = .388;
  dist[1] = .724;
  dist[2] = 1.000;
  dist[3] = 1.524;
  dist[4] = 5.200;
  dist[5] = 9.510;
}

// take the natural logs of the given list of values
void logarithms()
{
  int counter;
  
  for (counter = 0; counter < 6; counter++)
  {
    l_period[counter] = log(period[counter]);
    l_dist[counter]   = log(dist[counter]);
  }
}

void summations()
{
  int counter;
  
  for (counter = 0; counter < 6; counter++)
  {
    sum_l_dist        += l_dist[counter];
    sum_l_dist_2      += pow(l_dist[counter],2);
    sum_l_period      += l_period[counter];
    sum_l_period_dist += l_period[counter] * l_dist[counter];
  }
}

int main()
{
  initialize();
  logarithms();
  summations();
  
  // calculate slopes and intercepts using least-squares relationships
  double slope = (sum_l_period_dist - sum_l_dist * sum_l_period / 6.0) / (sum_l_dist_2 - sum_l_dist*sum_l_dist / 6.0);
  double intercept = (sum_l_period - sum_l_dist * slope) / 6.0;
  
  cout << "linear fit intercept is: " << intercept << endl;
  cout << "linear fit slope is: " << slope << endl;
  
  return 0;
}
