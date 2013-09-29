#include <iostream>
#include <cmath>
#include "kyle.h"       // Including the header files from each of us.
#include "arjun.h"      // (I've made blank header files for each of you
#include "jeremy.h"     // so we can practice, but it'd be good if you
                        // could either review or briefly look up how to
                        // work with header files.)

using namespace std;

int main()
{
  double number;

  // have the user enter a number
  cout << "Enter a number: ";
  cin >> number;

  // alter the number via our three functions and
  // display the results after each execution
  number = func_kyle(number);
  cout << "Kyle's number: " << number << endl;
  number = func_arjun(number);
  cout << "Arjun's number: " << number << endl;
  number = func_jeremy(number);
  cout << "Jeremy's number: " << number << endl;

  return 0;
}
