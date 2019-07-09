#include <fstream>
#include <sstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
#include <cfloat>

using namespace std;

int main() {
  long double var,varnew,div;
  int i;
  var = div = 10;
  cout.precision(32);
  for (i=1;i<=100;i++) {
    cout << i << " " << var << endl;
    div /= 10;
    varnew = var + div;
    if (varnew == var) break;
    else var = varnew;
  }

  return 0;
}
