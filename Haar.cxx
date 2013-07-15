#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include "omp.h"

#define MAX_RESOLUTION 4294967295 // max unsigned int guaranteed range

typedef unsigned int uint;
typedef long long ll;
typedef unsigned long long ull;

using namespace std;

// Function prototypes
pair< uint, uint > validateParams(int argc, char *argv[]);
double avg(const double a, const double b);
void pushforward(const uint hRes, const uint pRes);

double avg(const double a, const double b) return (double) (a+b)/2;

pair< uint, uint > validateParams(int argc, char *argv[]) {
  uint hRes, pRes;
  if (argc == 3) {
    hRes = argv[1];
    pRes = argv[2];
    if (hRes < 1 || pRes < 1) {
      cout << "Please enter natural number resolution values.\n";
      return make_pair(0, 0);
    }
  } else {
    cout << "Please call this program with two parameters, the first being"
         << "the resolution of the Haar measure binning, and the second being"
         << "the resolution of the binning of the pushforward measure.\n";
    return make_pair(0,0);
  }

  return make_pair(hRes, pRes);
}

void pushforward(const uint hRes, const uint pRes) {
  double pfwd[pRes]; // our bins in the image
  // we are going to interate over hRes^2 bins of [0,2pi]^2
#pragma omp parallel for schedule(static) shared(pfwd)
  for 

int main(int argc, char *argv[]) {
  // the first return value is the haar resolution; the second is the pushforward resolution.
  pair< uint, uint > params = validateParams(argc, argv);
  
  // We bin [0, 3] in an array with the number of bins equal to the pushforward resolution.
  pushforward(params.first, params.second);

  return 0;
}

