#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include "omp.h"

// sqrt of max unsigned long long guaranteed range
// also is max unsigned int guranteed range
#define MAX_RESOLUTION 4294967295

typedef unsigned int uint;
typedef long long ll;
typedef unsigned long long ull;

using namespace std;

// Function prototypes
uint validateParams(int argc, char *argv[]);
void pushforward(const uint hRes);
double trace(double theta1);
double haar(double theta1);
void print(double *tr, double *prob, const uint hRes);

uint validateParams(int argc, char *argv[]) {
  uint hRes;
  if (argc == 2) {
    hRes = argv[1];
    if (hRes < 1) {
      cout << "Please enter a natural number for the resolution.\n";
      return 0;
    }
    if (hRes >= MAX_RESOLUTION) {
      cout << "Please enter a natural number for the resolution"
	   << "less than " << MAX_RESOLUTION << "to prevent overflow.\n";
      return 0;
    }
  } else {
    cout << "Please call this program with one parameter,"
         << "the resolution of the Haar measure binning.";
    return 0;
  }

  return hRes;
}

double trace(double theta1) return 2*cos(theta1);
double haar(double theta1) return sin(theta1)*sin(theta1);

void pushforward(double *tr, double *prob, const uint hRes) {
  double binwidth = 2*M_PI/hRes;
  // we are going to interate over hRes^2 bins of [0,2pi]^2
#pragma omp parallel for schedule(static) shared(tr, prob)
  for (uint xbin = 0; xbin < hRes; ++xbin) {
    double theta1 = (xbin+0.5)*2*M_Pi/hRes; // our point in theta space
    tr[xbin] = trace(theta1);
    prob[xbin] = haar(theta1);
  }
  return;
}

void print(double *tr, double *prob, const uint hRes) {
  FILE *output;
  char output_filename[32];
  sprintf (output_filename, "haarOut");
  output = fopen ( output_filename, "w" );

  fprintf(output, "{");
  for (uint a = 0; a < hRes, ++a) {
    fprintf(output, "{%f,%f}", tr[a], prob[a]);
    if (a < (hRes-1)) fprintf(output, ",");
  }
  fprintf(output, "}");
  fclose ( output );
  return;
}

int main(int argc, char *argv[]) {
  // the first return value is the haar resolution; the second is the pushforward resolution.
  uint hRes = validateParams(argc, argv);
  if (!hRes) return 0;
  
  // We bin [0, 3] in an array with the number of bins equal to the pushforward resolution.
  double tr[hRes]; // the trace evaluated at a particular point
  double prob[hRes]; // the probability of that trace showing up
  pushforward(tr, prob, hRes);
  print(tr, prob, hRes);
  return 0;
}
