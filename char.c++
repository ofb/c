#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include "mpfr.h"
#include "mpc.h"
#include "mpi.h"
#include "omp.h"

using namespace std;
using namespace NTL;

pair< unsigned int, unsigned int > validateParams(int argc, char *argv[], unsigned int lambdas[]) {
  unsigned int primeIndex;
  unsigned int iterationCount;
  // first argument is not used (it's the name of the program); second is the prime
  // index to start at; third is the number of primes to iterate on; all remaining
  // arguments are lambda values. Therefore a correct call includes 4 or more parameters.
  if (argc > 3)
    {
      primeIndex = atoi(argv[1]); // converts char array into int
      if (primeIndex < 1) {
	cout << "Please enter a positive index.\n";
	return make_pair(0, 0);
      }
      iterationCount = atoi(argv[2]);
      if (iterationCount < 1) {
	cout << "Please enter a positive number of iterations.\n";
	return make_pair(0, 0);
      }
      // |argv| = argc
      // first argument is the string of the program called then prime
      // (first assignment), then number of iteratons (second assignment)
      // and now lambdas
      for (int i = 3; i < argc; ++i) {
	lambdas[i-3] = atoi(argv[i]);
      }
    }
  else
    {
      cout << "Please call this program like \'char index iterations lambda1 lambda2...\' "
	   << "where index is the index of the prime you wish to start at, "
	   << "iterations equals the number of iterations you wish to run, "
	   << "and the list of lambdas equal the values for the lambda parameters"
	   << "you wish to run.\n"
	   << "We only have log tables for the first 10,000 primes, "
	   << "so please choose a value for the prime index <= 10,000.\n";
      return make_pair(0, 0);
    }
  
  return make_pair(primeIndex, iterationCount);
}

void getLogs(unsigned long p, unsigned long table[]) {
  stringstream ss;
  ss << p;
  // maybe need to change this to an absolute path once we run on a compute node
  string path = string("./logtables/table") + ss.str();
  ifstream infile(path.c_str());
  unsigned long token;
  for (unsigned long i = 0; i < p-1; ++i) {
    infile >> token;
    table[i] = token;
  }
  return;
}

void assignZeta(mpc_t &zeta, unsigned long p) {
  mpfr_t pi;
  mpfr_init2(pi, 64);
  mpfr_const_pi(pi,MPFR_RNDN);

  mpc_t zetaConstruct;
  mpc_init2(zetaConstruct, 64);
  mpc_set_ui(zetaConstruct, 1, MPFR_RNDN); // set to 1
  mpc_mul_i(zetaConstruct, zetaConstruct, 1, MPFR_RNDN); // multiply by i
  mpc_mul_fr(zetaConstruct, zetaConstruct, pi, MPFR_RNDN); // mult by MPFR
  mpc_mul_si(zetaConstruct, zetaConstruct, 2, MPFR_RNDN); // mult by long int
  mpc_div_ui(zetaConstruct, zetaConstruct, p-1, MPFR_RNDN); // divide by unsigned long int
  mpc_exp(zetaConstruct, zetaConstruct, MPFR_RNDN);

  mpc_set(zeta, zetaConstruct, MPFR_RNDN);  
  // printf("zeta is assigned to ");
  // mpc_out_str(stdout, 10, 0, zeta, MPFR_RNDN);
  // cout << endl;
  mpfr_clear(pi);
  mpc_clear(zetaConstruct);
  return;
}

// Here is our hack for the fact that we can't have vectors of mpc_t
// and we can't pass variable multidimensional arrays to functions:
// a one-dimensional array of length equal to the number of lambdas
// times (p-2). Yuk.
void fillV(const unsigned int lambdaLength,
	   unsigned long p,
	   mpc_t &zeta,
	   const unsigned int lambdas[],
	   unsigned long logtable[],
	   mpc_t evalV[]) {
  // We will only need to evaluate chi (by raising zeta to an element
  // of F) p times, which equals the order of the multiplicative
  // group. Though log has image in 0,...,p-2, n*log has codomain equal
  // to all of F, 0,...,p-1, and is usually onto.
  mpc_t primZetaEval[p];
  for (unsigned long n = 0; n < p; ++n) {
    // initialize the mpc_t objects in our array
    mpc_init2(primZetaEval[n], 64);
    // This 'evaluates' chi
    mpc_pow_ui(primZetaEval[n], zeta, n, MPFR_RNDN);
  }
  // There are p-1 elements in the multiplicative group of characters,
  // but we aren't interested in evaluating the trivial character.
  // We are using our horrible hack of collapsing what should be a
  // 2D array into a 1D array, so each lambda has (p-2) characters
  // associated with it.
  for (unsigned long n = 0; n < (p-2)*lambdaLength; ++n) {
    mpc_init2(evalV[n], 64);
    mpc_set_si(evalV[n], 0, MPFR_RNDN);
  }

  // Remember: the multiplicative characters form a group of order p-1,
  // not p.
  // log needs mod p; char exponent needs mod p-1
  ZZ_p::init(conv< ZZ > (p));
  ZZ_pBak bak;
  ZZ_p logArgZZp; // mod p
  ZZ_p chiArgZZp; // mod p-1
  ZZ intermediateZZ;
  // lambdas account for a shift in the evaluation,
  // as does the choice of chi. Since the chis form
  // a multiplicative group, we only need to evaluate
  // a primitive root to get the data for all of them.
  unsigned int lambda;
  #pragma omp parallel for schedule(static) default(private) shared(p, evalV, primZetaEval)
  for (unsigned long p1 = 0; p1 < p; ++p1) {
    for (unsigned long p2 = 0; p2 < p; ++p2) {
      for (unsigned long c = 1; c < p-1; ++c) {
	for (unsigned int l = 0; l < lambdaLength; ++l) {
	  // have to convert from ulong to ZZ to ZZ_p since
	  // ZZ_p doesn't directly accept ulong conversions
	  conv(logArgZZp, conv< ZZ > (p1+p2+p1*p1*p2+p1*p2*p2+p1*p1*p2*p2+lambdas[l]));
	  if (IsZero(logArgZZp)) continue;
	  // We find n*Log(a+lambda).
	  // Remember that the logtable index is given by the element of
	  // the group that of which you want the log minus one.
	  // Implicit promotion from ulong to ZZ_p on both arguments of mul.
	  // Explicit namespace to indicate where this function is coming from.
	  NTL::mul(intermediateZZ,
		   conv< ZZ > (c),
		   conv< ZZ > (logtable[conv< unsigned long > (logArgZZp-1)])
		   );
	  bak.save(); // saves the p modulus before switching to mod p-1
	  ZZ_p::init(conv< ZZ > (p-1));
	  conv(chiArgZZp, intermediateZZ); // chiArgument is now mod p-1
	  // We look up the evaluation of chi at this point.
	  // the primZetaEval array is actually canonically indexed; i.e.
	  // zeta^n is in the nth spot.	  
	  #pragma omp critical (summing) {
	    mpc_add(evalV[(p-2)*l+(c-1)], evalV[(p-2)*l+(c-1)],
		    primZetaEval[conv< unsigned long > (chiArgZZp)], MPFR_RNDN);
	  }
	  bak.restore(); // restore the modulus to p
	}
      }
    }
  }

  // clean up
  for (unsigned long n = 0; n < p-1; ++n) {
    mpc_clear(primZetaEval[n]);
  }
  return;
}

void print2DV(const unsigned int lambdaLength, unsigned long p, mpc_t V[]) {
  FILE *output;
  char output_filename[ 32 ];
  sprintf (output_filename, "out%06lu.txt", p);
  output = fopen ( output_filename, "w" );

  fprintf(output, "{{%lu},{", p);
  // iterate on the lambdas
  for (unsigned int a = 0; a < lambdaLength; ++a) {
    fprintf(output, "{");
    for (unsigned long b = 0; b < (p-2); ++b) {
      mpfr_out_str(output, 10, 0, mpc_realref(V[a*(p-2)+b]), MPFR_RNDN);
      fprintf(output, "+(");
      mpfr_out_str(output, 10, 0, mpc_imagref(V[a*(p-2)+b]), MPFR_RNDN);
      fprintf(output, "I)");
      if (b < p-3) fprintf(output, ",");
    }
    fprintf(output, "}");
    if (a < lambdaLength-1) fprintf(output, ",");
  }
  fprintf(output, "}}");
  fclose ( output );
  return;
}

void clear2DV(const unsigned long arrayLength, mpc_t V[]) {
  for (unsigned long i = 0; i < arrayLength; ++i) {
    mpc_clear(V[i]);
  }
  return;
}

void pCharSum(const unsigned long primeIndex,
	      const unsigned int lambdaLength,
	      const unsigned int lambdas[]) {
  // NTL prime upper bound is > 100 000 so we're safe
  PrimeSeq s;
  unsigned long p;
  for (unsigned long i = 0; i < primeIndex; ++i)
    p = s.next();
  unsigned long logtable[p-1];
  getLogs(p, logtable);
  
  mpc_t zeta;
  mpc_init2(zeta, 64);
  assignZeta(zeta, p);

  mpc_t lambdaChiV[lambdaLength*(p-2)];
  fillV(lambdaLength, p, zeta, lambdas, logtable, lambdaChiV);
  print2DV(lambdaLength, p, lambdaChiV);
  
  // cleanup
  mpc_clear(zeta);
  clear2DV(lambdaLength*(p-2), lambdaChiV);
  return;
}  

int main(int argc, char *argv[]) {
  // rank is the MPI process id, nprocs is the number of MPI processes.
  int rank, nprocs;

  const unsigned int lambdaLength = argc-3;
  unsigned int lambdas[argc-3];
  // params.first is primeIndex. params.second is iterationCount.
  pair< unsigned int, unsigned int > params = validateParams(argc, argv, lambdas);
  if (params.first == 0) return 0;

  // initialize MPI and assign values to rank and nprocs.
  MPI_Init ( &argc, &argv );
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
  MPI_Comm_size ( MPI_COMM_WORLD, &nprocs );

  for (unsigned int n = params.first; n < params.first + params.second; ++n) {
    if ( n % nprocs == rank ) {
      pCharSum(n, lambdaLength, lambdas);
    }
  }

  // Sync the MPI processes and quit
  MPI_Finalize();
  return 0;
}
