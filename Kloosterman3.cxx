#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <NTL/ZZ.h>
#include "mpfr.h"
#include "mpc.h"
#include "mpi.h"
#include "omp.h"

#define MAX_PRIME_INDEX 180000

typedef unsigned long long ull;

using namespace std;
using namespace NTL;

pair< unsigned int, unsigned int > validateParams(int argc, char *argv[]) {
  unsigned int primeIndex;
  unsigned int iterationCount;
  // first argument is not used (it's the name of the program); second is the prime
  // index to start at; third is the number of primes to iterate on; all remaining
  // arguments are lambda values. Therefore a correct call includes 4 or more parameters.
  if (argc > 2)
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
    }
  else
    {
      cout << "Please call this program like \'Kloosterman3 index iterations\' "
	   << "where index is the index of the prime you wish to start at and "
	   << "iterations equals the number of iterations you wish to run."
	   << "Calculations in this program overflow for a prime > " << MAX_PRIME_INDEX << ","
	   << "so please choose a value for the prime index <= " << MAX_PRIME_INDEX << ".\n";
      return make_pair(0, 0);
    }
  if (primeIndex + iterationCount > MAX_PRIME_INDEX) {
	  cout << "Make sure that your prime index does not exceed " << MAX_PRIME_INDEX << "\n";
	  return make_pair(0, 0);
  }
  
  return make_pair(primeIndex, iterationCount);
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
  mpc_div_ui(zetaConstruct, zetaConstruct, p, MPFR_RNDN); // divide by unsigned long int
  mpc_exp(zetaConstruct, zetaConstruct, MPFR_RNDN);

  mpc_set(zeta, zetaConstruct, MPFR_RNDN);  
  mpfr_clear(pi);
  mpc_clear(zetaConstruct);
  return;
}

void fillV(unsigned long p,
	   mpc_t &zeta,
	   mpc_t eval) {
  mpc_t primZetaEval[p];
  for (unsigned long n = 0; n < p; ++n) {
    // initialize the mpc_t objects in our array
    mpc_init2(primZetaEval[n], 64);
    // This 'evaluates' psi
    mpc_pow_ui(primZetaEval[n], zeta, n, MPFR_RNDN);
  }
  // There are p elements in the multiplicative group of characters,
  // but we aren't interested in evaluating the trivial character.
  mpc_init2(eval, 64);
  mpc_set_si(eval, 0, MPFR_RNDN);

  // Remember: the additive characters form a group of order p.
  // zeta exponent needs mod p.
#pragma omp parallel for schedule(static) shared(p, eval, primZetaEval)
  for (unsigned long p1 = 0; p1 < p; ++p1) {
    ull p1ull, p2ull, p3ull;
    unsigned long zetaPower;
    p1ull = (ull) p1;
    for (unsigned long p2 = 0; p2 < p; ++p2) {
      p2ull = (ull) p2;
      for (unsigned long p3 = 0; p3 < p; ++p3) {
	p3ull = (ull) p3;
	if ( (((ull) p1ull*p2ull*p3ull) % (ull) p) != 1) continue;
        // our polynomial is p1+p2+p3
        // Remember, additive characters form a multiplicative group of order p.
        zetaPower = (p1ull+p2ull+p3ull) % (ull) p;
        // We look up the evaluation of psi at this point.
        // the primZetaEval array is actually canonically indexed; i.e.
        // zeta^n is in the nth spot.	  
	#pragma omp critical (summing)
	{
	  mpc_add(eval, eval, primZetaEval[zetaPower], MPFR_RNDN);
	}
      }
    }
  }

  // clean up
  for (unsigned long n = 0; n < p; ++n) {
    mpc_clear(primZetaEval[n]);
  }
  return;
}

void printSum(unsigned long p, mpc_t V) {
  FILE *output;
  char output_filename[ 32 ];
  sprintf (output_filename, "Kloosterman3/kloos%06lu", p);
  output = fopen ( output_filename, "w" );

  fprintf(output, "{{%lu},{", p);
  mpfr_out_str(output, 10, 0, mpc_realref(V), MPFR_RNDN);
  fprintf(output, "+(");
  mpfr_out_str(output, 10, 0, mpc_imagref(V), MPFR_RNDN);
  fprintf(output, "I)");
  fprintf(output, "}}");
  
  fclose ( output );
  return;
}

void pCharSum(const unsigned long primeIndex) {
  // NTL prime upper bound is > 100 000 so we're safe
  PrimeSeq s;
  unsigned long p;
  for (unsigned long i = 0; i < primeIndex; ++i)
    p = s.next();
  
  mpc_t zeta;
  mpc_init2(zeta, 64);
  assignZeta(zeta, p);

  fillV(p, zeta, sum);
  printSum(p, sum);
  
  // cleanup
  mpc_clear(sum)
  mpc_clear(zeta);
  return;
}

int main(int argc, char *argv[]) {
  // rank is the MPI process id, nprocs is the number of MPI processes.
  int rank, nprocs;

  // params.first is primeIndex. params.second is iterationCount.
  pair< unsigned int, unsigned int > params = validateParams(argc, argv);
  if (params.first == 0) return 0;

  // initialize MPI and assign values to rank and nprocs.
  MPI_Init ( &argc, &argv );
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
  MPI_Comm_size ( MPI_COMM_WORLD, &nprocs );

  // get name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  
  printf("Hello from processor %s, rank %d out of %d processors\n", processor_name, rank, nprocs);
  for (unsigned int n = params.first; n < params.first + params.second; ++n) {
    if ( n % nprocs == rank ) {
      pCharSum(n);
    }
  }

  // Sync the MPI processes and quit
  MPI_Finalize();
  return 0;
}
