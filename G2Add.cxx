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

typedef unsigned long long ull;

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
	   << "you wish to run.\n";
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

// Here is our hack for the fact that we can't have vectors of mpc_t
// and we can't pass variable multidimensional arrays to functions:
// a one-dimensional array of length equal to the number of lambdas
// times (p-1). Yuk.
void fillV(const unsigned int lambdaLength,
	   unsigned long p,
	   mpc_t &zeta,
	   const unsigned int lambdas[],
	   mpc_t evalV[]) {
  mpc_t primZetaEval[p];
  for (unsigned long n = 0; n < p; ++n) {
    // initialize the mpc_t objects in our array
    mpc_init2(primZetaEval[n], 64);
    // This 'evaluates' psi
    mpc_pow_ui(primZetaEval[n], zeta, n, MPFR_RNDN);
  }
  // There are p elements in the multiplicative group of characters,
  // but we aren't interested in evaluating the trivial character.
  // We are using our horrible hack of collapsing what should be a
  // 2D array into a 1D array, so each lambda has (p-1) characters
  // associated with it.
  for (unsigned long n = 0; n < (p-1)*lambdaLength; ++n) {
    mpc_init2(evalV[n], 64);
    mpc_set_si(evalV[n], 0, MPFR_RNDN);
  }

  // Remember: the multiplicative characters form a group of order p
  // zeta exponent needs mod p
  // lambdas account for a shift in the evaluation,
  // as does the choice of psi. Since the psis form
  // a multiplicative group, we only need to evaluate
  // a primitive root to get the data for all of them.
#pragma omp parallel for schedule(static) shared(lambdas, p, evalV, primZetaEval)
  for (unsigned long p1 = 0; p1 < p; ++p1) {
    ull logLookup;
    mpz_t psiArg, psiArgLambda, zetaPower;
    mpz_t p1big[6]; // to store powers of p1
    mpz_t p2big[4]; // to store powers of p2
    mpz_init(psiArg);
    mpz_init(psiArgLambda);
    mpz_init(zetaPower);
    mpz_init_set_ui(p1big[0], p1);
    for (int i = 1; i < 6; ++i) {
      mpz_init(p1big[i]);
      mpz_mul(p1big[i], p1big[0], p1big[i-1]);
    }
    for (int i = 0; i < 4; ++i)
      mpz_init(p2big[i]);
    for (unsigned long p2 = 0; p2 < p; ++p2) {
      mpz_set_ui(p2big[0], p2);
      for (int i = 1; i < 4; ++i)
	mpz_mul(p2big[i], p2big[0], p2big[i-1]);
      mpz_set_ui(psiArg, 0);
      // our polynomial is 1+y+xy+x^2y+x^3y+x^2y^2+lambdax^3y^2+x^4y^2+x^3y^3+x^4y^3+x^5y^3+x^6y^3+x^6y^4
      // the way the indices of the p1big and p2big arrays work, p1^i = p1big[i-1].
      mpz_add_ui(psiArg, psiArg, 1);
      mpz_add(psiArg, psiArg, p2big[0]);
      mpz_mod_ui(psiArg, psiArg, p);
      mpz_addmul(psiArg, p1big[0], p2big[0]);
      mpz_mod_ui(psiArg, psiArg, p);
      mpz_addmul(psiArg, p1big[1], p2big[0]);
      mpz_mod_ui(psiArg, psiArg, p);
      mpz_addmul(psiArg, p1big[2], p2big[0]);
      mpz_mod_ui(psiArg, psiArg, p);
      mpz_addmul(psiArg, p1big[1], p2big[1]);
      mpz_mod_ui(psiArg, psiArg, p);
      mpz_addmul(psiArg, p1big[3], p2big[1]);
      mpz_mod_ui(psiArg, psiArg, p);
      mpz_addmul(psiArg, p1big[2], p2big[2]);
      mpz_mod_ui(psiArg, psiArg, p);
      mpz_addmul(psiArg, p1big[3], p2big[2]);
      mpz_mod_ui(psiArg, psiArg, p);
      mpz_addmul(psiArg, p1big[4], p2big[2]);
      mpz_mod_ui(psiArg, psiArg, p);
      mpz_addmul(psiArg, p1big[5], p2big[2]);
      mpz_mod_ui(psiArg, psiArg, p);
      mpz_addmul(psiArg, p1big[5], p2big[3]);
      mpz_mod_ui(psiArg, psiArg, p);
      for (unsigned int l = 0; l < lambdaLength; ++l) {
	mpz_set_ui(psiArgLambda, 0);
	mpz_set_ui(zetaPower, 0);
	mpz_mul(psiArgLambda, p1big[2], p2big[1]);
	mpz_mul_ui(psiArgLambda, psiArgLambda, lambdas[l]);
	mpz_add(psiArgLambda, psiArg, psiArgLambda);
	mpz_mod_ui(psiArgLambda, psiArgLambda, p);
	// disregard the trivial character where c = 0 or p
        for (unsigned long c = 1; c < p; ++c) {
	  // Remember, additive characters form a multiplicative group of order p.
	  mpz_set(zetaPower, psiArgLambda);
	  mpz_mul_ui(zetaPower, zetaPower, c);
	  // Remember that zeta is a (p-1)th root of unity
	  mpz_mod_ui(zetaPower, zetaPower, p);
	  // We look up the evaluation of psi at this point.
	  // the primZetaEval array is actually canonically indexed; i.e.
	  // zeta^n is in the nth spot.	  
	  #pragma omp critical (summing)
	  {
	    mpc_add(evalV[(p-1)*l+(c-1)], evalV[(p-1)*l+(c-1)],
		    primZetaEval[mpz_get_ui(zetaPower)], MPFR_RNDN);
	  }
	}
      }
    }
    mpz_clear(psiArg);
    mpz_clear(psiArgLambda);
    mpz_clear(zetaPower);
    for (int i = 0; i < 6; ++i)
      mpz_clear(p1big[i]);
    for (int i = 0; i < 4; ++i)
      mpz_clear(p2big[i]);
  }

  // clean up
  for (unsigned long n = 0; n < p; ++n) {
    mpc_clear(primZetaEval[n]);
  }
  
  return;
}

void print2DV(const unsigned int lambdaLength, unsigned long p, mpc_t V[]) {
  FILE *output;
  char output_filename[ 32 ];
  sprintf (output_filename, "G2AddOut/g2add%06lu", p);
  output = fopen ( output_filename, "w" );

  fprintf(output, "{{%lu},{", p);
  // iterate on the lambdas
  for (unsigned int a = 0; a < lambdaLength; ++a) {
    fprintf(output, "{");
    for (unsigned long b = 0; b < (p-1); ++b) {
      mpfr_out_str(output, 10, 0, mpc_realref(V[a*(p-1)+b]), MPFR_RNDN);
      fprintf(output, "+(");
      mpfr_out_str(output, 10, 0, mpc_imagref(V[a*(p-1)+b]), MPFR_RNDN);
      fprintf(output, "I)");
      if (b < p-2) fprintf(output, ",");
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
  
  mpc_t zeta;
  mpc_init2(zeta, 64);
  assignZeta(zeta, p);

  mpc_t lambdaPsiV[lambdaLength*(p-1)];
  fillV(lambdaLength, p, zeta, lambdas, lambdaPsiV);
  print2DV(lambdaLength, p, lambdaPsiV);
  
  // cleanup
  mpc_clear(zeta);
  clear2DV(lambdaLength*(p-1), lambdaPsiV);
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

  // get name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  
  printf("Hello from processor %s, rank %d out of %d processors\n", processor_name, rank, nprocs);
  for (unsigned int n = params.first; n < params.first + params.second; ++n) {
    if ( n % nprocs == rank ) {
      pCharSum(n, lambdaLength, lambdas);
    }
  }

  // Sync the MPI processes and quit
  MPI_Finalize();
  return 0;
}
