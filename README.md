c
=
Compile and install the following libraries:
gmp-4.3.2 (ftp://ftp.gnu.org/gnu/gmp/gmp-4.3.2.tar.gz)
ntl-6.0.0 (http://www.shoup.net/ntl/ntl-6.0.0.tar.gz)
mpfr-3.1.2 (http://www.mpfr.org/mpfr-current/mpfr-3.1.2.tar.gz)
mpc-1.0.1 (http://www.multiprecision.org/mpc/download/mpc-1.0.1.tar.gz)

I have run into issues compiling NTL 6.0.0 with a version of GMP newer than 4.3.2.
We use MPI and OpenMP, so make sure they're installed in your environment as well.