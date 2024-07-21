//
// Random Number Generation Library
//
// Parts of this library is written by D. E. Knuth.
//

#include "random_number.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>
// #include <math.h>
// #include <stdlib.h>
// #include <time.h>
#include <numeric>
#include <type_traits>

namespace random_number {

// -----------------------------------------------------------------------------
//                                                            random_int_uniform
// -----------------------------------------------------------------------------
#pragma mark - Uniform Random Integer Generator

random_int_uniform::random_int_uniform (int32_t seed)
    : _ran_arr_dummy{-1}
    , _ran_arr_started{-1}
    , _ran_arr_ptr{&_ran_arr_dummy}
    , _ran_arr_buf{std::make_unique<int32_t[]> (QUALITY)}
    , _ran_x{std::make_unique<int32_t[]> (KK)} {
  ran_start (seed);
}

int
random_int_uniform::operator() () {
  return ran_arr_next();
}

void
random_int_uniform::ran_array (int32_t* aa, size_t n) {
  size_t i, j;
  for (j = 0; j < KK; j++)
    aa[j] = _ran_x[j];

  for (; j < n; j++)
    aa[j] = _mod_diff (aa[j - KK], aa[j - LL]);

  for (i = 0; i < LL; i++, j++)
    _ran_x[i] = _mod_diff (aa[j - KK], aa[j - LL]);

  for (; i < KK; i++, j++)
    _ran_x[i] = _mod_diff (aa[j - KK], _ran_x[i - LL]);
}

void
random_int_uniform::ran_start (int seed) {
  int     t, j;
  auto    x  = std::make_unique<int32_t[]> (KK + KK - 1);  // the preparation buffer
  int32_t ss = (seed + 2) & (MM - 2);

  for (j = 0; j < KK; j++) {
    x[j] = ss;  // bootstrap the buffer
    ss <<= 1;
    if (ss >= MM) ss -= MM - 2;  // cyclic shift 29 bits
  }

  x[1]++;  // make x[1] (and only x[1]) odd
  for (ss = seed & (MM - 1), t = TT - 1; t;) {
    for (j = KK - 1; j > 0; j--)
      x[j + j] = x[j], x[j + j - 1] = 0;  // "square"

    for (j = KK + KK - 2; j >= KK; j--)
      x[j - (KK - LL)]      = _mod_diff (x[j - (KK - LL)], x[j]),
                  x[j - KK] = _mod_diff (x[j - KK], x[j]);

    if (_is_odd (ss)) {  // "multiply by z"
      for (j = KK; j > 0; j--)
        x[j] = x[j - 1];

      x[0]  = x[KK];  // shift the buffer cyclically
      x[LL] = _mod_diff (x[LL], x[KK]);
    }
    if (ss) ss >>= 1;
    else t--;
  }

  for (j = 0; j < LL; j++)
    _ran_x[j + KK - LL] = x[j];

  for (; j < KK; j++)
    _ran_x[j - LL] = x[j];

  for (j = 0; j < 10; j++)
    ran_array (x.get(), KK + KK - 1);  // warm things up

  _ran_arr_ptr = &_ran_arr_started;
}

int
random_int_uniform::_ran_arr_cycle() {
  if (_ran_arr_ptr == &_ran_arr_dummy) ran_start (314159L);  // The user forgot to initialize

  ran_array (_ran_arr_buf.get(), QUALITY);
  _ran_arr_buf[KK] = -1;
  _ran_arr_ptr     = _ran_arr_buf.get() + 1;
  return _ran_arr_buf[0];
}

const int32_t&
random_int_uniform::ran_x (size_t i) const {
  return (i < KK ? _ran_x[i] : _ran_x[0]);
}

// -----------------------------------------------------------------------------
//                                                         random_double_uniform
// -----------------------------------------------------------------------------
#pragma mark - Uniform Random Fraction Generator

random_double_uniform::random_double_uniform (int32_t seed)
    : _ranf_arr_dummy{-1}
    , _ranf_arr_started{-1}
    , _ranf_arr_ptr{&_ranf_arr_dummy}
    , _ranf_arr_buf{std::make_unique<double[]> (QUALITY)}
    , _ran_u{std::make_unique<double[]> (KK)} {
  ranf_start (seed);
}

double
random_double_uniform::operator() () {
  return ranf_arr_next();
}

void
random_double_uniform::ranf_array (double* aa, size_t n) {
  // aa: destination
  // n : array length (must be at least KK)
  size_t i, j;
  for (j = 0; j < KK; j++)
    aa[j] = _ran_u[j];

  for (; j < n; j++)
    aa[j] = _mod_sum (aa[j - KK], aa[j - LL]);

  for (i = 0; i < LL; i++, j++)
    _ran_u[i] = _mod_sum (aa[j - KK], aa[j - LL]);

  for (; i < KK; i++, j++)
    _ran_u[i] = _mod_sum (aa[j - KK], _ran_u[i - LL]);
}

void
random_double_uniform::ranf_start (int32_t seed) {
  int              t, s, j;
  auto             u   = std::make_unique<double[]> (KK + KK - 1);
  constexpr double ulp = (1.0 / (1L << 30)) / (1L << 22);  // 2 to the -52
  double           ss  = 2.0 * ulp * ((seed & 0x3fffffff) + 2);

  for (j = 0; j < KK; j++) {
    u[j] = ss;  // bootstrap the buffer
    ss += ss;
    if (ss >= 1.0) ss -= 1.0 - 2 * ulp;  // cyclic shift of 51 bits
  }

  u[1] += ulp;  // make u[1] (and only u[1]) "odd"
  for (s = seed & 0x3fffffff, t = TT - 1; t;) {
    for (j = KK - 1; j > 0; j--)
      u[j + j] = u[j], u[j + j - 1] = 0.0;  // "square"

    for (j = KK + KK - 2; j >= KK; j--) {
      u[j - (KK - LL)] = _mod_sum (u[j - (KK - LL)], u[j]);
      u[j - KK]        = _mod_sum (u[j - KK], u[j]);
    }

    if (_is_odd (s)) {  // "multiply by z"
      for (j = KK; j > 0; j--)
        u[j] = u[j - 1];

      u[0]  = u[KK];  // shift the buffer cyclically
      u[LL] = _mod_sum (u[LL], u[KK]);
    }
    if (s) s >>= 1;
    else t--;
  }
  for (j = 0; j < LL; j++)
    _ran_u[j + KK - LL] = u[j];

  for (; j < KK; j++)
    _ran_u[j - LL] = u[j];

  for (j = 0; j < 10; j++)
    ranf_array (u.get(), KK + KK - 1);  // warm things up

  _ranf_arr_ptr = &_ranf_arr_started;
}

double
random_double_uniform::_ranf_arr_cycle() {
  if (_ranf_arr_ptr == &_ranf_arr_dummy) ranf_start (314159LL);  // the user forgot to initialize

  ranf_array (_ranf_arr_buf.get(), QUALITY);
  _ranf_arr_buf[KK] = -1;
  _ranf_arr_ptr     = _ranf_arr_buf.get() + 1;
  return _ranf_arr_buf[0];
}

const double&
random_double_uniform::ran_u (size_t i) const {
  return (i < KK ? _ran_u[i] : _ran_u[0]);
}

// -----------------------------------------------------------------------------
//                                              random_int_discrete_distribution
// -----------------------------------------------------------------------------
#pragma mark - Random Integer Generator with Discrete Distribution

void
random_int_discrete_distribution::print_table() const {
  for (const auto& entry : table) {
    std::cout << "i = " << entry.i << ", Ui = " << entry.acceptance << ", Ki = " << entry.alias
              << " - "
              << (entry.fill_level == level::exactly_full
                      ? "F"
                      : (entry.fill_level == level::overfull ? "O" : "U"))
              << '\n';
  }
}

random_int_discrete_distribution::random_int_discrete_distribution (
    const std::vector<double>& prob,
    const int32_t              seed
)
    : max_rnd_number{static_cast<int> (prob.size())}
    , random_uniform_double_generator{seed} {
  generate_table (prob);
}

void
random_int_discrete_distribution::generate_table (const std::vector<double>& prob) {
  // Initialize U(i) = n p(i)
  double sum = std::accumulate (prob.cbegin(), prob.cend(), 0.0);
  if (numeric::close_enough (sum, 0.0)) return;

  const int n = max_rnd_number;
  table.reserve (n);

  int i = 0;
  std::transform (
      prob.cbegin(),
      prob.cend(),
      std::back_inserter (table),
      [n, sum, &i] (const double p) {
        double u = n * p / sum;

        auto fill_level = level::exactly_full;
        if (u > 1.0) fill_level = level::overfull;
        if (u < 1.0) fill_level = level::underfull;

        i++;
        int k = (fill_level == level::exactly_full ? i : 0);

        return table_entry{i, u, k, fill_level};
      }
  );

  // Update the table
  for (int i = 0; i < n; ++i) {
    auto entry_i = std::find_if (table.begin(), table.end(), [] (const auto& e) {
      return e.fill_level == level::overfull;
    });

    auto entry_j = std::find_if (table.begin(), table.end(), [] (const auto& e) {
      return e.fill_level == level::underfull;
    });

    if (entry_j == table.end()) break;

    entry_j->alias = entry_i->i;
    entry_i->acceptance += (entry_j->acceptance - 1.);
    entry_j->fill_level = level::exactly_full;

    if (entry_i->acceptance < 1.0) entry_i->fill_level = level::underfull;

    if (numeric::close_enough (entry_i->acceptance, 1.0)) entry_i->fill_level = level::exactly_full;
  }
}

int
random_int_discrete_distribution::operator() () const {
  double       random_double = random_uniform_double_generator();
  const int    n             = max_rnd_number;
  const int    i             = std::floor (random_double * n) + 1;
  const double y             = random_double * n + 1 - i;

  if (y < table[i - 1].acceptance) return i;
  return table[i - 1].alias;
}

const std::vector<random_int_discrete_distribution::table_entry>&
random_int_discrete_distribution::tbl() const {
  return table;
}

// *****************************************************************************
//                                                              Utility routines
// *****************************************************************************

//// Forward Declarations of private functions
// void _brownian_bridge_bisect (int j, int J, double*, double*, double*);
//
// int dist_normal_polar_rejection (double* aa, int n) {
//     /*
//      This function generates n normally distributed random numbers.
//      The mean and standard deviation of the distribution is
//      0 and 1, respectively.
//
//      To see the detailed description of this algorithm,
//      consult D.E.Knuth's ``The Art of Computer Programming Volume 2''
//      section 3.4.1 on pp.119.
//
//      aa[]: array that will contain generated random sequence
//      n: int that tells the size of aa[]
//      */
//     int m, i, idx = 0;
//     double v1, v2, s, tmp;
//     int done = 0;
//
//     if (is_odd (n)) {
//         m = n + 1;
//     } else {
//         m = n;
//     }
//     auto _a = std::make_unique<double[]> (2 * m); // uniform deviates
//     // NOTE: _a is twice as large as _aa.
//     // see D.E.Knuth(pp.122) for the rationale.
//     auto _aa = std::make_unique<double[]> (m); // normal deviates
//
//     while (!done) {
//         // build an array of uniformly random numbers
//         // _ranf_start((int)time(NULL));
//         _ranf_start (static_cast<int32_t> (time (NULL)));
//         _ranf_array (_a.get(), 2 * m);
//
//         // calculate normally distributed numbers
//         // from the uniformly random numbers
//         for (i = 0; i < m; i++) {
//             v1 = _a[i * 2] * 2. - 1.;
//             v2 = _a[i * 2 + 1] * 2. - 1.;
//             s = v1 * v1 + v2 * v2;
//             if (s >= 1)
//                 continue;
//
//             tmp = sqrt (-2 * log (s) / s);
//             _aa[idx++] = v1 * tmp;
//             _aa[idx++] = v2 * tmp;
//
//             if (idx > m) {
//                 done = 1;
//                 break;
//             }
//         }
//     }
//     for (i = 0; i < n; i++) {
//         aa[i] = _aa[i];
//     }
//
//     return 0;
// }
//
// int dist_gamma (double* aa, int n, double a, double r) {
//     /*
//      int dist_gamma(double aa[], int n, double a, double r)
//      generates random deviates from gamma distribution
//
//      Function
//      Generates random deviates from the gamma distribution whose
//      density is
//
//      (A**R)/Gamma(R) * X**(R-1) * Exp(-A*X)
//
//      Arguments
//      a --> Location parameter of Gamma distribution
//      r --> Shape parameter of Gamma distribution
//
//      Method
//      Renames SGAMMA from TOMS as slightly modified by BWB to use RANF
//      instead of SUNIF.
//
//      For details see:
//      (Case R >= 1.0)
//      Ahrens, J.H. and Dieter, U.
//      Generating Gamma Variates by a
//      Modified Rejection Technique.
//      Comm. ACM, 25,1 (Jan. 1982), 47 - 54.
//
//      Algorithm GD
//      (Case 0.0 <= R <= 1.0)
//      Ahrens, J.H. and Dieter, U.
//      Computer Methods for Sampling from Gamma,
//      Beta, Poisson and Binomial Distributions.
//      Computing, 12 (1974), 223-246/
//      Adapted algorithm GS.
//      */
//
//     for (int i = 0; i < n; i++) {
//         aa[i] = gengam (a, r);
//     }
//
//     return 0;
// }
//
// int brownian_bridge_bisection (
//     double* tt,
//     double* zz,
//     double z_0,
//     double z_N,
//     int N
//) {
//     /*
//      Constructs a Brownian bridge from z_0 to z_N
//      this function uses bisection construction method that gives
//      ``better control'' over the generated sample paths.
//
//      tt[]: array of time history. from 0 to N, i.e., size = N+1.
//      zz[]: array of z_t (Wiener process). size = N+1.
//      z_0: starting point of z_t at tt[0].
//      z_N: destination of z_t at tt[N].
//      N: the number of steps in tt and zz.
//      */
//
//     // prepare data
//     // 1. build an array containing epsilon
//     auto eps = std::make_unique<double[]> (N + 1); // array of epsilon
//     _ranf_start ((int)time (NULL));
//     _ranf_array (eps.get(), N + 1); // normal, mean = 0, std = 1
//
//     // 2. put z_0 and z_N at both ends of array, zz
//     zz[0] = z_0;
//     zz[N] = z_N;
//
//     // Now, all the other elements of zz will be filled by
//     // recursive bisection routine.
//
//     // Call recursive bisection routine
//     _brownian_bridge_bisect (0, N, tt, zz, eps.get());
//
//     return 0;
// }
//
// void _brownian_bridge_bisect (
//     int j,
//     int J,
//     double* tt,
//     double* zz,
//     double* eps
//) {
//     /*
//      Recursive bisection routine called by brownian_bridge_bisection
//
//      int j;         left index, i.e., j- in the text
//      int J;         right index, i.e., j+ in the text
//      double* tt;     array of time
//      double* zz;     array of z
//      double* eps;     array of epsilon, standard normal random sequence
//      */
//
//     int jm;              // mid-way between j and J
//     double td, tld, trd; // time differences
//     // td: t_(j+) - t_(j-)
//     // tld: t_(jm) - t_(j-)
//     // trd: t_(j+) - t_(jm)
//
//     // Evaluate the condition to end the recursion
//     if (J == j + 1)
//         return;
//
//     // calculate index, jm, the point where bisection occurs
//     jm = (J - j) / 2 + j;
//
//     // calculate time differences
//     // this enables more efficient calculation of the bisection
//     td = tt[J] - tt[j];
//     tld = tt[jm] - tt[j];
//     trd = tt[J] - tt[jm];
//
//     // the result of the bisection
//     zz[jm] = (trd * zz[j] + tld * zz[J]) / td + sqrt (trd * tld / td) * eps[jm];
//
//     // now, call this bisection method (recursion) for the left part of jm
//     _brownian_bridge_bisect (j, jm, tt, zz, eps);
//
//     // now, call this bisection method (recursion) for the right part of jm
//     _brownian_bridge_bisect (jm, J, tt, zz, eps);
// }

}  // namespace random_number
