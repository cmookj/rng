//
// Random Number Generation Library
//
// Parts of this library is written by D. E. Knuth.
//

#ifndef _RANDOM_NUMBER_GENERATOR_H_
#define _RANDOM_NUMBER_GENERATOR_H_

#include <cmath>
#include <memory>
#include <vector>

namespace numeric {

template <typename T>
bool
close_enough (const T left, const T right, int max_ulps_diff = 1) {
  bool _close_enough = false;

  T lhs = std::min (left, right);
  T rhs = std::max (left, right);

  do {
    if (lhs >= rhs) {
      _close_enough = true;
      break;
    }
    lhs = std::nextafter (lhs, INFINITY);
  } while (max_ulps_diff--);

  return _close_enough;
}

}  // namespace numeric

namespace random_number {

// -----------------------------------------------------------------------------
//                                              Uniform random integer generator
// -----------------------------------------------------------------------------
class random_int_uniform final {
public:
  [[nodiscard]] random_int_uniform (int32_t);
  virtual ~random_int_uniform() = default;

  random_int_uniform (const random_int_uniform&) = delete;
  random_int_uniform (random_int_uniform&&)      = default;

  random_int_uniform&
  operator= (const random_int_uniform&) = delete;
  random_int_uniform&
  operator= (random_int_uniform&&) = default;

  // int operator() () const;

  // Do this before using _ran_array
  // Selector for different streams
  void
  ran_start (int32_t seed);

  // Put n new random numbers in array
  void
                 ran_array (int32_t*, size_t n);
  const int32_t& ran_x (size_t) const;

  int
  ran_arr_next () {
    return *_ran_arr_ptr >= 0 ? *_ran_arr_ptr++ : _ran_arr_cycle();
  }

private:
  static constexpr int KK      = 100;         // the long lag
  static constexpr int LL      = 37;          // the short lag
  static constexpr int MM      = (1L << 30);  // the modulus
  static constexpr int QUALITY = 1009;        // recommended quality level for high-res use
  static constexpr int TT      = 70;          // guaranteed separation between streams

  int32_t  _ran_arr_dummy, _ran_arr_started;
  int32_t* _ran_arr_ptr;  // the next random number, or -1

  std::unique_ptr<int32_t[]> _ran_arr_buf;
  std::unique_ptr<int32_t[]> _ran_x;

  int
  _ran_arr_cycle ();

  // Subtraction mod MM
  int
  _mod_diff (const int32_t x, const int32_t y) {
    return (x - y) & (MM - 1);
  }

  // Units bit of x
  bool
  _is_odd (const int32_t x) {
    return x & 1;
  }
};

// -----------------------------------------------------------------------------
//                                               Uniform random double generator
// -----------------------------------------------------------------------------
class random_double_uniform final {
public:
  [[nodiscard]] random_double_uniform (int32_t);
  virtual ~random_double_uniform() = default;

  random_double_uniform (const random_double_uniform&) = delete;
  random_double_uniform (random_double_uniform&&)      = default;

  random_double_uniform&
  operator= (const random_double_uniform&) = delete;
  random_double_uniform&
  operator= (random_double_uniform&&) = default;

  // double operator() () const;

  // Do this before using _ranf_array
  // Selector for different streams
  void
  ranf_start (int32_t seed);

  // Put n new random fractions in aa
  void
                ranf_array (double* aa, size_t n);
  const double& ran_u (size_t) const;

  double
  ranf_arr_next () {
    return *_ranf_arr_ptr >= 0 ? *_ranf_arr_ptr++ : _ranf_arr_cycle();
  }

private:
  static constexpr int KK      = 100;   // the long lag
  static constexpr int LL      = 37;    // the short lag
  static constexpr int QUALITY = 1009;  // recommended quality level for high-res use
  static constexpr int TT      = 70;    // guaranteed separation between streams

  std::unique_ptr<double[]> _ran_u        = std::make_unique<double[]> (KK);  // the generator state
  std::unique_ptr<double[]> _ranf_arr_buf = std::make_unique<double[]> (QUALITY);
  double                    _ranf_arr_dummy = -1.0, _ranf_arr_started = -1.0;
  double*                   _ranf_arr_ptr = &_ranf_arr_dummy;  // the next random fraction, or -1

  // (x+y) mod 1.0
  double
  _mod_sum (const double x, const double y) {
    return (x + y) - static_cast<int> (x + y);
  }

  double
  _ranf_arr_cycle ();

  // Units bit of x
  bool
  _is_odd (const int32_t x) {
    return x & 1;
  }
};

// -----------------------------------------------------------------------------
//                           Random integer generator with discrete distribution
// -----------------------------------------------------------------------------
class random_int_discrete_distribution final {
public:
  [[nodiscard]] random_int_discrete_distribution (
      const std::vector<double>& prob,
      const int32_t              seed
  );
  virtual ~random_int_discrete_distribution() = default;

  random_int_discrete_distribution (const random_int_discrete_distribution&) = delete;
  random_int_discrete_distribution (random_int_discrete_distribution&&)      = default;

  random_int_discrete_distribution&
  operator= (const random_int_discrete_distribution&) = delete;
  random_int_discrete_distribution&
  operator= (random_int_discrete_distribution&&) = delete;

  int
  operator() () const;

private:
  void
  generate_table (const std::vector<double>& prob);

  void
  print_table () const;

  const int max_rnd_number;
  enum class level { overfull, underfull, exactly_full };
  struct table_entry {
    int    i;
    double acceptance;  // U
    int    alias;       // K
    level  fill_level;
  };
  std::vector<table_entry> table;

  mutable random_double_uniform random_uniform_double_generator;

  // -------------------------------------
  //  Appendix: for debugging purpose only
  // -------------------------------------
public:
  const std::vector<table_entry>&
  tbl () const;
};

}  // namespace random_number

#endif
