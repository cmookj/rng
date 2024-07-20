//
// Random Number Generation Library
//
// Changmook Chun. (c) 2003--2023.
//
// Parts of this library is written by D. E. Knuth.
//

#ifndef __RANDOM_H_
#define __RANDOM_H_

#include <memory>

namespace random_number {

class random_int_uniform {
public:
    random_int_uniform(int32_t);
    virtual ~random_int_uniform() = default;
    
    random_int_uniform(const random_int_uniform&) = delete;
    random_int_uniform(random_int_uniform&&) = default;
    
    random_int_uniform& operator= (const random_int_uniform&) = delete;
    random_int_uniform& operator= (random_int_uniform&&) = default;
    
    // int operator() () const;
    
    // Do this before using _ran_array
    // Selector for different streams
    void ran_start (int32_t seed);

    // Put n new random numbers in array
    void ran_array (int32_t*, size_t n);
    const int32_t& ran_x(size_t) const;

    int ran_arr_next() {
        return *_ran_arr_ptr >= 0 ? *_ran_arr_ptr++: _ran_arr_cycle();
    }

protected:
    static constexpr int KK = 100;        // the long lag
    static constexpr int LL = 37;         // the short lag
    static constexpr int MM = (1L << 30); // the modulus
    static constexpr int QUALITY = 1009;  // recommended quality level for high-res use
    static constexpr int TT = 70;         // guaranteed separation between streams
    
    int32_t _ran_arr_dummy, _ran_arr_started;
    int32_t* _ran_arr_ptr; // the next random number, or -1

    std::unique_ptr<int32_t[]> _ran_arr_buf;
    std::unique_ptr<int32_t[]> _ran_x;

    int _ran_arr_cycle ();
    
    // Subtraction mod MM
    int _mod_diff (const int32_t x, const int32_t y) { return (x - y) & (MM - 1); }

    // Units bit of x
    bool _is_odd (const int32_t x) { return x & 1; }
};

class random_double_uniform {
public:
    random_double_uniform(int32_t);
    virtual ~random_double_uniform() = default;
    
    random_double_uniform(const random_double_uniform&) = delete;
    random_double_uniform(random_double_uniform&&) = default;
    
    random_double_uniform& operator= (const random_double_uniform&) = delete;
    random_double_uniform& operator= (random_double_uniform&&) = default;
    
    // double operator() () const;

    // Do this before using _ranf_array
    // Selector for different streams
    void ranf_start (int32_t seed);

    // Put n new random fractions in aa
    void ranf_array (double* aa, size_t n);
    const double& ran_u(size_t) const;
    
    double ranf_arr_next() {
        return *_ranf_arr_ptr >= 0 ? *_ranf_arr_ptr++ : _ranf_arr_cycle();
    }

protected:
    static constexpr int KK = 100;        // the long lag
    static constexpr int LL = 37;         // the short lag
    static constexpr int QUALITY = 1009;  // recommended quality level for high-res use
    static constexpr int TT = 70;         // guaranteed separation between streams
    
    std::unique_ptr<double[]> _ran_u = std::make_unique<double[]> (KK); // the generator state
    std::unique_ptr<double[]> _ranf_arr_buf = std::make_unique<double[]> (QUALITY);
    double _ranf_arr_dummy = -1.0, _ranf_arr_started=-1.0;
    double* _ranf_arr_ptr = &_ranf_arr_dummy; // the next random fraction, or -1

    // (x+y) mod 1.0
    double _mod_sum (const double x, const double y) {
        return (x + y) - static_cast<int> (x + y);
    }

    double _ranf_arr_cycle ();
    
    // Units bit of x
    bool _is_odd (const int32_t x) { return x & 1; }
};

//// normal dist.
//int dist_normal_polar_rejection (double*, int n);
//
//// gamma dist.
//int dist_gamma (double*, int n, double a, double r);
//
//// Brownian bridge
//int brownian_bridge_bisection (
//    double*, double*, double z_0, double z_N, int N
//);

} // namespace RandomNumber

#endif // __RANDOM_H_
