#include "random_number.hpp"

#include <algorithm>
#include <vector>

#include <gtest/gtest.h>

constexpr size_t dim = 2009;

TEST (TestRandomInteger, Uniform) {
  std::vector<int32_t> aa_int (dim);

  random_number::random_int_uniform rng_int{310952L};
  for (int m = 0; m < 2009; ++m)
    rng_int.ran_array (aa_int.data(), 1009);

  EXPECT_EQ (rng_int.ran_x (0), 995235265);

  rng_int.ran_start (310952L);
  for (int m = 0; m < 1009; ++m)
    rng_int.ran_array (aa_int.data(), 2009);

  EXPECT_EQ (rng_int.ran_x (0), 995235265);
}

TEST (TestRandomDouble, Uniform) {
  std::vector<double>                  aa_dbl (dim);
  random_number::random_double_uniform rng_dbl{310952L};
  for (int m = 0; m < 2009; ++m)
    rng_dbl.ranf_array (aa_dbl.data(), 1009);

  EXPECT_FLOAT_EQ (rng_dbl.ran_u (0), 0.36410514377569680455);

  rng_dbl.ranf_start (310952L);
  for (int m = 0; m < 1009; ++m)
    rng_dbl.ranf_array (aa_dbl.data(), 2009);

  EXPECT_FLOAT_EQ (rng_dbl.ran_u (0), 0.36410514377569680455);
}
