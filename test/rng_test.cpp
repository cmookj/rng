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

TEST (TestRandomInteger, DiscreteDistribution) {
  std::vector<double>                             prob{0.25, 0.3, 0.1, 0.2, 0.15};
  random_number::random_int_discrete_distribution rng{prob, 310952L};

  EXPECT_EQ (rng.tbl()[0].acceptance, 0.75);
  EXPECT_EQ (rng.tbl()[1].acceptance, 1.0);
  EXPECT_EQ (rng.tbl()[2].acceptance, 0.5);
  EXPECT_EQ (rng.tbl()[3].acceptance, 1.0);
  EXPECT_EQ (rng.tbl()[4].acceptance, 0.75);

  EXPECT_EQ (rng.tbl()[0].alias, 2);
  EXPECT_EQ (rng.tbl()[2].alias, 1);
  EXPECT_EQ (rng.tbl()[4].alias, 2);

  const int        max_count = 30000000;
  std::vector<int> histogram (prob.size(), 0);
  for (int i = 0; i < max_count; ++i) {
    int random_number = rng();
    histogram[random_number - 1]++;
  }

  for (int i = 0; i < histogram.size(); ++i) {
    EXPECT_NEAR (prob[i], double (histogram[i]) / double (max_count), 0.005);
  }
}
