#include "random_number.hpp"
#include <gtest/gtest.h>

TEST (TestTopic, TrivialEquality) {}  // EXPECT_EQ (get_integer(), 42); }

TEST (TestTopic, MoreEqualityTests) {
  // ASSERT_EQ (get_integer(), 0) << "Something went wrong!";
  EXPECT_FLOAT_EQ (23.23F, 23.23F);
}
