cmake_minimum_required (VERSION 3.16..3.24)
project (rng VERSION 0.0.1
  DESCRIPTION "Random number generator"
                    LANGUAGES CXX)

add_library (rng STATIC
    ${PROJECT_SOURCE_DIR}/src/random_number.cpp)

target_include_directories (rng PUBLIC
    ${PROJECT_SOURCE_DIR}/include)

