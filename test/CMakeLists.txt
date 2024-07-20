# BUILD_TESTING variable is created by include(CTest)
# It is set to ON by default
if (BUILD_TESTING)
    add_executable (rng_test rng_test.cpp)
    target_link_libraries (rng_test PRIVATE 
        rng
        GTest::gtest_main)

    target_include_directories (rng_test PRIVATE 
        ${PROJECT_SOURCE_DIR}/rng/include)

    include (GoogleTest)
    # Finds all the Google tests associated with the executable
    gtest_discover_tests (rng_test)
endif()

