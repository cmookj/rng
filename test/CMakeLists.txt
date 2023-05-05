# BUILD_TESTING variable is created by include(CTest)
# It is set to ON by default
if (BUILD_TESTING)
    add_executable (my_test my_test.cpp)
    target_link_libraries (my_test PRIVATE 
        my_library
        GTest::gtest_main)

    target_include_directories (my_test PRIVATE 
        ${PROJECT_SOURCE_DIR}/my_project/include)

    include (GoogleTest)
    # Finds all the Google tests associated with the executable
    gtest_discover_tests (my_test)
endif()

