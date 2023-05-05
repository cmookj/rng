cmake_minimum_required (VERSION 3.16..3.24)
project (my_project VERSION 0.0.1
                    DESCRIPTION "A test project"
                    LANGUAGES CXX)

add_library (my_library STATIC
    ${PROJECT_SOURCE_DIR}/src/my_src.cpp)

target_include_directories (my_library PUBLIC
    ${PROJECT_SOURCE_DIR}/include)

