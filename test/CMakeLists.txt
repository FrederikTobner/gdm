cmake_minimum_required(VERSION 3.24)

enable_testing()


include(../cmake/CodeCoverage.cmake)
append_coverage_compiler_flags()

set(GDM_TESTS ${PROJECT_NAME}-Tests)

#Fetches google-test framework from the github repo
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.14.0
)

set(SOURCE_PATH ${PROJECT_SOURCE_DIR}/src)

file(GLOB_RECURSE TEST_SOURCES *.cpp *.hpp)

# Builds tests
add_executable(${GDM_TESTS} 
        ${TEST_SOURCES} ${GDM_SOURCES})

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

#Fetches google test
FetchContent_MakeAvailable(googletest)

#Links googletest and googlemock to the test executable
target_link_libraries(${GDM_TESTS} GTest::gtest_main GTest::gmock_main)

include(GoogleTest)

gtest_discover_tests(${GDM_TESTS})
