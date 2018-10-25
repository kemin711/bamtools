# (c) 2017 Kemin Zhou
include_directories(SYSTEM /usr/local/include)
include_directories(${INSTALL_PREFIX}/include)
include_directories(${CMAKE_INSTALL_PREFIX}/include)
link_directories(${INSTALL_PREFIX}/lib)
link_directories(${CMAKE_INSTALL_PREFIX}/lib)

#TODO: this make file is not working with parallel make yet
# needs to spend one or two days to learn the intricacies of
# cmake

#find_package(PkgConfig)
#pkg_check_modules(GTEST REQUIRED gtest>=1.7.0)
#pkg_check_modules(GMOCK REQUIRED gmock>=1.7.0)

#include_directories(
#	${GTEST_INCLUDE_DIRS}
#   ${GMOCK_INCLUDE_DIRS}
#)
enable_testing()
add_executable(testbamalignment testbamalignment.cpp)
#add_dependencies(testbamalignment ${BamTools_SOURCE_DIR}/lib/libbamtools.a)
target_link_libraries(testbamalignment gtest gtest_main pthread ${BamTools_SOURCE_DIR}/lib/libbamtools.a)
add_test(testbamalignment, testbamalignment)
