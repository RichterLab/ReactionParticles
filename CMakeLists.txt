cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(ReactionParticle CXX)

# Default to Release
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Require C++ 11
set (CMAKE_CXX_STANDARD 11)

# Setup Include Directories
include_directories(".")

# CUDA
option( BUILD_CUDA "Build CUDA code" OFF)
if (BUILD_CUDA)
  find_package(CUDA REQUIRED)

  # Setup CUDA Compilation for Maximum compatability
  set( CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode arch=compute_20,code=compute_20 -gencode arch=compute_30,code=compute_30 -gencode arch=compute_35,code=compute_35 -gencode arch=compute_50,code=compute_50")
  set( CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -DBUILD_CUDA -std=c++11" )

  set_source_files_properties( "main.cpp" PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ )
  set_source_files_properties( "reaction.cpp" PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ )
  CUDA_ADD_EXECUTABLE( reaction-particle "main.cpp" "reaction.cpp")
else(BUILD_CUDA)
    add_executable (reaction-particle "main.cpp" "reaction.cpp")
endif(BUILD_CUDA)

# Build Executable
install ( TARGETS reaction-particle RUNTIME DESTINATION "bin")

option( BUILD_TESTS "Build tests" OFF)
if (BUILD_TESTS)
    # Download and unpack googletest at configure time
    configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" . WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest )
    execute_process(COMMAND ${CMAKE_COMMAND} --build . WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest )

    # Prevent GoogleTest from overriding our compiler/linker options
    # when building with Visual Studio
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    # Add googletest directly to our build. This adds
    # the following targets: gtest, gtest_main, gmock
    # and gmock_main
    add_subdirectory(test/google ${CMAKE_BINARY_DIR}/googletest-build)

    # The gtest/gmock targets carry header search path
    # dependencies automatically when using CMake 2.8.11 or
    # later. Otherwise we have to add them here ourselves.
    if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
    endif()

    include_directories("${CMAKE_SOURCE_DIR}" "test/")
    add_executable (reaction-particle-test "test/main.cpp" "test/reaction.cpp" "reaction.cpp")

    option( BUILD_TEST_COVERAGE "Enable code coverage calculation" OFF)
    if( BUILD_TEST_COVERAGE )
        set_target_properties(reaction-particle-test PROPERTIES COMPILE_FLAGS "-coverage" )
        set_target_properties(reaction-particle-test PROPERTIES LINK_FLAGS "-coverage" )
    endif( BUILD_TEST_COVERAGE )

    target_link_libraries (reaction-particle-test gtest)
endif (BUILD_TESTS)