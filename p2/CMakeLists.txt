project(poisson)

cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0012 NEW)

enable_language(C)
enable_language(Fortran)

# CMake Options to enable/disable MPI and OpenMP
option(ENABLE_OPENMP "Enable OpenMP support?" ON)
option(ENABLE_MPI    "Enable MPI support?" ON)

if(ENABLE_MPI)
  find_package(MPI)
endif()

if(MPI_FOUND)
  add_definitions(-DHAVE_MPI=1)
  set(INCLUDES ${INCLUDES} ${MPI_INCLUDE_PATH})
  set(DEPLIBS_C ${DEPLIBS_C} ${MPI_C_LIBRARIES})
endif()

if(ENABLE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    add_definitions(-DHAVE_OPENMP=1)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  endif()
endif()

# YouCompleteMe stuff
if(EXISTS "${CMAKE_SOURCE_DIR}/.ycm")
  set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")
  if(EXISTS "${CMAKE_BINARY_DIR}/compile_commands.json")
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_if_different
      ${CMAKE_BINARY_DIR}/compile_commands.json
	  ${CMAKE_SOURCE_DIR}/compile_commands.json
    )
  endif()
endif()

# Use C99
set(CMAKE_C_FLAGS "-std=c99 ${CMAKE_C_FLAGS}")

# Send the include directories to the compiler
include_directories(${INCLUDES})

# Add a library with common sources
add_library(common fst.f)

# The C based poisson application
add_executable(poisson poisson.c list.c transpose.c)

# Link to libraries
target_link_libraries(poisson common ${DEPLIBS_C})

# We want to link this as a C program
set_target_properties(poisson PROPERTIES LINKER_LANGUAGE C)
