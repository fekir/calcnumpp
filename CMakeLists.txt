cmake_minimum_required(VERSION 3.2)

project(calcnum CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(catch REQUIRED)
include_directories(SYSTEM ${CATCH_INCLUDE_DIRS})
include(compiler_settings)

# set standard, without extensions
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


include_directories("include")

add_executable(${PROJECT_NAME}Test
	main.cpp
	# generic includes
	include/utils.hpp
	include/statistics.hpp
	include/intervals.hpp

	# root finding algorithms
	include/bisec.hpp
	include/newton.hpp
	# --- test

	# test generic methods
	test/utils.cpp
	test/intervals.cpp

	# test root finding algorithms
	test/bisec.cpp
	test/newton.cpp
)
