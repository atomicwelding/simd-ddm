cmake_minimum_required(VERSION 3.11.2)
set(CMAKE_VERBOSE_MAKEFILE OFF)
project(ddm CXX)

# generate json file for completion
set(CMAKE_EXPORT_COMPILE_COMMANDS 1) 

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Necessary libraries
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/scripts)

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

find_package(Boost REQUIRED COMPONENTS program_options)
find_package(FFTW REQUIRED)


###########################################################################
# Custom targets to ease the switching between debug or release compilation
###########################################################################

set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -g")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -mtune=native -DNDEBUG")
add_custom_target(debug
	COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
	COMMAND +${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
	COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
)
add_custom_target(release
	COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
	COMMAND +${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
	COMMENT "Switch CMAKE_BUILD_TYPE to Release"
)


#############################################
# Custom target for showing the project stats
#############################################

file(GLOB_RECURSE PROJECTFILES src/*.cpp src/*.hpp)
file(GLOB_RECURSE EXCLUDE json.h)
list(REMOVE_ITEM PROJECTFILES ${EXCLUDE})

add_custom_target(stats
	COMMAND cloc ${PROJECTFILES}
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src
	COMMENT "Statistics of the project:"
)


#################################################
# Definition of library directories and cpp files
#################################################

set(DDMLIB_DIR src/includes)
file(GLOB DDMLIB_CPP ${DDMLIB_DIR}/*.cpp)

INCLUDE_DIRECTORIES(
	${DDMLIB_DIR}
)


#################
# The main target
#################

add_executable(${PROJECT_NAME} src/main.cpp)
set_property(
	TARGET ${PROJECT_NAME}
	PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
)
add_custom_target(run
	COMMAND ./${PROJECT_NAME}
	DEPENDS ${PROJECT_NAME}
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	COMMENT "Run with ${CMAKE_BUILD_TYPE} configuration"
)


#####################
# The library targets
#####################

add_library(Ddm ${DDMLIB_CPP})

###########################################
# We link the main target and the libraries
###########################################

target_link_libraries(${PROJECT_NAME} Ddm)
target_link_libraries(Ddm ${FFTW_LIBRARIES} ${Boost_LIBRARIES})
