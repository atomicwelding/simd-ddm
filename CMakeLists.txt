cmake_minimum_required(VERSION 3.11.2)
set(CMAKE_VERBOSE_MAKEFILE OFF)
project(simd-ddm CXX)

# generate json file for completion
set(CMAKE_EXPORT_COMPILE_COMMANDS 1) 

if(DEFINED ENV{MYROOT})
	set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{MYROOT})
	link_directories($ENV{MYROOT}/lib)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Necessary libraries
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/scripts)

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

find_package(Boost REQUIRED COMPONENTS program_options filesystem)
if(WIN32)
	find_package(FFTW3 CONFIG REQUIRED)
	find_package(FFTW3f CONFIG REQUIRED)
	set(FFTW_LIBRARIES ${FFTW3_LIBRARIES} ${FFTW3f_LIBRARIES})
	set(FFTW_INCLUDES ${FFTW3_INCLUDE_DIRS})
	link_directories(${FFTW3f_LIBRARY_DIRS})
	if(FFTW_LIBRARIES)
		message("-- Found FFTW3 libraries: ${FFTW_LIBRARIES}")
	endif()
else()
	find_package(FFTW REQUIRED)
endif()
find_package(TinyTIFFShared REQUIRED)
find_package(TIFF REQUIRED)
find_package(GSL REQUIRED)


###########################################################################
# Custom targets to ease the switching between debug or release compilation
###########################################################################

set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -g")
if(WIN32)
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fp=precise -march=cascadelake -DNDEBUG")
else()
	set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -mtune=native -DNDEBUG")
endif()
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
	${Boost_INCLUDE_DIRS}
	${FFTW_INCLUDES}
	${TIFF_INCLUDE_DIRS}
	${GSL_INCLUDE_DIRS}
	$ENV{MYROOT}/include
)


#################
# The main target
#################

add_executable(${PROJECT_NAME} src/main.cpp)
if(WIN32)
	set_property(
		TARGET ${PROJECT_NAME}
		PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
	)
else()
	set_property(
		TARGET ${PROJECT_NAME}
		PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
	)
endif()

#####################
# The library targets
#####################

add_library(Ddm ${DDMLIB_CPP})

###########################################
# We link the main target and the libraries
###########################################

target_link_libraries(${PROJECT_NAME} Ddm)
target_link_libraries(Ddm ${FFTW_LIBRARIES} ${Boost_LIBRARIES} ${TIFF_LIBRARIES} TinyTIFFShared GSL::gsl OpenMP::OpenMP_CXX)
