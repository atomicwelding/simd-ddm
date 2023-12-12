CC = clang++
CFLAGS = -std=c++20 -Wall -L/usr/lib/x86_64-linux-gnu/ -L/usr/local/lib -fopenmp=libomp -Qunused-arguments -O3 -mavx2 -Wmisleading-indentation
LDFLAGS = -fopenmp -lfftw3f_omp -lfftw3 -lfftw3f -lm -lboost_program_options -lTinyTIFFShared_Release -lgsl

SRCS = src/main.cpp src/includes/utils.cpp src/includes/app.cpp src/includes/stack.cpp src/includes/timer.cpp src/includes/curve_fit.cpp src/includes/ddm.cpp src/includes/delays.cpp src/includes/fit.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = main

.PHONY: all clean

all: $(TARGET) clean

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET) $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS)
