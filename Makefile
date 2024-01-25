CC = clang++

SRC = src/main.cpp \
      $(wildcard src/includes/*.cpp)

OBJS = $(SRC:.cpp=.o)

CFLAGS = -std=c++20 \
         -Wall -O2 -Wmisleading-indentation -fopenmp=libomp

LDFLAGS = -fopenmp -lfftw3f_omp -lfftw3 -lfftw3f -lm -lboost_program_options -lTinyTIFFShared_Release -lgsl

# Platform specific
ARMFLAGS =  -I/opt/homebrew/include/ -I/usr/local/include/ \
            -L/opt/homebrew/lib/  -L/usr/local/lib/ \
            -rpath /usr/local/lib/

x86_64FLAGS = -mavx2 \
-L/usr/lib/x86_64-linux-gnu/ -L/usr/local/lib/

TARGET = armbuild 

.PHONY: all x86_64build armbuild  clean

all : $(TARGET) clean

x86_64build: $(OBJS)
	$(CC) $(CFLAGS) $(x86_64FLAGS) -o main $(OBJS) $(LDFLAGS)

armbuild: $(OBJS)
	$(CC) $(CFLAGS) $(ARMFLAGS) -o main $(OBJS) $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $(ARMFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS)

