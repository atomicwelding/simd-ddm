CC = clang++
CFLAGS = -std=c++20 -Wall -L/usr/lib/x86_64-linux-gnu/
LDFLAGS = -lfftw3 -lboost_program_options

SRCS = main.cpp app.cpp utils/stack.cpp utils/frame.cpp
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
