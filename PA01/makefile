CC  = g++
CFLAGS = -c -g -Wall
SRCS = generator.cpp
OBJS = ${SRCS:.cpp=.o}

Generator: $(OBJS)
	$(CC) -o $@ $(OBJS)

clean:
	-rm -f Generator *.o

.SUFFIXES: .cpp .c .o

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@