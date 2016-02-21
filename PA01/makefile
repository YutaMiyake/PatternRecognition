CC  = g++ -std=c++11
CFLAGS = -c -g -Wall
SRCS = Generator.cpp BayesianClassifier.cpp driver.cpp
OBJS = ${SRCS:.cpp=.o}

TestDriver: $(OBJS)
	$(CC) -o $@ $(OBJS)

clean:
	-rm -f TestDriver *.o

.SUFFIXES: .cpp .c .o

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@