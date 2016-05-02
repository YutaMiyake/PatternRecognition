CC  = g++ -std=c++11
CFLAGS = -c -g -Wall
SRCS = Driver.cpp Debugger.cpp libsvm-3.21/svm.cpp Matrix.cpp ML.cpp BayesianClassifier.cpp
OBJS = ${SRCS:.cpp=.o}

TestDriver: $(OBJS)
	$(CC) -o $@ $(OBJS)

clean:
	-rm -f TestDriver *.o

.SUFFIXES: .cpp .c .o

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@