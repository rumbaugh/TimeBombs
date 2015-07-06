CFLAGS = -O2 -Wall -Wextra -ansi -pedantic -I/home/rumbaugh/git/RJObject -I/usr/local/include
LIBS = -L/home/rumbaugh/git/RJObject -lrjobject -L/usr/local/lib -ldnest3 \
       -lgsl -lgslcblas \
       -lboost_system -lboost_thread-mt

default:
	g++ $(CFLAGS) -c *.cpp
	g++ -o main *.o $(LIBS)
	rm -f *.o
