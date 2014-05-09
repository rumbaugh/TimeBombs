CFLAGS = -O2 -Wall -Wextra -ansi -pedantic -I${RJOBJECT_DIR}
LIBS =  -L${RJOBJECT_DIR} -lrjobject -ldnest3 -lgsl -lgslcblas \
	-lboost_system -lboost_thread

default:
	g++ $(CFLAGS) -c *.cpp
	g++ -o main *.o $(LIBS)
	rm -f *.o
