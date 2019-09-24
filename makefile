PROJECT=junction_test
CC=g++-9
CC+=-std=c++17
CFLAGS=-g3
CFLAGS+=-Wall


${PROJECT}: ${PROJECT}.o
	${CC} ${CFLAGS}  -o ${PROJECT} ${PROJECT}.o
${PROJECT}.o: ${PROJECT}.cpp junction.hpp
	${CC} ${CFLAGS}  -c ${PROJECT}.cpp
clean:
	rm -rf *.o
