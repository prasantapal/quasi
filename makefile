PROJECT=particle
CC=g++
CC+=-std=c++17
CFLAGS=-g3
CFLAGS+=-Wall
CFLAGS+=-mmacosx-version-min=10.8


${PROJECT}: ${PROJECT}.o
	${CC} ${CFLAGS}  -o ${PROJECT} ${PROJECT}.o
${PROJECT}.o: ${PROJECT}.cpp
	${CC} ${CFLAGS}  -c ${PROJECT}.cpp
clean:
	rm -rf *.o
