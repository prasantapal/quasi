PROJECT=quasi
CC=g++-9
CC+=-std=c++2a
CFLAGS =-g3
CFLAGS +=-Wall
CFLAGS +=-I/Users/prasantapal/softwares/lib/jsoncpp/include
LIBS += -ljsoncpp


${PROJECT}: ${PROJECT}.o
	${CC} ${CFLAGS}  -o ${PROJECT} ${PROJECT}.o ${LIBS}
${PROJECT}.o: ${PROJECT}.cpp 
	${CC} ${CFLAGS}  -c ${PROJECT}.cpp
clean:
	rm -rf *.o
