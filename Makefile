all: Matrix.o mersenne.o MultiViewCoClust.o main.o
	g++ -g -O3 -o multiViewCoClu Matrix.o mersenne.o MultiViewCoClust.o main.o

Matrix.o: Matrix.cpp
			g++ -g -c Matrix.cpp

mersenne.o: mersenne.cpp
			g++ -g -c mersenne.cpp

MultiViewCoClust.o : MultiViewCoClust.cpp
			g++ -g -c MultiViewCoClust.cpp
main.o:	main.cpp
			g++ -g -c main.cpp

clean: multiViewCoClu *.o
		rm multiViewCoClu
		rm *.o