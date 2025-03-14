main: main.o WormRS18.o Worm2D.o WormBody.o NervousSystem.o StretchReceptor18.o Muscles.o TSearch.o random.o argUtils.o jsonUtils.o
	g++ -pthread -o main main.o WormRS18.o Worm2D.o WormBody.o argUtils.o NervousSystem.o StretchReceptor18.o Muscles.o jsonUtils.o TSearch.o random.o 
random.o: random.cpp random.h VectorMatrix.h
	g++ -c -O3 -flto random.cpp
jsonUtils.o: jsonUtils.cpp jsonUtils.h
	g++ -c -O3 -std=c++11 -I/opt/homebrew/Cellar/nlohmann-json/3.11.3/include  -flto jsonUtils.cpp	
TSearch.o: TSearch.cpp TSearch.h
	g++ -c -O3 -flto TSearch.cpp
argUtils.o: ../argUtils.cpp ../argUtils.h
	g++ -c -O3 -flto ../argUtils.cpp
Worm2D.o: Worm2D.cpp Worm2D.h
	g++ -c -O3 -flto Worm2D.cpp
WormRS18.o: WormRS18.cpp WormRS18.h
	g++ -c -O3 -flto WormRS18.cpp
WormBody.o: WormBody.cpp WormBody.h
	g++ -c -O3 -flto WormBody.cpp
NervousSystem.o: NervousSystem.cpp NervousSystem.h VectorMatrix.h random.h
	g++ -c -O3 -flto NervousSystem.cpp
StretchReceptor18.o: StretchReceptor18.cpp StretchReceptor18.h
	g++ -c -O3 -flto StretchReceptor18.cpp
Muscles.o: Muscles.cpp Muscles.h VectorMatrix.h random.h
	g++ -c -O3 -flto Muscles.cpp
main.o: main.cpp WormRS18.h Worm2D.h WormBody.h StretchReceptor18.h Muscles.h TSearch.h
	g++ -c -O3 -flto main.cpp
clean:
	rm -f *.o main
