CC = g++ -Wall
PYTHON_CONFIG ?= python3-config

LIBS := $(shell $(PYTHON_CONFIG) --embed --libs)
LDFLAGS := $(shell $(PYTHON_CONFIG) --embed --ldflags)
CXXFLAGS := $(shell $(PYTHON_CONFIG) --cflags)

EXTRA_FLAGS :=   -fPIE -O3  -Iinc
main: main.o c302NervousSystem.o owSignalSimulator.o owSignalSimulatorForWorm2D.o WormRS18.o Worm21.o WormCE.o Segment21.o Worm2D.o WormBody.o Evolution21.o Evolution.o EvolutionCE.o EvolutionRS18.o NervousSystem.o StretchReceptor18.o StretchReceptorCE.o Muscles.o TSearch.o random.o jsonUtils.o
	$(CC) $(CXXFLAGS) $(LDFLAGS) $(EXTRA_FLAGS) -pthread -o main main.o c302NervousSystem.o owSignalSimulator.o owSignalSimulatorForWorm2D.o WormRS18.o Worm21.o WormCE.o Segment21.o Evolution21.o Evolution.o EvolutionCE.o EvolutionRS18.o Worm2D.o WormBody.o NervousSystem.o StretchReceptor18.o StretchReceptorCE.o Muscles.o jsonUtils.o TSearch.o random.o $(LIBS)
c302NervousSystem.o: ../neuromlLocal/c302NervousSystem.cpp ../neuromlLocal/owSignalSimulatorForWorm2D.h 
	$(CC) -c -O3 $(CXXFLAGS) $(LDFLAGS) $(EXTRA_FLAGS) ../neuromlLocal/c302NervousSystem.cpp 
	@echo ' '
owSignalSimulator.o: ../neuromlLocal/owSignalSimulator.cpp ../neuromlLocal/owSignalSimulator.h ../neuromlLocal/owINeuronSimulator.h
	$(CC) -c -O3 $(CXXFLAGS) $(LDFLAGS) $(EXTRA_FLAGS) ../neuromlLocal/owSignalSimulator.cpp
	@echo ' '
owSignalSimulatorForWorm2D.o: ../neuromlLocal/owSignalSimulatorForWorm2D.cpp ../neuromlLocal/owSignalSimulatorForWorm2D.h ../neuromlLocal/owSignalSimulator.h
	$(CC) -c -O3 $(CXXFLAGS) $(LDFLAGS) $(EXTRA_FLAGS) ../neuromlLocal/owSignalSimulatorForWorm2D.cpp
	@echo ' '
random.o: random.cpp random.h #VectorMatrix.h
	$(CC) -c -O3 -flto random.cpp
jsonUtils.o: jsonUtils.cpp jsonUtils.h
	$(CC) -c -O3 -std=c++11 -I/opt/homebrew/Cellar/nlohmann-json/3.11.3/include  -flto jsonUtils.cpp	
TSearch.o: TSearch.cpp TSearch.h
	$(CC) -c -O3 -flto TSearch.cpp
#argUtils.o: ../argUtils.cpp ../argUtils.h
#	g++ -c -O3 -std=c++11 -I/opt/homebrew/Cellar/nlohmann-json/3.11.3/include -flto ../argUtils.cpp
Worm2D.o: Worm2D.cpp Worm2D.h
	$(CC) -c -O3 -std=c++11 -I/opt/homebrew/Cellar/nlohmann-json/3.11.3/include -flto Worm2D.cpp
EvolutionRS18.o: EvolutionRS18.cpp EvolutionRS18.h
	$(CC) -c -O3 -std=c++11 -I/opt/homebrew/Cellar/nlohmann-json/3.11.3/include -flto EvolutionRS18.cpp
EvolutionCE.o: EvolutionCE.cpp EvolutionCE.h
	$(CC) -c -O3 -std=c++11 -I/opt/homebrew/Cellar/nlohmann-json/3.11.3/include -flto EvolutionCE.cpp
Evolution21.o: Evolution21.cpp Evolution21.h
	$(CC) -c -O3 -std=c++11 -I/opt/homebrew/Cellar/nlohmann-json/3.11.3/include -flto Evolution21.cpp
Evolution.o: Evolution.cpp Evolution.h
	$(CC) -c -O3 -std=c++11 -I/opt/homebrew/Cellar/nlohmann-json/3.11.3/include -flto Evolution.cpp
WormRS18.o: WormRS18.cpp WormRS18.h
	$(CC) -c -O3 -std=c++11 -I/opt/homebrew/Cellar/nlohmann-json/3.11.3/include -flto WormRS18.cpp
WormCE.o: WormCE.cpp WormCE.h
	$(CC) -c -O3 $(CXXFLAGS) $(LDFLAGS) $(EXTRA_FLAGS) -std=c++11 -I/opt/homebrew/Cellar/nlohmann-json/3.11.3/include -flto WormCE.cpp
Worm21.o: Worm21.cpp Worm21.h
	$(CC) -c -O3 -std=c++11 -I/opt/homebrew/Cellar/nlohmann-json/3.11.3/include -flto Worm21.cpp
Segment21.o: Segment21.cpp Segment21.h
	$(CC) -c -O3 -std=c++11 -flto Segment21.cpp	
WormBody.o: WormBody.cpp WormBody.h
	$(CC) -c -O3 -flto WormBody.cpp
NervousSystem.o: NervousSystem.cpp NervousSystem.h #VectorMatrix.h random.h ../neuromlLocal/NervousSystemBase.h
	$(CC) -c -O3 -std=c++11 -flto NervousSystem.cpp
StretchReceptor18.o: StretchReceptor18.cpp StretchReceptor18.h
	$(CC) -c -O3 -std=c++11 -flto StretchReceptor18.cpp
StretchReceptorCE.o: StretchReceptorCE.cpp StretchReceptorCE.h
	$(CC) -c -O3 -std=c++11 -flto StretchReceptorCE.cpp	
Muscles.o: Muscles.cpp Muscles.h #VectorMatrix.h random.h
	$(CC) -c -O3 -flto Muscles.cpp
main.o: main.cpp #WormRS18.h Worm2D.h WormBody.h StretchReceptor18.h Muscles.h TSearch.h
	$(CC) -c -O3 -std=c++11 -I/opt/homebrew/Cellar/nlohmann-json/3.11.3/include -flto main.cpp
clean:
	rm -f *.o main
