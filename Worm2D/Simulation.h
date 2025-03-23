#include "Worm2D.h"

class Simulation
{
public:
    Simulation(string directoryName_):directoryName(directoryName_)
    {
        actfile.open(rename_file("act.dat"));
    }
    void runSimulation(Worm2D & w);


    ~Simulation(){actfile.close();}
    
private:
double Duration, StepSize;
string rename_file(string filename){return directoryName + "/" + filename;}
string directoryName;
int skip_steps;
ofstream actfile;

};
