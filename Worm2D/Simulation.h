#include "Worm2D.h"
#include "Evolution.h"
class Simulation
{
public:
    Simulation(const evoPars & ep1):
    directoryName(ep1.directoryName),
    Duration(ep1.Duration),
    StepSize(ep1.StepSize),
    skip_steps(ep1.skip_steps)
    {
        actfile.open(rename_file("act.dat"));
    }
    void runSimulation(Worm2D & w);


    ~Simulation(){actfile.close();}
    
private:
const double Duration, StepSize;
string rename_file(string filename);
string directoryName;
const int skip_steps;
ofstream actfile;

};
