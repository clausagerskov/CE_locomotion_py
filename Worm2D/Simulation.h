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
        curvfile.open(rename_file("sim_curv.dat"));
        bodyfile.open(rename_file("sim_body.dat"));
        actfile.open(rename_file("sim_act.dat"));

    }
    void runSimulation(Worm2D & w);


    ~Simulation(){actfile.close(); curvfile.close(); bodyfile.close();}
    
private:
const double Duration, StepSize;
string rename_file(string filename);
string directoryName;
const int skip_steps;
ofstream actfile, curvfile, bodyfile;

};
