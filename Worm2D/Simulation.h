#include "Worm2D.h"
#include "Evolution.h"

struct simPars{
    string directoryName;
    int skip_steps;
    double Duration;       //
    double Transient;       //
    double StepSize;
};


class Simulation
{
public:

    Simulation(const simPars & sp1):sp(sp1)
    {
        curvfile.open(rename_file("sim_curv.dat"));
        bodyfile.open(rename_file("sim_body.dat"));
        actfile.open(rename_file("sim_act.dat"));
        velfile.open(rename_file("sim_vel.dat"));

    }
    void runSimulation(Worm2D & w);


    ~Simulation(){actfile.close(); curvfile.close(); bodyfile.close(); velfile.close();}
    
private:
const simPars sp;
string rename_file(string filename);
ofstream actfile, curvfile, bodyfile, velfile;

};
