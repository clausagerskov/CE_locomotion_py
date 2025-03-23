#include "Simulation.h"

string Simulation::rename_file(string filename){return directoryName + "/" + filename;}

void Simulation::runSimulation(Worm2D & w)
{
    
    for (double t = 0.0; t <= Duration; t += StepSize){
        w.Step(StepSize, 1);
        w.DumpActState(actfile, skip_steps);
    }
    
}
