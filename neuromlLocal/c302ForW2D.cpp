#include "c302ForW2D.h"
//#include "owSignalSimulator.h"

//const bool skipCalc = 1;
const std::string defaultSimClassName = "Worm2DNRNSimulation";

c302ForW2D::c302ForW2D(const std::string & simFileName):
simulation(new SignalSimulatorForWorm2D(simFileName,defaultSimClassName,0.005))
{SetPopStructure("DA DB DD VD VA VB", 10);}

//c302NervousSystem::c302NervousSystem():
//simulation(new SignalSimulatorForWorm2D("neuromlLocal.main_sim",defaultSimClassName,0.005)){}

c302ForW2D::c302ForW2D():
simulation(new SignalSimulatorForWorm2D("main_sim",defaultSimClassName,"neuromlLocal",0.005))
{SetPopStructure("DA DB DD VD VA VB", 10);}

c302ForW2D::c302ForW2D(const std::string & simFileName, 
const std::string & simDirName):simulation(new SignalSimulatorForWorm2D(simFileName,
defaultSimClassName,simDirName,0.005)){SetPopStructure("DA DB DD VD VA VB", 10);}



void c302ForW2D::SetPopStructure(const std::string & popStruct, int popSize)
{
    simulation->strOneValFunc("set_up", popStruct, popSize);
}

void c302ForW2D::SetNeuronExternalInput(int i, double value)
{
simulation->oneValFunc("set_neuron_input",i-1,value);
}

double c302ForW2D::NeuronOutput(int i)
{
return output_value[i-1];
}

void c302ForW2D::EulerStep(double stepsize)
{       
    output_value = simulation->run();
}

