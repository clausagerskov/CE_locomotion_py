#include "c302NervousSystem.h"
//#include "owSignalSimulator.h"

const bool skipCalc = 1;
const std::string defaultSimClassName = "Worm2DNRNSimulation";

c302NervousSystem::c302NervousSystem(const std::string & simFileName):
simulation(new SignalSimulatorForWorm2D(simFileName,defaultSimClassName,0.005)){}

//c302NervousSystem::c302NervousSystem():
//simulation(new SignalSimulatorForWorm2D("neuromlLocal.main_sim",defaultSimClassName,0.005)){}
c302NervousSystem::c302NervousSystem():
simulation(new SignalSimulatorForWorm2D("main_sim",defaultSimClassName,"neuromlLocal",0.005)){}

c302NervousSystem::c302NervousSystem(const std::string & simFileName, 
const std::string & simDirName):simulation(new SignalSimulatorForWorm2D(simFileName,
defaultSimClassName,simDirName,0.005)){}


/* void c302NervousSystem::setSimulator(const std::string & simFileName, 
const std::string & simClassName, float timeStep)
{   
    if (simulation!=nullptr) delete simulation;
    simulation = new SignalSimulatorForWorm2D(simFileName, simClassName, timeStep);
}

void c302NervousSystem::setSimulator(const std::string & simDirName,
const std::string & simFileName, const std::string & simClassName, 
float timeStep)
{   
    if (simulation!=nullptr) delete simulation;
    simulation = new SignalSimulatorForWorm2D(simFileName, simClassName, timeStep);
} */


void c302NervousSystem::SetNeuronExternalInput(int i, double value)
{
simulation->oneValFunc("set_neuron_input",i-1,value);
}

double c302NervousSystem::NeuronOutput(int i)
{
return output_value[i-1];
}

void c302NervousSystem::EulerStep(double stepsize)
{       
    output_value = simulation->run();
}

void c302NervousSystem::SetChemicalSynapseWeight(int from, int to, double value)
{
    if (skipCalc) return;
    simulation->twoValFunc("set_synaptic_weight",from-1,to-1,value);
}

void c302NervousSystem::SetNeuronBias(int i, double value)
{
    if (skipCalc) return;
    simulation->strValFunc("set_neuron_parameter", "bias", i-1, value);
}

void c302NervousSystem::SetNeuronGain(int i, double value)
{   
    if (skipCalc) return;
    simulation->strValFunc("set_neuron_parameter", "gain", i-1, value);
}

void c302NervousSystem::SetNeuronTimeConstant(int i, double value)
{   
    if (skipCalc) return;
    simulation->strValFunc("set_neuron_parameter", "tau", i-1, value);
}

double  c302NervousSystem::NeuronTimeConstant(int i)
{
    if (skipCalc) return 0;
    return simulation->getValFunc("get_neuron_parameter", "tau", i-1);
}

double c302NervousSystem::NeuronBias(int i)
{
    if (skipCalc) return 0;
    return simulation->getValFunc("get_neuron_parameter", "bias", i-1);
}

double c302NervousSystem::NeuronState(int i)
{
    if (skipCalc) return 0;
    return simulation->getValFunc("get_neuron_parameter", "state", i-1);
}
double c302NervousSystem::ChemicalSynapseWeight(int from, int to)
{
    if (skipCalc) return 0;
    return simulation->getTwoValFunc("get_synaptic_weight", from-1, to-1);

}
/* void c302NervousSystem::RandomizeCircuitState(double lb, double ub) {return;}
void c302NervousSystem::RandomizeCircuitState(double lb, double ub, RandomState &rs) {return;}
void c302NervousSystem::RandomizeCircuitOutput(double lb, double ub) {return;}
void c302NervousSystem::RandomizeCircuitOutput(double lb, double ub, RandomState &rs) {return;} */

void c302NervousSystem::SetCircuitSize(int newsize, int maxchemconns, int maxelecconns) {return;}
void c302NervousSystem::SetNeuronOutput(int i, double value) {return;}
double c302NervousSystem::ElectricalSynapseWeight(int from, int to) {return 0;}
void c302NervousSystem::SetElectricalSynapseWeight(int n1, int n2, double value) {return;} 

c302NervousSystem::~c302NervousSystem()
{if (simulation!=nullptr) delete simulation;}