#include "NervousSystemBase.h"
//#include "../NervousSystem.h"
#include "owSignalSimulatorForWorm2D.h"
#include "neuroml_utils.h"
//#include "../random.h"

class c302NervousSystem : public NervousSystemBase
{

public:

c302NervousSystem(const std::string & simFileName);
c302NervousSystem();
c302NervousSystem(const std::string & simFileName, const std::string & simDirName);

/* void setSimulator(
const std::string & simFileName = "neuromlLocal.main_sim", 
const std::string & simClassName = "Worm2DNRNSimulation", 
float timeStep=0.005);
void setSimulator(const std::string & simDirName,
const std::string & simFileName = "neuromlLocal.main_sim", 
const std::string & simClassName = "Worm2DNRNSimulation", 
float timeStep=0.005); */


//required for worm step

void SetNeuronExternalInput(int i, double value);
double NeuronOutput(int i);
void EulerStep(double );

//required for population structure in main_sim

void SetPopStructure(const std::string & popStruct, int popSize);

//required for nervous system setup

void SetChemicalSynapseWeight(int from, int to, double value);
void SetNeuronBias(int i, double value);
double NeuronBias(int i);
void SetNeuronGain(int i, double value);
void SetNeuronTimeConstant(int i, double value);
double NeuronTimeConstant(int i);
double NeuronState(int i);
double ChemicalSynapseWeight(int from, int to);


void SetCircuitSize(int newsize, int maxchemconns, int maxelecconns) ;
void SetNeuronOutput(int i, double value);
double ElectricalSynapseWeight(int from, int to);
void SetElectricalSynapseWeight(int n1, int n2, double value);


~c302NervousSystem();

// int CircuitSize(void)  {}
// void SetNeuronState(int i, double value) {}
 //void SetNeuronBias(int i, double value) {}
// double NeuronGain(int i) {}
// double NeuronExternalInput(int i) {}
 //void SetChemicalSynapseWeight(int from, int to, double value) {}
// void InternalSetElectricalSynapseWeight(int from, int to, double value) {}

    void RandomizeCircuitState(double lb, double ub);
    void RandomizeCircuitState(double lb, double ub, RandomState &rs);
    void RandomizeCircuitOutput(double lb, double ub);
    void RandomizeCircuitOutput(double lb, double ub, RandomState &rs);
 

 const std::vector<float> & getOutputValues() const {return output_value;}
 ostream & writeOutputValues(ostream & os) {return writeVector(os,output_value);}
 bool skipCalc = 1;

private:

SignalSimulatorForWorm2D *simulation;
std::vector<float> output_value;


};