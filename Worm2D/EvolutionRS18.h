#include "TSearch.h"
#include "VectorMatrix.h"
#include "Evolution.h"


class EvolutionRS18:public Evolution
{
    public:
    

    EvolutionRS18(int argc, const char* argv[])
    :Evolution(argc,argv, {".", 42, RANK_BASED, GENETIC_ALGORITHM, 
        96, 1000, 0.1, 0.5, UNIFORM, 
        1.1, 0.04, 1, 0, 1, 0, 50.0, 10.0, 0.01, 23, 30}, 30
    ),speedoutput(atoi(getParameter(argc,argv,"--speed_output", "0"))),
    evo_seed(atoi(getParameter(argc,argv,"--evo_seed", "0"))){}

    void GenPhenMapping(TVector<double> &gen, TVector<double> &phen);
    double EvaluationFunction(TVector<double> &v, RandomState &rs);
    double EvaluationFunctionOrig(TVector<double> &v, RandomState &rs);
    double EvaluationFunctionNoOut(TVector<double> &v, RandomState &rs);
    void RunSimulation(TVector<double> &v, RandomState &rs);
    void configure();
    
    protected:
    void addExtraParsToJson(json & j);

    
    private:


    const double fps = 25.0;
    const int skip = (int) (1/(evoPars1.StepSize*fps));

// Genotype -> Phenotype Mapping (Ventral cord)
    const double	  BiasRange			         	= 15.0;
    const double    SCRange                 = 15.0;
    const double    CSRange                 = 15.0;
    const double    TauMin                  = 0.5; //
    const double    TauMax                  = 2.0;

    const double    ESRange                 = 2.0;

    const double    SRmax                   = 200.0;
    const double    NMJmax                  = 1.0;

// (Head)
const double    HCSRange                = 15.0;

// Fitness
const double    AvgSpeed = 0.00022;              // Average speed of the worm in meters per seconds
const double    BBCfit = AvgSpeed*evoPars1.Duration;

// Size of genotype (VC)
    const bool speedoutput;
    const bool evo_seed;

};