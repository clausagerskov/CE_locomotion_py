#include "TSearch.h"
#include "VectorMatrix.h"
#include "Evolution.h"

class EvolutionRS18:public Evolution
{
    public:
    //EvolutionRS18(const SuppliedArgs2018 & sa)
    //:Evolution(getEvoPars(sa),30),speedoutput(sa.speedoutput),evo_seed(sa.evo_seed){}
    
    /* EvolutionRS18(int argc, const char* argv[])
    :Evolution(argc,argv),speedoutput(getParameter<bool>(argc,argv,"--speed_output", 0)),
    evo_seed(getParameter<bool>(argc,argv,"--evo_seed", 0)){} */

    EvolutionRS18(int argc, const char* argv[])
    :Evolution(argc,argv),speedoutput(atoi(getParameter(argc,argv,"--speed_output", "0"))),
    evo_seed(atoi(getParameter(argc,argv,"--evo_seed", "0"))){}

    void GenPhenMapping(TVector<double> &gen, TVector<double> &phen);
    double EvaluationFunction(TVector<double> &v, RandomState &rs);
    double EvaluationFunctionOrig(TVector<double> &v, RandomState &rs);
    double EvaluationFunctionNoOut(TVector<double> &v, RandomState &rs);
    void RunSimulation(TVector<double> &v, RandomState &rs);
    void configure();
    //SuppliedArgs2018& supArgs1 = static_cast<SuppliedArgs2018&>(*supArgs1_ptr);
    //SuppliedArgs2018* const makeArgsPtr() {cout << "set pointer" << endl; return new SuppliedArgs2018();}
    //TSearch* const makeTSearchPtr() {return new TSearch(VectSize);}

    //SuppliedArgs2018& supArgs1;
    //evoPars getEvoPars(const SuppliedArgs2018 & sa);
    evoPars getDefaultPars();
    //const int VectSize = 30;
    private:

    //double Duration = 50.0;           // Seconds
    //double Transient = 10.0;          // 
    //double StepSize = 0.01;
    //int N_curvs = 23;                 // Number of cuvature points

// Used for Dumping: Frame rate for recording datais set to 50 frames per second

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