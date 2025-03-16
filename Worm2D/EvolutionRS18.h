#include "TSearch.h"
#include "VectorMatrix.h"
#include "Evolution.h"

class EvolutionRS18:public Evolution
{
    public:
    EvolutionRS18():Evolution(new SuppliedArgs2018(),new TSearch(30)),
    supArgs1(static_cast<SuppliedArgs2018&>(*supArgs1_ptr))
    {
        evoPars1 = {supArgs1.randomseed, RANK_BASED, GENETIC_ALGORITHM, 
            supArgs1.pop_size, supArgs1.max_gens, 0.1, 0.5, UNIFORM, 
            1.1, 0.04, 1, 0, 1};

        cout << "const called " << s->VectorSize() << endl;
    
    }
    void GenPhenMapping(TVector<double> &gen, TVector<double> &phen);
    double EvaluationFunction(TVector<double> &v, RandomState &rs);
    double EvaluationFunctionOrig(TVector<double> &v, RandomState &rs);
    double EvaluationFunctionNoOut(TVector<double> &v, RandomState &rs);
    void configure();
    //SuppliedArgs2018& supArgs1 = static_cast<SuppliedArgs2018&>(*supArgs1_ptr);
    //SuppliedArgs2018* const makeArgsPtr() {cout << "set pointer" << endl; return new SuppliedArgs2018();}
    //TSearch* const makeTSearchPtr() {return new TSearch(VectSize);}

    SuppliedArgs2018& supArgs1;

    const int VectSize = 30;
    private:

    double Duration = 50.0;           // Seconds
    double Transient = 10.0;          // 
    double StepSize = 0.01;
    int N_curvs = 23;                 // Number of cuvature points

// Used for Dumping: Frame rate for recording datais set to 50 frames per second
    double fps = 25.0;
    int skip = (int) (1/(StepSize*fps));

// Genotype -> Phenotype Mapping (Ventral cord)
    double	  BiasRange			         	= 15.0;
    double    SCRange                 = 15.0;
    double    CSRange                 = 15.0;
    double    TauMin                  = 0.5; //
    double    TauMax                  = 2.0;

    double    ESRange                 = 2.0;

    double    SRmax                   = 200.0;
    double    NMJmax                  = 1.0;

// (Head)
    double    HCSRange                = 15.0;

// Fitness
    double    AvgSpeed = 0.00022;              // Average speed of the worm in meters per seconds
    double    BBCfit = AvgSpeed*Duration;

// Size of genotype (VC)
    

};