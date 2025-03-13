#include "TSearch.h"
#include "VectorMatrix.h"
#include "Worm2D.h"

class Evolution
{
    public:
    virtual void GenPhenMapping(TVector<double> &gen, TVector<double> &phen) = 0;
    virtual double EvaluationFunction(TVector<double> &v, RandomState &rs) = 0;

};

class EvolutionRS18:public Evolution
{
    public:
    void GenPhenMapping(TVector<double> &gen, TVector<double> &phen);
    double EvaluationFunction(TVector<double> &v, RandomState &rs);

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
    int VectSize = 30;

};