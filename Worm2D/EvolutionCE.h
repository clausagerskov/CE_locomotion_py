#include "TSearch.h"
#include "VectorMatrix.h"
#include "Evolution.h"

class EvolutionCE:public Evolution
{
    public:
    EvolutionCE(const SuppliedArgs & sa):Evolution(getEvoPars(sa),17){}

    void GenPhenMapping(TVector<double> &gen, TVector<double> &phen);
    double EvaluationFunction(TVector<double> &v, RandomState &rs);
    double Evaluation(TVector<double> &v, RandomState &rs, int direction);
    evoPars getEvoPars(const SuppliedArgs & sa);
    double save_traces(TVector<double> &v, RandomState &rs);
    void RunSimulation(TVector<double> &v, RandomState &rs);

    
    //SuppliedArgs2018& supArgs1 = static_cast<SuppliedArgs2018&>(*supArgs1_ptr);
    //SuppliedArgs2018* const makeArgsPtr() {cout << "set pointer" << endl; return new SuppliedArgs2018();}
    //TSearch* const makeTSearchPtr() {return new TSearch(VectSize);}

    //SuppliedArgs& supArgs1;
     // Size of genotype
    //const int VectSize = 17;

    
    private:
    int skip_steps = 10;
    
    
    // Integration parameters
    const int Duration = 24;
    const double Transient = 8.0;
    const double StepSize = 0.005;
    const int N_curvs = 23;
    
    // Fitness traj
    const double    AvgSpeed = 0.0001; //0.00022;              // Average speed of the worm in meters per seconds
    const double    BBCfit = AvgSpeed*Duration;
    
    // Genotype -> Phenotype Mapping Ranges
    const double    BiasRange               = 16.0; //15.0;
    const double    SCRange                 = 16.0; //15.0;
    const double    CSRange                 = 16.0; //15.0;
    const double    ESRange                 = 2.0;
    const double    SRmax                   = 200.0;
    const double    NMJmax                  = 0.8; //1.2;
    const double    NMJmin                  = 0.0;
    
    const int SR_A = 1;
    const int SR_B = 2;
    
   
};