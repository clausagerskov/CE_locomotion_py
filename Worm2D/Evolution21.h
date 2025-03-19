#include "TSearch.h"
#include "VectorMatrix.h"
#include "Evolution.h"

class Evolution21:public Evolution
{
    public:
    Evolution21(const SuppliedArgs2021 & sa):Evolution(getEvoPars(sa),44,getSimPars(sa)){}

    void GenPhenMapping(TVector<double> &gen, TVector<double> &phen);
    double EvaluationFunction(TVector<double> &v, RandomState &rs);
    double EvaluationFunction1(TVector<double> &v, RandomState &rs);
    double EvaluationFunction2(TVector<double> &v, RandomState &rs);
    double EvaluationFunction2Output(TVector<double> &v, RandomState &rs);
    void RunSimulation(TVector<double> &v, RandomState &rs);
    
    int finish_Bosc(int Generation,double BestPerf,double AvgPerf,double PerfVar);
   

    //SuppliedArgs2018& supArgs1 = static_cast<SuppliedArgs2018&>(*supArgs1_ptr);
    //SuppliedArgs2018* const makeArgsPtr() {cout << "set pointer" << endl; return new SuppliedArgs2018();}
    //TSearch* const makeTSearchPtr() {return new TSearch(VectSize);}

    //SuppliedArgs2021& supArgs1;
     // Size of genotype

   // const int VectSize = 44;
    
    
    protected:
    simPars getSimPars(const SuppliedArgs &);
    evoPars getEvoPars(const SuppliedArgs2021 & sa);
    int skip_steps = 10;

     // Integration parameters
     const double Duration = 40.0;       //
     const double Transient = 10.0;       //
     const double StepSize = 0.005;
     const int N_curvs = 23;             // Number of cuvature points
     
     double OSCT = 0.25 * Duration; // Cap for oscillation evaluation
     const double agarfreq = 0.44;
     
     // Genotype -> Phenotype Mapping (Ventral cord)
     const double	BiasRange				= 15.0;
     const double    SCRange                 = 15.0;
     const double    CSRange                 = 15.0;
     const double    TauMin                 = 0.1;
     const double    TauMax                 = 2.5;
     const double    ESRange                 = 2.0;
     const double    NMJmax                  = 1.2;
     const double    IIRange                 = 15.0;
     
     // Fitness
     const double    AvgSpeed = 0.00022;             // Average speed of the worm in meters per seconds
     const double    BBCfit = AvgSpeed*Duration;
     void configure_p2();
   
};