#include "TSearch.h"
#include "VectorMatrix.h"
#include "../argUtils.h"
#include <functional>

template <typename T>
struct Callback;

template <typename Ret, typename... Params>
struct Callback<Ret(Params...)> {
   template <typename... Args> 
   static Ret callback(Args... args) {                    
      return func(args...);  
   }
   static std::function<Ret(Params...)> func; 
};

template <typename Ret, typename... Params>
std::function<Ret(Params...)> Callback<Ret(Params...)>::func;

struct evoPars;

class Evolution
{
    public:
    //Evolution():supArgs1_ptr(makeArgsPtr()),s(makeTSearchPtr()){}
    virtual void GenPhenMapping(TVector<double> &gen, TVector<double> &phen) = 0;
    virtual double EvaluationFunction(TVector<double> &v, RandomState &rs) = 0;
    
    virtual void configure();
    //virtual SuppliedArgs* const makeArgsPtr(){return NULL;}
    //virtual TSearch* const makeTSearchPtr(){return NULL;}
    //void EvolutionaryRunDisplay(int Generation, double BestPerf, double AvgPerf, double PerfVar);

    virtual ~Evolution()
    {
      //evolfile.close();
    if (supArgs1_ptr) delete supArgs1_ptr; 
    if (s) delete s;
    //if (phenotype) delete phenotype;
    }

    protected:
    void configure_p1();
    void configure_p2();
    void EvolutionaryRunDisplay(int Generation, double BestPerf, double AvgPerf, double PerfVar);
    void ResultsDisplay(TSearch &s);
    Evolution(SuppliedArgs* sa, TSearch* t):supArgs1_ptr(sa),s(t){}
    
    SuppliedArgs* const supArgs1_ptr;
    TSearch* const s; //(VectSize);
    //TVector<double> * phenotype; // (1, VectSize);

   evoPars evoPars1;
    private:

    
    ofstream evolfile;
    //string bestfilename; 
   //evoPars evoPars1;

};

