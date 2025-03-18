#pragma once
#include "TSearch.h"
#include "VectorMatrix.h"
#include "../argUtils.h"
#include <functional>
#include <iomanip> 

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

struct evoPars{
   string directoryName;
   long randomseed;
   TSelectionMode SelectionMode;
   TReproductionMode ReproductionMode;
   int PopulationSize;
   int MaxGenerations;
   double MutationVariance;
   double CrossoverProbability;
   TCrossoverMode CrossoverMode;
   double MaxExpectedOffspring;
   double ElitistFraction;
   int SearchConstraint;
   int CheckpointInterval;
   bool ReEvaluationFlag;
};


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
      evolfile.close();
    //if (supArgs1_ptr) delete supArgs1_ptr; 
      if (s) delete s;
    //if (phenotype) delete phenotype;
    }

    protected:

    string rename_file(string filename){return evoPars1.directoryName + "/" + filename;}
    void configure_p1();
    void configure_p2();
    void EvolutionaryRunDisplay(int Generation, double BestPerf, double AvgPerf, double PerfVar);
    void ResultsDisplay(TSearch &s);
    
    
    Evolution(const evoPars & ep, TSearch* t)
    :s(t),evoPars1(ep)//,evoPars1(getEvoPars(sa))//bestfilename(rename_file("best.gen.dat"))
    {
     
      evolfile.open(rename_file("fitness.dat"));
      evolfile << setprecision(10);
    }
    
    //SuppliedArgs* const supArgs1_ptr;
    TSearch* const s; //(VectSize);
    //TVector<double> * phenotype; // (1, VectSize);

    const evoPars evoPars1;
    
    private:

    
    ofstream evolfile;
    //const string bestfilename; 
    

   //evoPars evoPars1;

};

