#include "Evolution.h"
#include <iomanip> 

struct evoPars{
    long randomseed;
    int SelectionMode = RANK_BASED;
    int ReproductionMode = GENETIC_ALGORITHM;
    int PopulationSize;
    int MaxGenerations;
    double MutationVariance;
    double CrossoverProbability;
    int CrossoverMode = UNIFORM;
    double MaxExpectedOffspring;
    double ElitistFraction;
    int SearchConstraint;
    int CheckpointInterval;
    bool ReEvaluationFlag;
 };


void Evolution::EvolutionaryRunDisplay(int Generation, double BestPerf, double AvgPerf, double PerfVar)
{
    evolfile << Generation << " " << BestPerf << " " << AvgPerf << " " << PerfVar << endl;
}

//void (Evolution::*RD)(TSearch &s) = &Evolution::ResultsDisplay;
//typedef void (Evolution::*RD)(TSearch &s);
//RD rd1 = &Evolution::ResultsDisplay;

//void (*ResultsDisplay)(TSearch &s) = NULL;

//ResultsDisplay = &Evolution::ResultsDisplay;

void Evolution::ResultsDisplay(TSearch &s)
{
    TVector<double> bestVector;
    ofstream BestIndividualFile;

    bestVector = s.BestIndividual();
    BestIndividualFile.open(supArgs1_ptr->rename_file("best.gen.dat"));
    //BestIndividualFile.open(bestfilename);
    BestIndividualFile << setprecision(32);
    BestIndividualFile << bestVector << endl;
    BestIndividualFile.close();
}

void Evolution::configure_p1()
{
    evolfile.open(supArgs1_ptr->rename_file("fitness.dat"));
    s->SetRandomSeed(supArgs1_ptr->randomseed);

    {typedef void (*callback_t)(int, double, double, double);
    Callback<void(int, double, double, double)>::func 
    = std::bind(&Evolution::EvolutionaryRunDisplay, this, 
        std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    callback_t func = static_cast<callback_t>(Callback<void(int, double, double, double)>::callback); 
    s->SetPopulationStatisticsDisplayFunction(func);
    }

    {typedef void (*callback_t)(TSearch&);
    Callback<void(TSearch&)>::func = std::bind(&Evolution::ResultsDisplay, this, std::placeholders::_1);
    callback_t func = static_cast<callback_t>(Callback<void(TSearch&)>::callback); 
    s->SetSearchResultsDisplayFunction(func);
    }

  
    s->SetSelectionMode(RANK_BASED);             //{FITNESS_PROPORTIONATE,RANK_BASED}
    s->SetReproductionMode(GENETIC_ALGORITHM);	// {HILL_CLIMBING, GENETIC_ALGORITHM}
    s->SetPopulationSize(supArgs1_ptr->pop_size); //96
    s->SetMaxGenerations(supArgs1_ptr->max_gens); //1000
    s->SetMutationVariance(0.1);                // For 71 parameters, an estimated avg change of 0.25 for weights (mapped to 15).
    s->SetCrossoverProbability(0.5);
    s->SetCrossoverMode(UNIFORM);              //{UNIFORM, TWO_POINT}
    s->SetMaxExpectedOffspring(1.1);
    s->SetElitistFraction(0.04);
    s->SetSearchConstraint(1);
    s->SetCheckpointInterval(0);
    s->SetReEvaluationFlag(1);

}


void Evolution::configure_p2()
{
    
    s->SetSearchTerminationFunction(NULL);

    {typedef double (*callback_t)(TVector<double> &, RandomState &);
    Callback<double(TVector<double> &, RandomState &)>::func = std::bind(&Evolution::EvaluationFunction, this, 
            std::placeholders::_1, std::placeholders::_2);
    callback_t func = static_cast<callback_t>(Callback<double(TVector<double> &, RandomState &)>::callback);
    s->SetEvaluationFunction(func);}


    //s->SetEvaluationFunction(EvaluationFunction);
    //s->SetEvaluationFunction(func);
    s->ExecuteSearch();
    evolfile.close();
}

void Evolution::configure()
{
    configure_p1();
    configure_p2();
}

   


