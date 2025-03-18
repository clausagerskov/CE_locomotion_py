#include "Evolution.h"
#include <iomanip> 



evoPars Evolution::getEvoPars(const SuppliedArgs & sa)
{
    return {sa.output_dir_name, sa.randomseed, RANK_BASED, GENETIC_ALGORITHM, 
        sa.pop_size, sa.max_gens, 0.1, 0.5, UNIFORM, 
        1.1, 0.04, 1, 0, 0};

}

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
    BestIndividualFile.open(rename_file("best.gen.dat"));
    //BestIndividualFile.open(bestfilename);
    BestIndividualFile << setprecision(32);
    BestIndividualFile << bestVector << endl;
    BestIndividualFile.close();
}

void Evolution::configure_p1()
{
   
    s->SetRandomSeed(evoPars1.randomseed);

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


    s->SetSelectionMode(evoPars1.SelectionMode);             //{FITNESS_PROPORTIONATE,RANK_BASED}
    s->SetReproductionMode(evoPars1.ReproductionMode);	// {HILL_CLIMBING, GENETIC_ALGORITHM}
    s->SetPopulationSize(evoPars1.PopulationSize); //96
    s->SetMaxGenerations(evoPars1.MaxGenerations); //1000
    s->SetMutationVariance(evoPars1.MutationVariance);                // For 71 parameters, an estimated avg change of 0.25 for weights (mapped to 15).
    s->SetCrossoverProbability(evoPars1.CrossoverProbability);
    s->SetCrossoverMode(evoPars1.CrossoverMode);              //{UNIFORM, TWO_POINT}
    s->SetMaxExpectedOffspring(evoPars1.MaxExpectedOffspring);
    s->SetElitistFraction(evoPars1.ElitistFraction);
    s->SetSearchConstraint(evoPars1.SearchConstraint);
    s->SetCheckpointInterval(evoPars1.CheckpointInterval);
    s->SetReEvaluationFlag(evoPars1.ReEvaluationFlag);

}


void Evolution::configure_p2()
{
    
    s->SetSearchTerminationFunction(NULL);

    {typedef double (*callback_t)(TVector<double> &, RandomState &);
    Callback<double(TVector<double> &, RandomState &)>::func = std::bind(&Evolution::EvaluationFunction, this, 
            std::placeholders::_1, std::placeholders::_2);
    callback_t func = static_cast<callback_t>(Callback<double(TVector<double> &, RandomState &)>::callback);
    s->SetEvaluationFunction(func);}

    s->ExecuteSearch();
  
}

void Evolution::configure()
{
    configure_p1();
    configure_p2();
}

   


