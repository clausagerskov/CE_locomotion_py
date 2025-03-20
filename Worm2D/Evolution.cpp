#include "Evolution.h"
//#include <string>
#include <string.h>
#include <sys/stat.h>
//#include <stdio.h>
#include <iostream>

const char* getParameter(int argc, const char* argv[], string parName, const char* defaultval)
{    
   const char* retval = defaultval;
   if (((argc-1) % 2) != 0)
   {cout << "The arguments are not configured correctly." << endl;exit(1);}
   for (int arg = 1; arg<argc; arg+=2) 
   if (strcmp(argv[arg],parName.c_str())==0) {retval = argv[arg+1];break;}
   return retval;
}

void Evolution::addParsToJson(json & j)
{  
    
    doubIntParamsHead par1pars = evoPars1.getParams();
    appendToJson<double>(j[par1pars.parDoub.head],par1pars.parDoub);
    appendToJson<long>(j[par1pars.parInt.head],par1pars.parInt);

    addExtraParsToJson(j);
}

evoPars Evolution::setPars(int argc, const char* argv[], evoPars ep1){


    if (((argc-1) % 2) != 0)
     {cout << "The arguments are not configured correctly." << endl;exit(1);}
    
    bool seed_flag = 1;
    
    for (int arg = 1; arg<argc; arg+=2)
    { 
    if (strcmp(argv[arg],"--maxgens")==0) ep1.MaxGenerations = atoi(argv[arg+1]);
    //if (strcmp(argv[arg],"--doevol")==0) do_evol = atoi(argv[arg+1]);
    //if (strcmp(argv[arg],"--dorandinit")==0) simRandomInit = atoi(argv[arg+1]);
    //if (strcmp(argv[arg],"--skipOrigSim")==0) skipOrigSim = atoi(argv[arg+1]);
    //if (strcmp(argv[arg],"--donml")==0) do_nml = atoi(argv[arg+1]);

    if (strcmp(argv[arg],"--folder")==0) {
      ep1.directoryName= argv[arg+1];
      struct stat sb;
      if (stat(ep1.directoryName.c_str(), &sb) != 0) 
      {cout << "Directory doesn't exist." << endl;exit(1);}
    }

    if (seed_flag){ 
    if (strcmp(argv[arg],"-R")==0) ep1.randomseed = atoi(argv[arg+1]);
    if (strcmp(argv[arg],"-r")==0) 
    {long randomseed1 = static_cast<long>(time(NULL));
           ep1.randomseed = randomseed1 + atoi(argv[arg+1]);
    }
    seed_flag = 0;
    }
    //if (strcmp(argv[arg], "--modelname")==0) model_name = argv[arg+1];
    if (strcmp(argv[arg],"-p")==0) ep1.PopulationSize = atoi(argv[arg+1]);
    if (strcmp(argv[arg],"-d")==0) ep1.Duration = atoi(argv[arg+1]);
  
    //if (strcmp(argv[arg],"--nervous")==0) nervousSystemNameForSim = argv[arg+1];
    
    }
    return ep1;

}

void Evolution::EvolutionaryRunDisplay(int Generation, double BestPerf, double AvgPerf, double PerfVar)
{
    evolfile << Generation << " " << BestPerf << " " << AvgPerf << " " << PerfVar << endl;
}



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

   


