//#include "VectorMatrix.h"
#include "WormRS18.h"
#include "WormCE.h"
#include "Worm21.h"
#include "../argUtils.h"
#include <iomanip>  // cout precision
#include "EvolutionRS18.h"
#include "EvolutionCE.h"
#include "Evolution21.h"


//SuppliedArgs2018 supArgs1;

int main (int argc, const char* argv[])
{
 
    std::cout << std::setprecision(10);

    string model_name = "";
    for (int arg = 1; arg<argc; arg+=2)
    {
        if (strcmp(argv[arg], "--modelname")==0) 
        {model_name = argv[arg+1];break;}
    }
    if (model_name == "")
    {
    cout << "Model name required for Worm2D. Exiting." << endl;
    return 0;
    }
    SuppliedArgs* sa = 0;
    Evolution* er = 0;

    long randomseed = static_cast<long>(time(NULL));
    if (model_name == "CE"){
        sa = new SuppliedArgs;
        sa->setArgs(argc,argv,randomseed);
        er = new EvolutionCE(*sa);
    }
    if (model_name == "RS18"){
        sa = new SuppliedArgs2018;
        sa->setArgs(argc,argv,randomseed);
        er = new EvolutionRS18(static_cast<SuppliedArgs2018&>(*sa));
    }
    if (model_name == "Net21"){
        cout << "making Net21" << endl;
        sa = new SuppliedArgs2021;
        sa->setArgs(argc,argv,randomseed);
        er = new Evolution21(static_cast<SuppliedArgs2021&>(*sa));
    }

    ofstream seedfile;
    seedfile.open(sa->rename_file("seed.dat"));
    seedfile << sa->randomseed << endl;
    seedfile.close();

    InitializeBodyConstants();
    if (sa->do_evol) 
    {
        cout << "config" << endl;
        er->configure();
    }
    
    
    RandomState rs;
    rs.SetRandomSeed(sa->randomseed);

    ifstream BestIndividualFile;
    TVector<double> bestVector(1, er->VectSize);
    BestIndividualFile.open(sa->rename_file("best.gen.dat"));
    BestIndividualFile >> bestVector;
    BestIndividualFile.close();
    
    er->RunSimulation(bestVector, rs);
    WormIzq* w = 0;

    TVector<double> phenotype(1, er->VectSize);
    er->GenPhenMapping(bestVector, phenotype);

    if (model_name == "CE") w = new WormCE(phenotype,0);
    if (model_name == "RS18") w = new Worm18(phenotype,0);
    if (model_name == "Net21") w = new Worm21(phenotype);

  
    w->InitializeState(rs);
    ofstream json_out(sa->rename_file("worm_data.json"));
    w->writeJsonFile(json_out);
    json_out.close();
    


    return 0;
}