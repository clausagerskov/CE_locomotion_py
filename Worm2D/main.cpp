//#include "VectorMatrix.h"
#include "WormRS18.h"
#include "WormCE.h"
#include "Worm21.h"
#include "Worm2DCE.h"
//#include "../argUtils.h"
#include <iomanip>  // cout precision
#include "EvolutionRS18.h"
#include "EvolutionCE.h"
#include "Evolution21.h"
#include "jsonUtils.h"

using json = nlohmann::json;



int main (int argc, const char* argv[])
{
 
    std::cout << std::setprecision(10);
    string model_name =  getParameter(argc,argv,"--modelname","");
    if (model_name == "")
    {
    cout << "Model name required for Worm2D. Exiting." << endl;
    return 0;
    }

    Evolution* er = 0;
    if (model_name == "CE") er = new EvolutionCE(argc,argv);
    if (model_name == "RS18") er = new EvolutionRS18(argc,argv);
    if (model_name == "Net21") er = new Evolution21(argc,argv);


    ofstream seedfile;
    seedfile.open(er->rename_file("seed.dat"));
    seedfile << er->itsEvoPars().randomseed << endl;
    seedfile.close();

    InitializeBodyConstants();
    if (atoi(getParameter(argc,argv,"--doevol","0"))) er->configure();
    
    
    
    RandomState rs;
    rs.SetRandomSeed(er->itsEvoPars().randomseed);

    ifstream BestIndividualFile;
    TVector<double> bestVector(1, er->itsEvoPars().VectSize);
    BestIndividualFile.open(er->rename_file("best.gen.dat"));
    BestIndividualFile >> bestVector;
    BestIndividualFile.close();
    
    er->RunSimulation(bestVector, rs);
    WormIzq* w = 0;

    TVector<double> phenotype(1, er->itsEvoPars().VectSize);
    er->GenPhenMapping(bestVector, phenotype);

    if (model_name == "CE") w = new WormCE(phenotype,0);
    if (model_name == "RS18") w = new Worm18(phenotype,0);
    if (model_name == "Net21") w = new Worm21(phenotype);

  
    w->InitializeState(rs);
    ofstream json_out(er->rename_file("worm_data.json"));

    json j;
    w->addParsToJson(j);
    er->addParsToJson(j);
    json_out << std::setw(4) << j << std::endl;

    //w->writeJsonFile(json_out);
    json_out.close();
    if (model_name == "CE"){

        Worm2DCE w(j);
        
    }
    delete er;
    delete w;

    return 0;
}