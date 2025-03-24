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
#include "Simulation.h"

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

    
    InitializeBodyConstants();
    if (atoi(getParameter(argc,argv,"--doevol","0"))) 
    {
        er->configure();
        ofstream seedfile;
        seedfile.open(er->rename_file("seed.dat"));
        seedfile << er->itsEvoPars().randomseed << endl;
        seedfile.close();
    }

    //get vector of best individual
    ifstream BestIndividualFile;
    TVector<double> bestVector(1, er->itsEvoPars().VectSize);
    BestIndividualFile.open(er->rename_file("best.gen.dat"));
    BestIndividualFile >> bestVector;
    BestIndividualFile.close();
    
    WormIzq* w = 0;
    TVector<double> phenotype(1, er->itsEvoPars().VectSize);
    er->GenPhenMapping(bestVector, phenotype);

    if (model_name == "CE") w = new WormCE(phenotype,0);
    if (model_name == "RS18") w = new Worm18(phenotype,0);
    if (model_name == "Net21") w = new Worm21(phenotype);


    // write worm_data.json if ran an evolution
    if (atoi(getParameter(argc,argv,"--doevol","0"))) {

        
        RandomState rs;
        rs.SetRandomSeed(er->itsEvoPars().randomseed);
        w->InitializeState(rs);
        ofstream json_out(er->rename_file("worm_data.json"));
    
        json j;
        w->addParsToJson(j);
        er->addParsToJson(j);
        json_out << std::setw(4) << j << std::endl;
        json_out.close();
        
    }

    bool do_nml =  atoi(getParameter(argc,argv,"--donml","0"));


    //run simulation with possibly different seed
    
    const int simrandseed =  atoi(getParameter(argc,argv,"-R","-1"));
    if (simrandseed == -1) {cout << "Seed not set properly. Exiting." << endl; return 0;}
    RandomState rs;
    rs.SetRandomSeed(simrandseed);

    if (!do_nml){
    
    er->RunSimulation(bestVector, rs);
    w->InitializeState(rs);
    Simulation s1(er->itsEvoPars());
    s1.runSimulation(*w);
    }
    else{
    if (model_name == "CE"){

        ifstream json_in(er->rename_file("worm_data.json"));
        json j;
        json_in >> j;

        {Worm2DCE w(j);
        er->RunSimulation(w, rs);
        }     

        {Worm2DCE w(j);
        w.InitializeState(rs);
        Simulation s1(er->itsEvoPars());
        s1.runSimulation(w);}


        json_in.close();
    }
    }

    delete er;
    delete w;

    return 0;
}