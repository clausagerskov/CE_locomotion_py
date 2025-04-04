//#include "VectorMatrix.h"
#include "WormRS18.h"
#include "WormCE.h"
#include "Worm21.h"
//#include "Worm2DCE.h"
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

    bool do_evol = atoi(getParameter(argc,argv,"--doevol","0"));
    if (do_evol) 
    {
        er->configure();
        ofstream seedfile;
        seedfile.open(er->rename_file("seed.dat"));
        seedfile << er->itsEvoPars().randomseed << endl;
        seedfile.close();
    }

    //get vector of best individual
   
    TVector<double> bestVector(1, er->itsEvoPars().VectSize);

    {ifstream BestIndividualFile;
    BestIndividualFile.open(er->rename_file("best.gen.dat"));
    BestIndividualFile >> bestVector;
    BestIndividualFile.close();}
    
    TVector<double> phenotype(1, er->itsEvoPars().VectSize);
    er->GenPhenMapping(bestVector, phenotype);

    bool do_json = 1;

    // write worm_data.json if ran an evolution
    if (do_json) {

        Worm2D* w = 0;
        

        if (model_name == "CE") w = new WormCE(phenotype,0);
        if (model_name == "RS18") w = new Worm18(phenotype,0);
        if (model_name == "Net21") w = new Worm21(phenotype);


        RandomState rs;
        rs.SetRandomSeed(er->itsEvoPars().randomseed);
        w->InitializeState(rs); 

        ofstream nsdump(er->rename_file("NSdump.dat"));
        nsdump << dynamic_cast<NervousSystem&>(w->itsNS());
        nsdump.close();


        ofstream json_out(er->rename_file("worm_data.json"));    
        json j;
        w->addParsToJson(j);
        er->addParsToJson(j);
        json_out << std::setw(4) << j << std::endl;
        json_out.close();
        delete w;
        
    }

    bool do_nml =  atoi(getParameter(argc,argv,"--donml","0"));

    //run simulation with possibly different seed
    
    const int simrandseed =  atoi(getParameter(argc,argv,"-R","-1"));
    if (simrandseed == -1) {cout << "Seed not set properly. Exiting." << endl; return 0;}
    

    if (!do_nml){
    
    //er->RunSimulation(bestVector, rs);

    Worm2D* w = 0;

   
    cout << "making worm" << endl;

    if (model_name == "CE") w = new WormCE(phenotype,0);
    if (model_name == "RS18") w = new Worm18(phenotype,0);
    if (model_name == "Net21") w = new Worm21(phenotype);

    cout << "making simulation" << endl;
    {RandomState rs;
    rs.SetRandomSeed(simrandseed);
    er->RunSimulation(*w, rs);}

    {ofstream phenfile(er->rename_file("phenotype.dat"));
    w->DumpParams(phenfile);
    phenfile.close();}

    {RandomState rs;
    rs.SetRandomSeed(simrandseed);
    w->InitializeState(rs);}
    w->initForSimulation();
    simPars sp1 = {er->itsEvoPars().directoryName,
        er->itsEvoPars().skip_steps, 60, 50, er->itsEvoPars().StepSize};
    Simulation s1(sp1);
    s1.runSimulation(*w);

    delete w;


    }
    else{

    ifstream json_in(er->rename_file("worm_data.json"));
    json j;
    json_in >> j;

    Worm2D* w = 0;
    if (model_name == "CE") w = new Worm2DCE(j);
    if (model_name == "Net21") w = new Worm2D21(j);

    {RandomState rs;
    rs.SetRandomSeed(simrandseed);
    er->RunSimulation(*w, rs);}

    {RandomState rs;
    rs.SetRandomSeed(simrandseed);
    w->InitializeState(rs);
    w->initForSimulation();
    simPars sp1 = {er->itsEvoPars().directoryName,
        er->itsEvoPars().skip_steps, 60, 50, er->itsEvoPars().StepSize};
    Simulation s1(sp1);
    s1.runSimulation(*w);}

    delete w;
    json_in.close();
    }

    delete er;
    

    return 0;
}