//#include "VectorMatrix.h"
#include "WormRS18.h"
#include "../argUtils.h"
#include <iomanip>  // cout precision
#include "EvolutionRS18.h"

//SuppliedArgs2018 supArgs1;

int main (int argc, const char* argv[])
{
 
    std::cout << std::setprecision(10);
    long randomseed = static_cast<long>(time(NULL));

    if (argc == 2) randomseed += atoi(argv[1]);
    SuppliedArgs2018 sa1;
    sa1.randomseed = randomseed;
    if (argc>2) if (!sa1.setArgs(argc,argv,randomseed)) return 0;

    EvolutionRS18 er18(sa1);
    

    ofstream seedfile;
    seedfile.open(sa1.rename_file("seed.dat"));
    seedfile << sa1.randomseed << endl;
    seedfile.close();
    

    InitializeBodyConstants();
    
    er18.configure();

    //er18.output = 1;
   // er18.speedoutput = 1;
    RandomState rs;
    //long seed = static_cast<long>(time(NULL));
    //rs.SetRandomSeed(seed);
    rs.SetRandomSeed(sa1.randomseed);

    std::cout << std::setprecision(10);

    // Code to run simulation:
    InitializeBodyConstants();

    ifstream BestIndividualFile;
    TVector<double> bestVector(1, er18.VectSize);
    BestIndividualFile.open(sa1.rename_file("best.gen.dat"));
    BestIndividualFile >> bestVector;
    BestIndividualFile.close();
    
    er18.EvaluationFunctionOrig(bestVector, rs);

    TVector<double> phenotype(1, er18.VectSize);
    er18.GenPhenMapping(bestVector, phenotype);

    Worm18 w(phenotype, 0);
    w.InitializeState(rs);
    ofstream json_out(sa1.rename_file("worm_data.json"));
    w.writeJsonFile(json_out);
    json_out.close();
    

    
    //supArgs1.writeMessage();
    //int	VectSize = 30;
    //TVector<double> phenotype(1, VectSize);
    //Worm18 w(phenotype, 0);
    //w.writeJsonFile();

    return 0;
}