#include "VectorMatrix.h"
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
    
    EvolutionRS18 er18;
    SuppliedArgs2018 & supArgs1 = er18.supArgs1;
    supArgs1.randomseed = randomseed;
    if (argc>2) if (!supArgs1.setArgs(argc,argv,randomseed)) return 0;

    ofstream seedfile;
    seedfile.open(supArgs1.rename_file("seed.dat"));
    seedfile << supArgs1.randomseed << endl;
    seedfile.close();
    

    InitializeBodyConstants();
    
    er18.configure();
    return 0;
    //supArgs1.writeMessage();
    //int	VectSize = 30;
    //TVector<double> phenotype(1, VectSize);
    //Worm18 w(phenotype, 0);
    //w.writeJsonFile();

    return 0;
}