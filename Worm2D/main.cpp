#include "VectorMatrix.h"
#include "WormRS18.h"
#include "../argUtils.h"
#include <iomanip>  // cout precision

//SuppliedArgs2018 supArgs1;

int main (int argc, const char* argv[])
{
 
    std::cout << std::setprecision(10);
    long randomseed = static_cast<long>(time(NULL));

    if (argc == 2) randomseed += atoi(argv[1]);
    
    //supArgs1.max_gens = 40;
    //supArgs1.pop_size = 56;

    SuppliedArgs2018 supArgs1;
    supArgs1.randomseed = randomseed;
    if (argc>2) if (!supArgs1.setArgs(argc,argv,randomseed)) return 0;



    //supArgs1.writeMessage();
    int	VectSize = 30;
    TVector<double> phenotype(1, VectSize);
    Worm18 w(phenotype, 0);
    w.writeJsonFile();

    return 0;
}