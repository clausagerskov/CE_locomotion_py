#include "VectorMatrix.h"
#include "WormRS18.h"

int main (int argc, const char* argv[])
{
    int	VectSize = 30;
    TVector<double> phenotype(1, VectSize);
    Worm18 w(phenotype, 0);

    return 0;
}