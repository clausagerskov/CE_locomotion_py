#include "TSearch.h"
#include "VectorMatrix.h"
#include "../argUtils.h"

class Evolution
{
    public:
    Evolution():supArgs1_ptr(nullptr){}
    virtual void GenPhenMapping(TVector<double> &gen, TVector<double> &phen) = 0;
    virtual double EvaluationFunction(TVector<double> &v, RandomState &rs) = 0;

    
    protected:
    SuppliedArgs* supArgs1_ptr;

};
