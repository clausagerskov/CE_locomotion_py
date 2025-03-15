#include "TSearch.h"
#include "VectorMatrix.h"
#include "Muscles.h"
#include "WormBody.h"
#include "NervousSystem.h"
//#include <nlohmann/json.hpp>
#include "jsonUtils.h"


//using json = nlohmann::json;

#define PI 3.14159265

class Worm2D {
    
    public:
    virtual void addParsToJson(json & j) = 0;

    };

struct wormIzqParams
{
    
    int N_neuronsperunit;
    int N_muscles;
    double T_muscle;
    int N_units;

    const doubIntParamsHead getParams() const
    {
        doubIntParamsHead var1;
        var1.parDoub.head = "Worm";
        var1.parInt.head = "Worm";
        var1.parDoub.names = {"T_muscle"};
        var1.parDoub.vals = {T_muscle};
        var1.parInt.names = {"N_neuronsperunit", "N_muscles", "N_units"};
        var1.parInt.vals = {N_neuronsperunit, N_muscles, N_units};
        return var1;
    }

};

class WormIzq: public Worm2D
{
public:
    WormIzq(wormIzqParams par1_);
    //WormIzq(wormIzqParams par1_, const NervousSystemBase & n);
    virtual void InitializeState(RandomState &rs) = 0;
    virtual void Step(double StepSize, double output) = 0;
    void addParsToJson(json & j);
    void writeJsonFile(ofstream & json_out);

    double CoMx();
    double CoMy();
    void Curvature(TVector<double> &c);
    double Orientation();
    void AngleCurvature(TVector<double> &c);

    virtual void DumpActState(ofstream &ofs, int skips) = 0;
    virtual void DumpParams(ofstream &ofs) = 0;
    void DumpBodyState(ofstream &ofs, int skips);
    
    ~WormIzq(){if (n_ptr) delete n_ptr;}

    protected:
    
    int nn(int neuronNumber, int unitNumber);
    double t; // Time

    Muscles m;
    WormBody b;
    NervousSystemBase * n_ptr;
    const wormIzqParams par1;

    // Body segment name conventions
    const int Head = 1;
    const int Tail = N_segments;
    private:

    virtual void setUp();
    
    virtual const vector<string> getCellNames() = 0;
    virtual void addExtraParsToJson(json & j) = 0;
    //virtual Params<double> getWormParams() = 0;
    virtual vector<doubIntParamsHead> getWormParams() = 0;
};