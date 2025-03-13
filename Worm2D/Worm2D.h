#include "TSearch.h"
#include "VectorMatrix.h"
#include "Muscles.h"
#include "WormBody.h"
#include "NervousSystem.h"
#include <nlohmann/json.hpp>



using json = nlohmann::json;

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
};

class WormIzq: public Worm2D
{
public:
    WormIzq(wormIzqParams par1_);
    //WormIzq(wormIzqParams par1_, const NervousSystemBase & n);
    virtual void InitializeState(RandomState &rs) = 0;
    virtual void Step(double StepSize, double output) = 0;
    void addParsToJson(json & j);

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

    double t; // Time

    Muscles m;
    WormBody b;
    NervousSystemBase * n_ptr;
    const wormIzqParams par1;

    // Body segment name conventions
    const int Head = 1;
    const int Tail = N_segments;

    void setUp();
    int nn(int neuronNumber, int unitNumber);
};