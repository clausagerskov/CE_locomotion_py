#include "TSearch.h"
#include "VectorMatrix.h"
#include "Muscles.h"
#include "WormBody.h"
#include "NervousSystem.h"

#define PI 3.14159265

class Worm2D {
    public:
       
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
    virtual void InitializeState(RandomState &rs) = 0;
    virtual void Step(double StepSize, double output) = 0;
    int nn(int neuronNumber, int unitNumber);
    double CoMx();
    double CoMy();
    void Curvature(TVector<double> &c);
    double Orientation();
    void AngleCurvature(TVector<double> &c);

    virtual void DumpActState(ofstream &ofs, int skips) = 0;
    virtual void DumpParams(ofstream &ofs) = 0;
    void DumpBodyState(ofstream &ofs, int skips);

    protected:

    double t; // Time

    Muscles m;
    WormBody b;
    NervousSystem n;
    const wormIzqParams par1;

    // Body segment name conventions
    const int Head = 1;
    const int Tail = N_segments;

};