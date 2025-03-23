#pragma once
//#include "../TSearch.h"
//#include "../VectorMatrix.h"
//#include "../Muscles.h"
//#include "../WormBody.h"
//#include "../NervousSystem.h"
//#include <nlohmann/json.hpp>
#include "jsonUtils.h"
#include "../neuromlLocal/NSBaseForW2D.h"
//using json = nlohmann::json;

#define PI 3.14159265

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



class Worm2D {
    
    public:
    //virtual void InitializeState(RandomState &rs) = 0;
    //virtual void DumpBodyState(ofstream &ofs, int skips) = 0;
    //virtual void DumpCurvature(ofstream &ofs, int skips) = 0;
    virtual void DumpActState(ofstream &ofs, int skips) = 0;
    void DumpBodyState(ofstream &ofs, int skips);

    virtual void Step(double StepSize, double output) = 0;
   
    virtual void InitializeState(RandomState &rs) = 0;

    double CoMx();
    double CoMy();
    void Curvature(TVector<double> &c);
    double Orientation();
    void AngleCurvature(TVector<double> &c);

     ~Worm2D(){if (n_ptr) delete n_ptr;}
    protected:

    Worm2D(wormIzqParams par1_, NSForW2D * n_ptr_);
    int nn(int neuronNumber, int unitNumber);
    void setUp();
    Muscles m;
    WormBody b;
    NSForW2D * const n_ptr;
    const wormIzqParams par1;
    double t; // Time

    const int Head = 1;
    const int Tail = N_segments;
};




class WormIzq : public Worm2D
{
public:
    //WormIzq(wormIzqParams par1_);
    //WormIzq(wormIzqParams par1_, const NervousSystemBase & n);
    //virtual void InitializeState(RandomState &rs) = 0;
    
    void addParsToJson(json & j);
    void writeJsonFile(ofstream & json_out);
    
    virtual void DumpParams(ofstream &ofs) = 0;

    
   // ~WormIzq(){if (n_ptr) delete n_ptr;}

    NervousSystem & n;

    protected:
    WormIzq(wormIzqParams par1);
   


    // Body segment name conventions
    
    private:

    //virtual void setUp();
    
    virtual const vector<string> getCellNames() = 0;
    virtual void addExtraParsToJson(json & j) = 0;
    //virtual Params<double> getWormParams() = 0;
    virtual vector<doubIntParamsHead> getWormParams() = 0;
};