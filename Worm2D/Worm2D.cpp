#include "Worm2D.h"


WormIzq::WormIzq(wormIzqParams par1_):par1(par1_),n_ptr(new NervousSystem)
{
    setUp();
}

/* WormIzq::WormIzq(wormIzqParams par1_, const NervousSystemBase & n):par1(par1_),n_ptr(n.clone())
{
    setUp();
} */

void WormIzq::setUp()
{
    m.SetMuscleParams(par1.N_muscles, par1.T_muscle);

}

int WormIzq::nn(int neuronNumber, int unitNumber)
{
    return neuronNumber+((unitNumber-1)*par1.N_neuronsperunit);
}

void WormIzq::InitializeState(RandomState &rs)
{
    t = 0.0;
    b.InitializeBodyState();
    m.InitializeMuscleState();
}


double WormIzq::CoMx()
{
    double temp = 0.0;
    for (int i = 1; i <= N_rods; i++) {
        temp += b.X(i);
    }
    return temp/N_rods;
}

double WormIzq::CoMy()
{
    double temp = 0.0;
    for (int i = 1; i <= N_rods; i++) {
        temp += b.Y(i);
    }
    return temp/N_rods;
}

void WormIzq::Curvature(TVector<double> &c)
{
    double dx1,dy1,dx2,dy2,a,a1,a2,seg;
    int k=1;

    for (int i = 3; i < N_segments-1; i+=2)
    {
        dx1 = b.X(i) - b.X(i-2);
        dy1 = b.Y(i) - b.Y(i-2);
        dx2 = b.X(i+2) - b.X(i);
        dy2 = b.Y(i+2) - b.Y(i);

        a1 = atan2(dy1,dx1);
        a2 = atan2(dy2,dx2);

        if (a1 > PI/2 and a2 < -PI/2)
            a = (a1 - 2*PI) - a2;
        else
            if (a1 < -PI/2 and a2 > PI/2)
                a = a1 - (a2 - 2*PI);
            else
                a = a1-a2;

        seg = sqrt(pow(b.X(i-2)-b.X(i+2),2) + pow(b.Y(i-2)-b.Y(i+2),2));
        c(k) = (2*sin(a)/seg)/1000;
        k++;
    }
}

double WormIzq::Orientation()
{
    return atan2(b.Y(Head)-b.Y(Tail),b.X(Head)-b.X(Tail));
}

void WormIzq::AngleCurvature(TVector<double> &c)
{
  double dx1,dy1,dx2,dy2,a,a1,a2,seg;
  int k=1;

  for (int i = 3; i < N_segments-1; i+=2)
  {
    dx1 = b.X(i) - b.X(i-2);
    dy1 = b.Y(i) - b.Y(i-2);
    dx2 = b.X(i+2) - b.X(i);
    dy2 = b.Y(i+2) - b.Y(i);

    a1 = atan2(dy1,dx1);
    a2 = atan2(dy2,dx2);

    if (a1 > PI/2 and a2 < -PI/2)
    a = (a1 - 2*PI) - a2;
    else
    if (a1 < -PI/2 and a2 > PI/2)
    a = a1 - (a2 - 2*PI);
    else
    a = a1-a2;
    c(k) = a;
    k++;
  }
}

void WormIzq::DumpBodyState(ofstream &ofs, int skips)
{
    static int tt = skips;

    if (++tt >= skips) {
        tt = 0;

        ofs << t;
        // Body
        for (int i = 1; i <= N_rods; i++)
        {
            ofs <<  " " << b.X(i) << " " << b.Y(i) << " " << b.Phi(i);
        }
        ofs << "\n";
    }
}

void WormIzq::addParsToJson(json & j)
{
    appendBodyToJson(j, b);
    appendMuscleToJson(j,m);
    string nsHead = "Nervous system";
    getNSJson(static_cast<NervousSystem&>(*n_ptr), j[nsHead]);
    appendCellNamesToJson(j[nsHead], getCellNames(), par1.N_units);
    Params<double> par = getWormParams();
    appendToJson<double>(j["Worm"],par);
    
    addExtraParsToJson(j);
}
