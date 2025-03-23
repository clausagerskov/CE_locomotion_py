
#include "jsonUtils.h"
#include <iomanip>
//#include "../argUtils.h"

using json = nlohmann::json;

struct toFromWeight{
    
    toFromWeight(weightentry w_val, int to_val){w=w_val;to=to_val;}
    toFromWeight(){}
    weightentry w;
    int to;
};


void to_json(json & j, const weightentry & w)
{
  j = json{{"from", w.from}, {"weight", w.weight}};
}

void to_json(json & j, const toFromWeight & w)
{
  j = json{{"to", w.to}, {"from", w.w.from}, {"weight", w.w.weight}};
}

void from_json(const json& j, toFromWeight & w) 
{
        j.at("to").get_to(w.to);
        j.at("from").get_to(w.w.from);
        j.at("weight").get_to(w.w.weight);
}





Params<double> getBodyParams(WormBody& b)
{

Params<double> par;

par.names = {"Medium", "L_worm", "R_min", "C_agar_par_total", 
"C_agar_perp_total", "C_water_par_total", "C_water_perp_total", "kappa_L", "kappa_D", 
"kappa_M0", "beta_L", "beta_D", "beta_M0", "delta_M"};

par.vals = {Medium, L_worm, R_min, C_agar_par_total, C_agar_perp_total, 
C_water_par_total, C_water_perp_total, kappa_L, kappa_D, kappa_M0, beta_L, 
beta_D, beta_M0, delta_M};

par.messages_inds.resize(par.vals.size());

for (size_t i=0;i<par.messages_inds.size();i++) par.messages_inds[i]=i;

par.messages = {    
"Normalized medium drag coefficient (0 = water, 1 = agar)",
"Length of worm in m",
"Minor radius of prolate ellipse body in m",
"Total tangential drag coefficient for agar in kg/s",
"Total rod normal drag coefficient in agar in kg/s",
"Total rod tangential drag coefficient for water in kg/s",
"Total rod normal drag coefficient for water in kg/s",
"Lateral spring constant in kg/s",
"Diagonal spring constant in kg/s",
"Baseline active muscle spring constant in kg/s",
"Lateral passive damping constant in s",
"Diagonal passive damping constant in s",
"Baseline active damping constant in s",
"Rest muscle length scaling constant"

};

return par;

}

Params<int> getBodyParamsInts(WormBody& b)
{
Params<int> par;
par.names = {"N_segments"};
par.vals = {N_segments};
par.messages = {"Number of body segments on each side, dorsal and ventral"};
par.messages_inds  = {0};
return par;
}

void appendBodyToJson(json & j, WormBody& b)
{
{Params<double> par = getBodyParams(b);
    appendToJson<double>(j["Body"],par);}
{Params<int> par = getBodyParamsInts(b);
    appendToJson<int>(j["Body"],par);}
} 



Params<double> getMusclesParamsDouble(Muscles & m)
{
Params<double> par;
par.names = {"T_muscle"};
par.vals = {m.T_muscle};
return par;
}

Params<int> getMusclesParamsInt(Muscles & m)
{
Params<int> par;
par.names = {"Nmuscles"};
par.vals = {m.Nmuscles};
return par;
}


void appendMuscleToJson(json & j, Muscles & m)
{
{Params<double> par = getMusclesParamsDouble(m);
appendToJson<double>(j["Muscle"],par);}
{Params<int> par = getMusclesParamsInt(m);
appendToJson<int>(j["Muscle"],par);}
}



Params< vector<double> > getNervousSysParamsDoubleNH(NervousSystem& c)
{
Params< vector<double> > par;
par.names = {"taus", "biases", "gains", "outputs", "states", "paststates", "Rtaus", "externalinputs"};
par.vals = {
getVector<double>(c.taus, c.size), 
getVector<double>(c.biases, c.size), 
getVector<double>(c.gains, c.size),
getVector<double>(c.outputs, c.size),
getVector<double>(c.states, c.size),
getVector<double>(c.paststates, c.size),
getVector<double>(c.Rtaus, c.size),
getVector<double>(c.externalinputs, c.size),
};
return par;
}

Params<int> getNervousSysParamsIntNH(NervousSystem& c)
{
Params<int> par;    
par.names = {"size", "maxchemcons", "maxelecconns"};
par.vals = {c.size, c.maxchemconns, c.maxelecconns};
return par;
}

Params< vector<int> > getNervousSysVecInt(NervousSystem& c)
{
Params< vector<int> > par;
par.names = {"NumChemicalConns", "NumElectricalConns"};
par.vals = {
getVector<int>(c.NumChemicalConns, c.size), 
getVector<int>(c.NumElectricalConns, c.size), 
};
return par;
}

void appendMatrixToJson(json & j, TMatrix<weightentry> & vec, TVector<int> & sizes, int tot_size)
{    
    vector<toFromWeight> newvec;
    for (int i=1; i<=tot_size; i++){    
        for (int j=1; j<=sizes[i]; j++) { 
            toFromWeight tv(vec[i][j], i);
            newvec.push_back(tv);}        
    }
    j["value"] = newvec;

}


void appendNSToJson(json & j, NervousSystem& c)
{
    j["Chemical weights"]["message"] = "chemical weights in sparse format";
    appendMatrixToJson(j["Chemical weights"], c.chemicalweights, c.NumChemicalConns, c.size);
    appendMatrixToJson(j["Electrical weights"], c.electricalweights, c.NumElectricalConns, c.size);
    j["Electrical weights"]["message"] = "electrical weights in sparse format";
}


void appendAllNSJson( json & j, NervousSystem & n)
{

{Params<vector<double> > parvec = getNervousSysParamsDoubleNH(n);
appendToJson<vector<double> >(j,parvec);}
        
{Params<int> parvec = getNervousSysParamsIntNH(n);
appendToJson<int>(j,parvec);}

{Params< vector<int> > parvec = getNervousSysVecInt(n);
appendToJson<vector<int> >(j,parvec);}

appendNSToJson(j, n);

}

Params< vector<string> > getNervousSysCellNames(const vector<string> & cell_names, int n_units)
{
Params< vector<string> > par;
par.names = {"Cell name"};
vector<string> cell_names_all;
for (int i=0;i<n_units;i++) cell_names_all.insert(cell_names_all.end(),cell_names.begin(),cell_names.end());
par.vals = {cell_names_all};
return par;
}

void appendCellNamesToJson(json & j, const vector<string> & cell_names, const int & num_reps)
{
    Params< vector<string> > parvec = getNervousSysCellNames(cell_names, num_reps);
    appendToJson<vector<string> >(j,parvec);
}
