#include "../c302NervousSystem.h"
#include <fstream>

using namespace std;


int main (int argc, const char* argv[])
{

c302NervousSystem n("main_sim", "parent");

n.SetPopStructure("DA DB DD VD VA VB", 10);
n.skipCalc = 1;

ofstream fout("testc302NervousSystem-output.dat");
for (int i=0;i<10000;i++){
if (i==3000) n.SetChemicalSynapseWeight(3, 1, 1);

if (i==6000) n.SetNeuronExternalInput(53, 1);

if (i==9000) n.SetNeuronBias(35, 1);


n.EulerStep(1);
fout << i;
n.writeOutputValues(fout);
fout << endl;
//cout << n.NeuronBias(55) << endl;
//cout << n.NeuronOutput(1) << endl;
//cout << n.NeuronOutput(2) << endl;
//cout << n.NeuronOutput(3) << endl;
}
fout.close();

return 0;
}