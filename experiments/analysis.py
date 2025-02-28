import sys
sys.path.append('..')

from run_main import make_directory

import os
import json
from collections import defaultdict


def getCellNames(network_json_data):
    return network_json_data["Nervous system"]["Cell name"]["value"]

def getNervousSystemVal(network_json_data, val):
    return network_json_data["Nervous system"][val]["value"]

def getPopNamesUnsorted(network_json_data, pop_num=6):
    cell_names = getCellNames(network_json_data)
    return cell_names[:pop_num]

def check_equal(list):
    return all(i == list[0] for i in list)

def getWeight(weights,pre_cell_ind,post_cell_ind):
    return [weight['weight'] for weight in weights if ((weight['from']==pre_cell_ind) &  (weight['to']==post_cell_ind))][0]
    
def getUniques(weights,cell_names):
    new_weights = {}
    for weight in weights:
        kk = cell_names[weight['from']-1] + cell_names[weight['to']-1]
        if kk in new_weights:
            if new_weights[kk]!=weight['weight']: 
                print("Weights not equal!")
                sys.exit(1)
        else:
            new_weights[kk]=weight['weight']
    return new_weights    
    

path_list = [] 
outFolderBases = ["varyEvolSeeds","varyEvolSeeds1"]
current = os.path.dirname(os.path.realpath(__file__))  # location of this file!
for outFolderBase in outFolderBases:
    path = current + "/" + outFolderBase
    dir_list = sorted(os.listdir(path)) 
    path_list += [path + '/' + dir for dir in dir_list]

print(len(path_list))
#sys.exit(1)

""" con_list = [
DA->VD: 15.0622
 DB->VD: 6.38036
 VD->VA: -4.39552
 VD->VB: -10.1964
 VA->DD: 15.0622
 VA->VD: 0
 VB->DD: 6.38036
 VB->VD: 0
] """


weights_list = []
for dir in path_list:
    json_file = dir + '/worm_data.json'
    with open(json_file, "r") as file:
        worm_data = json.load(file)
    cell_names = getCellNames(worm_data)    
    chem_weights = getNervousSystemVal(worm_data, 'Chemical weights')
    weights_list.append(getUniques(chem_weights,cell_names))
    


all_weights = defaultdict(list)

for d in weights_list: 
    for key, value in d.items():
        all_weights[key].append(value)    
print(all_weights)

make_directory('results')
import matplotlib.pyplot as plt   
import seaborn as sns
for key, val in all_weights.items():
#plt.hist(all_weights['VDVB'], density=True, bins=30)
    sns.displot(val, bins=10, kde=True)
    plt.savefig('results/' + key + '.png')