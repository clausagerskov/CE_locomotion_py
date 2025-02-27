
import os
import json

def getCellNames(network_json_data):
    return network_json_data["Nervous system"]["Cell name"]["value"]


def getPopNamesUnsorted(network_json_data, pop_num=6):
    cell_names = getCellNames(network_json_data)
    return cell_names[:pop_num]


outFolderBase = "varyEvolSeeds"
current = os.path.dirname(os.path.realpath(__file__))  # location of this file!
path = current + "/" + outFolderBase
dir_list = os.listdir(path) 

pre_cell = 'DA'
post_cell = 'VD'

for dir in dir_list:
    json_file = path + '/' + dir + '/worm_data.json'
    with open(json_file, "r") as file:
        worm_data = json.load(file)
    pop_names = getPopNamesUnsorted(worm_data)
    pre_cell_ind = []