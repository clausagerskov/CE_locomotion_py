import sys
import os

sys.path.append("..")

from run_main import run
from run_main import make_directory

outFolderBase = "izq_runs_nets"
# popSize = 10
if not make_directory(outFolderBase):
    sys.exit(1)

path_list = []
out_path_list = []
# outFolderBases = ["varyEvolSeeds", "varyEvolSeeds1", "varyEvolSeeds2", "varyEvolSeeds3"]
inputFolderBase = "izq_selected"
current = os.path.dirname(os.path.realpath(__file__))  # location of this file!

path = current + "/" + inputFolderBase
out_path = current + "/" + outFolderBase
dir_list = sorted(os.listdir(path))
path_list += [path + "/" + dir for dir in dir_list]
out_path_list += [out_path + "/" + dir for dir in dir_list]

for input_folder, output_folder in zip(path_list, out_path_list):
    run(
        outputFolderName=output_folder,
        inputFolderName=input_folder,
        # maxGens=300,
        # popSize=100,
        modelName="Net21",
        modelFolder="../Worm2D",
        doEvol=False,
        overwrite=True,
    )
