import sys

sys.path.append("..")

from run_main import run
from run_main import make_directory

outFolderBase = "varyEvolSeedsNet21_3"
# popSize = 10
if not make_directory(outFolderBase):
    sys.exit(1)

for ind in range(50):
    outputFolderName = outFolderBase + "/run_" + str(ind)
    run(
        outputFolderName=outputFolderName,
        # popSize=popSize,
        # mainProcessName="../main",
        # doNML=True,
        maxGens=300,
        popSize=100,
        modelName="Net21",
        modelFolder="../Worm2D",
        doEvol=True,
        overwrite=True,
    )
