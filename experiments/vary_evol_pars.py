import sys

sys.path.append("..")

from run_main import run
from run_main import make_directory

outFolderBase = "varyEvolSeedsNML1"
popSize = 10
if not make_directory(outFolderBase):
    sys.exit(1)

for ind in range(500):
    outputFolderName = outFolderBase + "/run_" + str(ind)
    run(
        outputFolderName=outputFolderName,
        popSize=popSize,
        mainProcessName="../main",
        doNML=True,
    )
