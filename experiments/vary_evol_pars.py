import sys
sys.path.append('..')

from run_main import run
from run_main import make_directory

outFolderBase = 'varyEvolSeeds'
popSize = 10
make_directory(outFolderBase)

for ind in range(5):
    outputFolderName = outFolderBase + '/run_' + str(ind)
    run(outputFolderName=outputFolderName, popSize=popSize, mainProcessName = '../main')
