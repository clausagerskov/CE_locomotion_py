from run_main import run

run(
    maxGens=100,
    popSize=16,
    RandSeed=881375,
    modelName="RS18",
    modelFolder="Worm2D",
    outputFolderName="exampleRunRS18W2D",
    doEvol=True,
    overwrite=True,
)
