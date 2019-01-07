# Chaste VertexComparison2013 project

A collection of test script files and their results, which use the central chaste cancer code to recreate the results presented in our manuscript ["Implementing vertex dynamics models of cell populations in biology within a consistent computational framework"](https://doi.org/10.1016/j.pbiomolbio.2013.09.003)

There are two directories - `src` and `test`. 

The `src` folder contains a number of `.hpp` and `.cpp` files defining classes that are used in the tests:

```
ModifiedWelikyOsterForce
MotileCellForce
MphaseGrowthTargerAreaModifer
NagaiHondaDifferentialAdhesionForce
ObstructionBoundaryCondition
SimpleVolumeBasedStochasticCellCycleModel
SimpleWntUniformDistCellCycleModel
VertexAngleForce
VertexDiffusionForce
```

The `test` folder contains `.hpp` files defining simulations that are described in the manuscript:

```
TestVertexActiveMigration
TestVertexCellSorting
TestVertexMonolayers
TestVertexRestrictedGeometry
TestVertexWoundHeal
```

Please note that this code has been updated to work with Chaste version 2018.1 and hence differs from the code described in the published manuscript (which worked with Chaste version 3.1).
