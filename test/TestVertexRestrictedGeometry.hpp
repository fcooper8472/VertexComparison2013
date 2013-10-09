/*

Copyright (c) 2005-2013, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTVERTEXRESTRICTEDGEOMETRY_HPP_
#define TESTVERTEXRESTRICTEDGEOMETRY_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "SimpleWntUniformDistCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "NagaiHondaMPhaseGrowthForce.hpp"
#include "ModifiedWelikyOsterForce.hpp"
#include "VertexAngleForce.hpp"
#include "PlaneBasedCellKiller.hpp"

#include "StochasticDurationCellCycleModel.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "ObstructionBoundaryCondition.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "Debug.hpp"


/**********************************************
 * THIS TEST WORKS WITH RELEASE 3.0 OF CHASTE *
 **********************************************/

class TestVertexRestrictedGeometry: public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    void TestVertexWithObstuction() throw (Exception)
    {
        // Create a regular compressed 2D vertex-based mesh, of size 12 by 1 elements
        CylindricalHoneycombVertexMeshGenerator generator(12, 1, true, 0.6);

        // Impose periodicity at the left- and right-hand boundaries of the domain
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        // Set up a vector of cells with given cell-cycle model
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<SimpleWntUniformDistCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,p_mesh->GetNumElements(), p_transit_type);

        // Create a vertex-based cell population and set up a cell-based simulation
        VertexBasedCellPopulation<2> population(*p_mesh, cells);
        OffLatticeSimulation<2> simulator(population);
        simulator.SetOutputDirectory("TestVertexObstruction");
        simulator.SetDt(0.001);
        simulator.SetEndTime(100.0);
        simulator.SetSamplingTimestepMultiple(1000);
        simulator.SetOutputNodeVelocities(true);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaMPhaseGrowthForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);          // lambda
        p_force->SetMatureCellTargetArea(1.0);                           // A_0
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);       // beta
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(5.0);      // gamma_cell
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0); // gamma_boundary
        simulator.AddForce(p_force);

        // Create cell killer for y = 12 and pass to the simulation
        c_vector<double, 2> point1 = zero_vector<double>(2);
        point1(1) = 12.0;
        c_vector<double, 2> normal1 = zero_vector<double>(2);
        normal1(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&population, point1, normal1));
        simulator.AddCellKiller(p_killer);

        // Create boundary condition y > 0
        c_vector<double, 2> point2 = zero_vector<double>(2);
        c_vector<double, 2> normal2 = zero_vector<double>(2);
        normal2(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_boundary_condition_1, (&population, point2, normal2));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_1);

        // Create obstruction boundary condition
        MAKE_PTR_ARGS(ObstructionBoundaryCondition<2>, p_boundary_condition_2, (&population));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_2);

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(population);
        WntConcentration<2>::Instance()->SetCryptLength(3); // So only cells at very bottom divide
//        WntConcentration<2>::Instance()->SetCryptLength(40); // So all cells  divide

        // Run simulation
        simulator.Solve();

        // Tidy up
        WntConcentration<2>::Instance()->Destroy();
    }
};

#endif /*TESTVERTEXRESTRICTEDGEOMETRY_HPP_*/
