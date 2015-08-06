/*

Copyright (c) 2005-2015, University of Oxford.
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

#ifndef TESTVERTEXMONOLAYERS_HPP_
#define TESTVERTEXMONOLAYERS_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "VolumeTrackingModifier.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "SimpleVolumeBasedStochasticCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "NagaiHondaForce.hpp"
#include "MPhaseGrowthTargetAreaModifier.hpp"
#include "WelikyOsterForce.hpp"
#include "ModifiedWelikyOsterForce.hpp"
#include "VertexAngleForce.hpp"
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/**********************************************
 * THIS CODE WORKS WITH RELEASE 3.3 OF CHASTE *
 **********************************************/

/**
 * Test suite defining a vertex dynamics simulation of a growing monolayer.
 */
class TestVertexMonolayers: public AbstractCellBasedTestSuite
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

    /**
     * Helper method. Populates a vector of cells of specified size.
     * Each cell has a simple stochastic area-dependent cell-cycle
     * model: it enters quiescence if its area drops below a specified
     * threshold. 
     */
    void SetUpCellsWithStochasticAreaDependentCellCycleModel(std::vector<CellPtr>& rCells, unsigned numCells, double quiescentVolumeFraction)
    {
        rCells.clear();
        rCells.reserve(numCells);

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        for (unsigned i=0; i<numCells; i++)
        {
            // Create a cell-cycle model and set contact inhibition parameter
        	SimpleVolumeBasedStochasticCellCycleModel* p_cell_cycle_model = new SimpleVolumeBasedStochasticCellCycleModel;
            p_cell_cycle_model->SetDimension(2);
            p_cell_cycle_model->SetQuiescentVolumeFraction(quiescentVolumeFraction);
            p_cell_cycle_model->SetSDuration(1e-10);
            p_cell_cycle_model->SetG2Duration(1e-10);

            // Create a cell using the cell-cycle model and a 'wild-type' cell mutation state
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);

            // Give the cell a random birth time
            double birth_time = -p_cell_cycle_model->GetAverageTransitCellCycleTime()*RandomNumberGenerator::Instance()->ranf();
            p_cell->SetBirthTime(birth_time);

            rCells.push_back(p_cell);
        }
    }

public:

    /**
     * Simulate a monolayer using the vertex dynamics model proposed by
     * T. Nagai and H. Honda ("A dynamic cell model for the formation of
     * epithelial tissues", Philosophical Magazine Part B 81:699-719).
     * 
     * Each of the vertex dynamics model parameter member variables are
     * rescaled such that mDampingConstantNormal takes the default value 1,
     * whereas Nagai and Honda (who denote the parameter by nu) take the
     * value 0.01.
     */
    void TestNagaiHondaMonolayer() throw (Exception)
    {
        // Create a simple 2D vertex-based mesh with a single element
        HoneycombVertexMeshGenerator mesh_generator(1, 1);
        MutableVertexMesh<2,2>* p_mesh = mesh_generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        // Set up a vector of cells
        std::vector<CellPtr> cells;

        // Create 1 cell, with stochastic area-dependent cell-cycle model with phi=0.9 as 0.812261 is relaxed area of cell
        SetUpCellsWithStochasticAreaDependentCellCycleModel(cells, 1, 0.812261*0.9);

        // Create a vertex-based cell population
        VertexBasedCellPopulation<2> population(*p_mesh, cells);

        // Specify what to output from the simulation
        population.AddCellWriter<CellAgesWriter>();
        population.AddCellWriter<CellAncestorWriter>();
        population.AddCellWriter<CellIdWriter>();
        population.AddCellWriter<CellMutationStatesWriter>();
        population.AddCellWriter<CellProliferativePhasesWriter>();
        population.AddCellWriter<CellProliferativeTypesWriter>();
        population.AddCellWriter<CellVolumesWriter>();

        // Set up a cell-based simulation, output directory, time step and end time
        OffLatticeSimulation<2> simulator(population);
        simulator.SetOutputDirectory("TestMonolayer/NagaiHonda");
        simulator.SetDt(0.001);
        simulator.SetEndTime(200);

        // Only record results every 1000 time steps
        simulator.SetSamplingTimestepMultiple(1000);

        // Add volume-tracking modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55);          // lambda
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);       // beta
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(5.0);      // gamma_cell
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0); // gamma_boundary
        simulator.AddForce(p_force);

        MAKE_PTR(MPhaseGrowthTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulator.Solve();
    }

    /**
     * Simulate a monolayer using a modified version of the vertex dynamics
     * model proposed by M. Weliky and G. Oster ("The mechanical basis of
     * cell rearrangement. I. Epithelial morphogenesis during Fundulus epiboly",
     * Development 109:373-386), including an additional membrane stiffness
     * force.
     *
     * The default values for the two model parameter member variables are
     * our own best estimates, since they are not given in the Weliky and
     * Oster paper.
     */
    void TestModifiedWelikyOsterMonolayerWithVertexAngleForce() throw (Exception)
    {
//        // Create a simple 2D vertex-based mesh with a single element
//        HoneycombVertexMeshGenerator mesh_generator(1, 1);
//        MutableVertexMesh<2,2>* p_mesh = mesh_generator.GetMesh();
//        p_mesh->SetCellRearrangementThreshold(0.1);
//
//        // Set up a vector of cells
//        std::vector<CellPtr> cells;
//
//        // Create 1 cell, with stochastic area-dependent cell-cycle model with phi=0.9 as 0.812261 is relaxed area of cell
//        SetUpCellsWithStochasticAreaDependentCellCycleModel(cells, 1, 0.812261*0.9);
//
//        // Create a vertex-based cell population
//        VertexBasedCellPopulation<2> population(*p_mesh, cells);
//
//        // Set population to output all data to results files
//        population.SetOutputCellIdData(true);
//        population.SetOutputCellMutationStates(true);
//        population.SetOutputCellAncestors(true);
//        population.SetOutputCellProliferativeTypes(true);
//        population.SetOutputCellVariables(true);
//        population.SetOutputCellCyclePhases(true);
//        population.SetOutputCellAges(true);
//        population.SetOutputCellVolumes(true);
//
//        // Set up a cell-based simulation, output directory, time step and end time
//        OffLatticeSimulation<2> simulator(population);
//        simulator.SetOutputDirectory("TestMonolayer/ModifiedWelikyOsterWithVertexAngleForce");
//        simulator.SetDt(0.001);
//        simulator.SetEndTime(115);
//
//        // Only record results every 1000 time steps
//        simulator.SetSamplingTimestepMultiple(2000);
//
//        // Add volume-tracking modifier
//        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
//        simulator.AddSimulationModifier(p_modifier);
//
//        // Create a force law and pass it to the simulation
//        MAKE_PTR(ModifiedWelikyOsterForce<2>, p_force);
//        p_force->SetWelikyOsterAreaParameter(1.0);      // beta
//        p_force->SetWelikyOsterPerimeterParameter(1.0); // kappa
//        simulator.AddForce(p_force);
//
//        // Create another force law and pass it to the simulation
//        MAKE_PTR(VertexAngleForce<2>, p_force2);
//        p_force2->SetAngleRestrainingParameter(10.0); // gamma
//        simulator.AddForce(p_force2);
//
//        // Run simulation
//        simulator.Solve();
//
//        CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);
//
//        OffLatticeSimulation<2>* p_simulator1;
//        p_simulator1 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load("TestMonolayer/ModifiedWelikyOsterWithVertexAngleForce", 115);
//        p_simulator1->SetDt(0.0005);
//        p_simulator1->SetEndTime(140);
//        p_simulator1->Solve();
//        CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(p_simulator1);
//
//        OffLatticeSimulation<2>* p_simulator2;
//        p_simulator2 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load("TestMonolayer/ModifiedWelikyOsterWithVertexAngleForce", 140);
//        p_simulator2->SetDt(0.0001);
//        p_simulator2->SetEndTime(160);
//        p_simulator2->Solve();
//        CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(p_simulator2);

        OffLatticeSimulation<2>* p_simulator3;
        p_simulator3 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load("TestMonolayer/ModifiedWelikyOsterWithVertexAngleForce", 160);
        p_simulator3->SetDt(0.001);
        p_simulator3->SetEndTime(200);
        p_simulator3->Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(p_simulator3);
    }
};

#endif /*TESTVERTEXMONOLAYERS_HPP_*/
