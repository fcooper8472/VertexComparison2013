/*

Copyright (c) 2005-2019, University of Oxford.
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

/**********************************************
 * THIS CODE WORKS WITH RELEASE 3.3 OF CHASTE *
 **********************************************/

#ifndef SIMPLEVOLUMEBASEDSTOCHASTICCELLCYCLEMODEL_HPP_
#define SIMPLEVOLUMEBASEDSTOCHASTICCELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"

/**
 * Simple area-based  cell cycle model with a Uniform phase duration with a default distribution of U[10,14].
 *
 * The cell cycle model inherits functionality from AbstractSimpleCellCycleModel. In addition
 * it imposes a form of 'contact inhibition' whereby a cell (not in M phase), whose area
 * drops below some critical fraction of its equilibrium value, temporarily stops progressing
 * through the cell cycle.
 */
class SimpleVolumeBasedStochasticCellCycleModel : public AbstractSimpleCellCycleModel
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);
        archive & mQuiescentVolumeFraction;
    }

    /**
     * The fraction of a cell's equilibrium area below which the cell
     * temporarily stops progressing through the cell cycle.
     */
    double mQuiescentVolumeFraction;

    /**
     * Stochastically set the G1 duration.  Called on cell creation at
     * the start of a simulation, and for both parent and daughter
     * cells at cell division.
     */
    void SetG1Duration();

public:

    /**
     * Default constructor.
     */
    SimpleVolumeBasedStochasticCellCycleModel();

    /**
     * Overridden UpdateCellCyclePhase() method.
     */
    void UpdateCellCyclePhase();

    /**
     * Overridden builder method to create new instances of
     * the cell cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Set #mQuiescentVolumeFraction.
     *
     * @param quiescentVolumeFraction
     */
    void SetQuiescentVolumeFraction(double quiescentVolumeFraction);

    /**
     * Get #mQuiescentVolumeFraction.
     */
    double GetQuiescentVolumeFraction();
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(SimpleVolumeBasedStochasticCellCycleModel)

#endif /*SIMPLEVOLUMEBASEDSTOCHASTICCELLCYCLEMODEL_HPP_*/
