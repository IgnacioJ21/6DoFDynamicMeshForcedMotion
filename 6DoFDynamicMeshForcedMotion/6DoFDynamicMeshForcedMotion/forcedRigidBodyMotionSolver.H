/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::forcedRigidBodyMotionSolver

Group
   grpMeshMotionSolvers

Description
    6-DoF solid-body mesh motion solver for an fvMesh.

    Applies SLERP interpolation of movement as function of distance to
    the object surface.

SourceFiles
    forcedRigidBodyMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef forcedRigidBodyMotionSolver_H
#define forcedRigidBodyMotionSolver_H

#include "displacementMotionSolver.H"
//#include "forcedRigidBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class forcedRigidBodyMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class forcedRigidBodyMotionSolver
:
    public displacementMotionSolver
{
    // Private data

        //- Six DoF motion object
        //forcedRigidBodyMotion motion_;
        
        vector centreOfRotation;
        vector centreOfRotation0_;
        tensor Q0_;
        tensor Rz0_;
        tensor Ry0_;
        tensor Rx0_; //Initial rotation (probably 0 degrees = identity matrix)
        vector a_;
        vector aR_;
        vector w_;
        vector wR_;
        vector alpha0_;
        
       // List<scalar> timeList_; // List of times to be interpolated // Needs testing by Ignacio
       // List<symmTensor> transRotList_; // List of translations and rotations to be interpolated // Needs testing by Ignacio
        
        wordRes patches_;

        //- Patches to integrate forces
        const labelHashSet patchSet_;

        //- Inner morphing distance (limit of solid-body region)
        const scalar di_;

        //- Outer morphing distance (limit of linear interpolation region)
        const scalar do_;

        //- Test-mode in which only the gravitational body-force is applied
        bool test_;

        //- Reference density required by the forces object for
        //  incompressible calculations, required if rho == rhoInf
        scalar rhoInf_;

        //- Name of density field, optional unless used for an
        //  incompressible simulation, when this needs to be specified
        //  as rhoInf
        word rhoName_;

	// run up factor, its the time in seconds
	scalar tSoft_;

        //- Current interpolation scale (1 at patches, 0 at distance_)
        pointScalarField scale_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Name of the uniformVectorField for CofG displacement
        word cOfGdisplacement_;


    // Private Member Functions

        //- No copy construct
        forcedRigidBodyMotionSolver
        (
            const forcedRigidBodyMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const forcedRigidBodyMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("forcedRigidBodyMotion");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        forcedRigidBodyMotionSolver
        (
            const polyMesh&,
            const IOdictionary& dict
        );


    //- Destructor
    ~forcedRigidBodyMotionSolver() = default;


    // Member Functions

        //- Return the six DoF motion object
        //const forcedRigidBodyMotion& motion() const;

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Write state using stream options
        /*virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool valid
        ) const;*/

        //- Read dynamicMeshDict dictionary
        virtual bool read();
        
            tmp<pointField> transform
            (
                const pointField& initialPoints,
                const scalarField& scale,
                const Time& t
            );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
