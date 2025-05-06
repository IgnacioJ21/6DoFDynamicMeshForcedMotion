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

\*---------------------------------------------------------------------------*/

#include "forcedTimeSeriesRigidBodyMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(forcedTimeSeriesRigidBodyMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        forcedTimeSeriesRigidBodyMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::forcedTimeSeriesRigidBodyMotionSolver::forcedTimeSeriesRigidBodyMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    /*motion_
    (
        coeffDict(),
        IOobject
        (
            "forcedRigidBodyMotionState",
            mesh.time().timeName(),
            "uniform",
            mesh
        ).typeHeaderOk<IOdictionary>(true)
      ? IOdictionary
        (
            IOobject
            (
                "forcedRigidBodyMotionState",
                mesh.time().timeName(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        )
      : coeffDict(),
        mesh.time()
    ),*/
    centreOfRotation0_(coeffDict().get<vector>("initialCoR")),
    centreOfRotation(vector(0,0,0)),
    Q0_(),
    Rz0_(),
    Ry0_(),
    Rx0_(),
    patches_(coeffDict().get<wordRes>("patches")),
    patchSet_(mesh.boundaryMesh().patchSet(patches_)),
    di_(coeffDict().get<scalar>("innerDistance")),
    do_(coeffDict().get<scalar>("outerDistance")),
    a_(coeffDict().get<vector>("translationalAmplitude")),
    aR_(coeffDict().get<vector>("rotationalAmplitude")),
    w_(coeffDict().get<vector>("translationalFrequency")),
    wR_(coeffDict().get<vector>("rotationalFrequency")),
    alpha0_(coeffDict().get<vector>("initialAngle")),
    //timeList_(coeffDict().get<List<scalar>>("timeList")), // Needs testing by Ignacio
    //transRotList_(coeffDict().get<List<symmTensor>>("transRotList")), // Needs testing by Ignacio
    test_(coeffDict().getOrDefault("test", false)),
    rhoInf_(1.0),
    rhoName_(coeffDict().getOrDefault<word>("rho", "rho")),
    scale_
    (
        IOobject
        (
            "motionScale",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, Zero)
    ),
    curTimeIndex_(-1),
    cOfGdisplacement_(coeffDict().getOrDefault<word>("cOfGdisplacement", "none"))
{

    Rz0_[0] = cos(alpha0_[2]);
    Rz0_[1] = -sin(alpha0_[2]);
    Rz0_[2] = 0;
    Rz0_[3] = sin(alpha0_[2]);
    Rz0_[4] = cos(alpha0_[2]);
    Rz0_[5] = 0;
    Rz0_[6] = 0;
    Rz0_[7] = 0;
    Rz0_[8] = 1;
    
    Ry0_[0] = cos(alpha0_[1]);
    Ry0_[1] = 0;
    Ry0_[2] = sin(alpha0_[1]);
    Ry0_[3] = 0;
    Ry0_[4] = 1;
    Ry0_[5] = 0;
    Ry0_[6] = -sin(alpha0_[1]);
    Ry0_[7] = 0;
    Ry0_[8] = cos(alpha0_[1]);
    
    Rx0_[0] = 1;
    Rx0_[1] = 0;
    Rx0_[2] = 0;
    Rx0_[3] = 0;
    Rx0_[4] = cos(alpha0_[0]);
    Rx0_[5] = -sin(alpha0_[0]);
    Rx0_[6] = 0;
    Rx0_[7] = sin(alpha0_[0]);
    Rx0_[8] = cos(alpha0_[0]);
    
       
    Q0_=(Rz0_ & Ry0_ & Rx0_);


    if (rhoName_ == "rhoInf")
    {
        coeffDict().readEntry("rhoInf", rhoInf_);
    }

    // Calculate scaling factor everywhere
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const pointMesh& pMesh = pointMesh::New(mesh);

        pointPatchDist pDist(pMesh, patchSet_, points0());

        // Scaling: 1 up to di then linear down to 0 at do away from patches
        scale_.primitiveFieldRef() =
            min
            (
                max
                (
                    (do_ - pDist.primitiveField())/(do_ - di_),
                    scalar(0)
                ),
                scalar(1)
            );

        // Convert the scale function to a cosine
        scale_.primitiveFieldRef() =
            min
            (
                max
                (
                    0.5
                  - 0.5
                   *cos(scale_.primitiveField()
                   *Foam::constant::mathematical::pi),
                    scalar(0)
                ),
                scalar(1)
            );

        pointConstraints::New(pMesh).constrain(scale_);
        scale_.write();
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*const Foam::forcedRigidBodyMotion&
Foam::forcedTimeSeriesRigidBodyMotionSolver::motion() const
{
    return motion_;
}*/


Foam::tmp<Foam::pointField>
Foam::forcedTimeSeriesRigidBodyMotionSolver::curPoints() const
{
    tmp<pointField> newPoints
    (
        points0() + pointDisplacement_.primitiveField()
    );

    if (!moveAllCells())
    {
        tmp<pointField> ttransformedPts(new pointField(mesh().points()));
        pointField& transformedPts = ttransformedPts.ref();

        UIndirectList<point>(transformedPts, pointIDs()) =
            pointField(newPoints.ref(), pointIDs());

        return ttransformedPts;
    }

    return newPoints;
}


void Foam::forcedTimeSeriesRigidBodyMotionSolver::solve()
{
    const Time& t = mesh().time();

    if (mesh().nPoints() != points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << points0().size()
            << " points; in the current mesh there are " << mesh().nPoints()
            << " points." << exit(FatalError);
    }

    // Store the motion state at the beginning of the time-step
    bool firstIter = false;
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        //motion_.newTime();
        curTimeIndex_ = this->db().time().timeIndex();
        firstIter = true;
    }

    dimensionedVector g("g", dimAcceleration, Zero);

    if (db().time().foundObject<uniformDimensionedVectorField>("g"))
    {
        g = db().time().lookupObject<uniformDimensionedVectorField>("g");
    }
    else
    {
        coeffDict().readIfPresent("g", g);
    }

    // const scalar ramp = min(max((this->db().time().value() - 5)/10, 0), 1);
    const scalar ramp = 1.0;

    if (test_)
    {
        /*motion_.update
        (
            firstIter,
            ramp*(motion_.mass()*g.value()),
            ramp*(motion_.mass()*(motion_.momentArm() ^ g.value())),
            t.deltaTValue(),
            t.deltaT0Value()
        );*/
    }
    else
    {
        dictionary forcesDict;

        forcesDict.add("type", functionObjects::forces::typeName);
        forcesDict.add("patches", patches_);
        forcesDict.add("rhoInf", rhoInf_);
        forcesDict.add("rho", rhoName_);
        //forcesDict.add("CofR", motion_.centreOfRotation());
        forcesDict.add("CofR", centreOfRotation);

        vector oldPos = centreOfRotation;

        functionObjects::forces f("forces", db(), forcesDict);

        f.calcForcesMoments();

        /*motion_.update
        (
            firstIter,
            ramp*(f.forceEff() + motion_.mass()*g.value()),
            ramp
           *(
               f.momentEff()
             + motion_.mass()*(motion_.momentArm() ^ g.value())
            ),
            t.deltaTValue(),
            t.deltaT0Value()
        );*/

        /* if (cOfGdisplacement_ != "none")
        {
            if
            (
                db().time().foundObject<uniformDimensionedVectorField>
                (
                    cOfGdisplacement_
                )
            )
            {
                auto& disp =
                    db().time().lookupObjectRef<uniformDimensionedVectorField>
                    (
                        cOfGdisplacement_
                    );

                disp += (motion_.centreOfRotation() - oldPos);
            }
        } */
    }

    // Update the displacements
    pointDisplacement_.primitiveFieldRef() = transform(points0(), scale_, t) - points0();
        //motion_.transform(points0(), scale_) - points0();

    // Displacement has changed. Update boundary conditions
    pointConstraints::New
    (
        pointDisplacement_.mesh()
    ).constrainDisplacement(pointDisplacement_);
}


/*bool Foam::forcedTimeSeriesRigidBodyMotionSolver::writeObject
(
    IOstreamOption streamOpt,
    const bool valid
) const
{
    IOdictionary dict
    (
        IOobject
        (
            "forcedRigidBodyMotionState",
            mesh().time().timeName(),
            "uniform",
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    motion_.state().write(dict);
    return dict.regIOobject::write();
}
*/

bool Foam::forcedTimeSeriesRigidBodyMotionSolver::read()
{
    if (displacementMotionSolver::read())
    {
        //motion_.read(coeffDict());

        return true;
    }

    return false;
}

Foam::tmp<Foam::pointField> Foam::forcedTimeSeriesRigidBodyMotionSolver::transform
(
    const pointField& initialPoints,
    const scalarField& scale,
    const Time& t
) 
{

    //scalar a = 0.1;
    //scalar w = 3.14;
    
    //Added this interpolation for time series input. Needs testing by Ignacio
    /*int i = 0;
    while (t.value() < timeList_[i])
    {
        i++;
    }  
    
    scalar t0 = timeList_[i-1];
    scalar t1 = timeList_[i];
    symmTensor tR0 = transRotList_[i-1];
    symmTensor tR1 = transRotList_[i];    
    symmTensor currentPosition = tR0 + ((tR1 - tR0)/(t1 - t0)) * (t.value() - t0);*/

    vector centreOfRotation = centreOfRotation0_ + vector(a_[0]*Foam::sin(w_[0]*t.value()),a_[1]*Foam::sin(w_[1]*t.value()),a_[2]*Foam::sin(w_[2]*t.value())); //centreOfRotation0 + vector(0.001,0,0);
    //scalar beta0 = 0;
    scalar alpha = alpha0_[0] + aR_[0]*Foam::sin(wR_[0]*t.value());
    scalar beta = alpha0_[1] + aR_[1]*Foam::sin(wR_[1]*t.value());
    scalar gamma = alpha0_[2] + aR_[2]*Foam::sin(wR_[2]*t.value());
    
    tensor Rz;
    tensor Ry;
    tensor Rx;
    
    Rz[0] = cos(gamma);
    Rz[1] = -sin(gamma);
    Rz[2] = 0;
    Rz[3] = sin(gamma);
    Rz[4] = cos(gamma);
    Rz[5] = 0;
    Rz[6] = 0;
    Rz[7] = 0;
    Rz[8] = 1;
    
    Ry[0] = cos(beta);
    Ry[1] = 0;
    Ry[2] = sin(beta);
    Ry[3] = 0;
    Ry[4] = 1;
    Ry[5] = 0;
    Ry[6] = -sin(beta);
    Ry[7] = 0;
    Ry[8] = cos(beta);
    
    Rx[0] = 1;
    Rx[1] = 0;
    Rx[2] = 0;
    Rx[3] = 0;
    Rx[4] = cos(alpha);
    Rx[5] = -sin(alpha);
    Rx[6] = 0;
    Rx[7] = sin(alpha);
    Rx[8] = cos(alpha);
    
    tensor Q; // instantaneous rotation
    Q=(Rz & Ry & Rx);
    
    Info << "    CentreOfRotation: " << centreOfRotation << endl;
    Info << "    OrientationMatrix: " << Q << endl;
        // Calculate the transformation septerion from the initial state
    septernion s
    (
        centreOfRotation - centreOfRotation0_,
        quaternion(Q.T() & Q0_)
    );

    tmp<pointField> tpoints(new pointField(initialPoints));
    pointField& points = tpoints.ref();

    forAll(points, pointi)
    {
        // Move non-stationary points
        if (scale[pointi] > SMALL)
        {
            // Use solid-body motion where scale = 1
            if (scale[pointi] > 1 - SMALL)
            {
                points[pointi] = centreOfRotation + (Q & Q0_.T() & (initialPoints[pointi] - centreOfRotation0_)); //transform(initialPoints[pointi]);
            }
            // Slerp septernion interpolation
            else
            {
                septernion ss(slerp(septernion::I, s, scale[pointi]));

                points[pointi] =
                    centreOfRotation0_
                  + ss.invTransformPoint
                    (
                        initialPoints[pointi]
                      - centreOfRotation0_
                    );
            }
        }
    }

    return tpoints;
}


// ************************************************************************* //
