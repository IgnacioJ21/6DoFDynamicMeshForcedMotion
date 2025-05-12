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

#include "forcedRigidBodyMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "mathematicalConstants.H"
#include <fstream>
#include "error.H" 
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(forcedRigidBodyMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        forcedRigidBodyMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::forcedRigidBodyMotionSolver::forcedRigidBodyMotionSolver
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
    FileDisplacement_(vector(0,0,0)),
    Q0_(),
    Rz0_(),
    Ry0_(),
    Rx0_(),
    patches_(coeffDict().get<wordRes>("patches")),
    motionFromFileDictName_("ForcedMotionDict"),
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
    motionFromFile_(coeffDict().getOrDefault("motionFromFile", false)),
    rhoInf_(1.0),
    rhoName_(coeffDict().getOrDefault<word>("rho", "rho")),
    tSoft_(coeffDict().getOrDefault<scalar>("tSoft", 0)),  
    FileDisplacementDoF_(coeffDict().getOrDefault<scalar>("FileDisplacementDoF", 0)),
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

    // Read the motion file and store it 
    const Time& runTime = mesh.time();

    if (motionFromFile_)
    {
        IOdictionary dict(
            IOobject(
                "ForcedMotionDict",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        fileName motionFileName(runTime.constant()/ word(dict.lookup("file")));
        readMotionFile(motionFileName);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*const Foam::forcedRigidBodyMotion&
Foam::forcedRigidBodyMotionSolver::motion() const
{
    return motion_;
}*/


Foam::tmp<Foam::pointField>
Foam::forcedRigidBodyMotionSolver::curPoints() const
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

void Foam::forcedRigidBodyMotionSolver::readMotionFile(const fileName& motionFileName)
{
    std::ifstream motionFile(motionFileName.c_str());
    
    if (!motionFile.is_open()) {
        FatalErrorInFunction
            << "Unable to open motion file " << motionFileName
            << Foam::endl;
        Foam::exit(FatalError);
    }
    
    std::string line;
    
    // Skip the first line if it's a header
    std::getline(motionFile, line);

    // Read the rest of the file and store the times and motions
    while (std::getline(motionFile, line)) {
        std::istringstream stream(line);
        scalar time;
        scalar motion;
        
        stream >> time >> motion;  // Read time and motion values
        
        times.append(time);  // Append to the 'times' list
        motions.append(motion);  // Append to the 'motions' list    
    }
    //  Print only once after reading all lines
    Foam::Info << "Motion data loaded from file " << motionFileName << Foam::nl;
    Foam::Info << "Number of data points: " << times.size() << Foam::nl;
}


void Foam::forcedRigidBodyMotionSolver::solve()
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


/*bool Foam::forcedRigidBodyMotionSolver::writeObject
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

bool Foam::forcedRigidBodyMotionSolver::read()
{
    if (displacementMotionSolver::read())
    {
        //motion_.read(coeffDict());

        return true;
    }

    return false;
}

Foam::tmp<Foam::pointField> Foam::forcedRigidBodyMotionSolver::transform
(
    const pointField& initialPoints,
    const scalarField& scale,
    const Time& t
) 
{


    scalar PI_(Foam::constant::mathematical::pi);
    scalar factorRunUp(1.0);
    vector FileDisplacement(0, 0, 0);
    vector centreOfRotation(0, 0, 0);
    scalar alpha, beta, gamma;

    if (tSoft_ > 0.0)
    {
        factorRunUp = (Foam::sin(2*PI_/(2.0*tSoft_)*Foam::min(tSoft_, t.value())+270*PI_/180)+1)/2;

//	factorRunUp = (sin(2*pi/(2.0*tSoft)*min(tSoft,t)+270*pi/180)+1)/2;
    }

    if (motionFromFile_)
    {
        
        scalar currentTime = db().time().value();  // or this->db().time().value()

        // Initialize the motion value
        scalar currentMotion = 0.0;

        // Find the two nearest times surrounding the current time
        label i = 0;
        while (i < times.size() && times[i] < currentTime)
        {
            ++i;
        }

        // If currentTime is beyond the last time point
        if (i == times.size())
        {
            currentMotion = motions[times.size() - 1];  // Use the last motion value
        }
        else if (i == 0)
        {
            currentMotion = motions[0];  // Use the first motion value if currentTime is before the first time
        }
        else
        {
            // Interpolation between times[i-1] and times[i]
            scalar t1 = times[i - 1];
            scalar t2 = times[i];
            scalar m1 = motions[i - 1];
            scalar m2 = motions[i];

            // Linear interpolation formula
            currentMotion = m1 + (m2 - m1) * (currentTime - t1) / (t2 - t1);
        }

        Foam::Info << "At time = " << currentTime 
           << ", currentMotion = " << currentMotion << Foam::endl;

        // Assign to desired DoF
        scalar DoF = FileDisplacementDoF_ - 1;
        Foam::Info << "Current FileDisplacementDOF: " << DoF << Foam::endl;
        FileDisplacement[DoF] = currentMotion;
        Foam::Info << "Current FileDisplacement: " << FileDisplacement << Foam::endl;
        centreOfRotation = centreOfRotation0_ + vector(factorRunUp*a_[0]*Foam::sin(w_[0]*t.value()),factorRunUp*a_[1]*Foam::sin(w_[1]*t.value()),factorRunUp*a_[2]*Foam::sin(w_[2]*t.value())) + FileDisplacement; //centreOfRotation0 + vector(0.001,0,0);
        Foam::Info << "Current centreOfRotation: " << centreOfRotation << Foam::endl;
        //scalar beta0 = 0;
        alpha = alpha0_[0] + factorRunUp*aR_[0]*Foam::sin(wR_[0]*t.value());
        beta = alpha0_[1] + factorRunUp*aR_[1]*Foam::sin(wR_[1]*t.value());
        gamma = alpha0_[2] + factorRunUp*aR_[2]*Foam::sin(wR_[2]*t.value());   
    }
    else
    {
        centreOfRotation = centreOfRotation0_ + vector(factorRunUp*a_[0]*Foam::sin(w_[0]*t.value()),factorRunUp*a_[1]*Foam::sin(w_[1]*t.value()),factorRunUp*a_[2]*Foam::sin(w_[2]*t.value())); 
        //scalar beta0 = 0;
        alpha = alpha0_[0] + factorRunUp*aR_[0]*Foam::sin(wR_[0]*t.value());
        beta = alpha0_[1] + factorRunUp*aR_[1]*Foam::sin(wR_[1]*t.value());
        gamma = alpha0_[2] + factorRunUp*aR_[2]*Foam::sin(wR_[2]*t.value());
    }
        
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
