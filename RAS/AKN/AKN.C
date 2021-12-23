/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "AKN.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> AKN<BasicTurbulenceModel>::yStar() const
{
    return pow(this->nu()*epsilon_,scalar(0.25))*y_/this->nu();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> AKN<BasicTurbulenceModel>::Rt() const
{
    return sqr(k_)/(this->nu()*epsilon_);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> AKN<BasicTurbulenceModel>::fMu() const
{
    return sqr(scalar(1) - exp(-yStar()_/14.0))
	  *(scalar(1) + 5.0/pow(Rt() + SMALL,0.57)*exp(-sqr(Rt()/200.0)));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> AKN<BasicTurbulenceModel>::f2() const
{
    return (scalar(1) - 0.3*exp(-sqr(Rt()/0.5)))
	  *sqr(scalar(1) - exp(-yStar()/3.1));
}


template<class BasicTurbulenceModel>
void AKN<BasicTurbulenceModel>::correctNut() const
{
    this->nut_ = Cmu_*fMu()*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> AKN<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
	new fvScalarMatrix
	(
	    k_,
	    dimVolume*this->rho_.dimensions()*k_.dimensions()
	    /dimTime
	)
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> AKN<BasicTurbulenceModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
	new fvScalarMatrix
	(
	    epsilon_,
	    dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
	    /dimTime
	)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
AKN<BasicTurbulenceModel>::AKN
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceSclaarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
	type,
	alpha,
	rho,
	U,
	alphaRhoPhi,
	phi,
	transport,
	propertiesName
    ),

    Cmu_
    (
	dimensioned<scalar>::lookupOrAddToDict
	(
	    "Cmu",
	    coeffDict_,
	    0.09
	)
    ),

    C1_
    (
	dimensioned<scalar>::lookupOrAddToDict
	(
	    "C1",
	    coeffDict_,
	    1.5
	)
    ),

    C2_
    (
	dimensioned<scalar>::lookupOrAddToDict
	(
	    "C2",
	    coeffDict_,
	    1.9
	)
    ),

    sigmaK_
    (
	dimensioned<scalar>::lookupOrAddToDict
	(
	    "sigmaK",
	    coeffDict_,
	    1.4
	)
    ),

    sigmaEps_
    (
	dimensioned<scalar>::lookupOrAddToDict
	(
	    "sigmaEps",
	    coeffDict_,
	    1.4
	)
    ),

    k_
    (
	IOobject
	(
	    "k",
	    this->runTime_.timeName(),
	    this->mesh_,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	),
	this->mesh_
    ),

    epsilon_
    (
	IOobject
	(
	    "epsilon",
	    this->runTime_.timeName(),
	    this->mesh_,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	),
	this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
	this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool AKN<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
	Cmu_.readIfPresent(this->coeffDict());
	C1_.readIfPresent(this->coeffDict());
	C2_.readIfPresent(this->coeffDict());
	sigmaK_.readIfPresent(this->coeffDict());
	sigmaEps_.readIfPresent(this->coeffDict());

	return true;
    }
    else
    {
	return false;
    }
}


template<class BasicTurbulenceModel>
void AKN<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
	return;
    }

    // Local reference
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const volVectorField& U = this->U_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField::Internal divU
    (
	fvc::div(fvc::absolute(this->phi(), U))().v()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);

    volScalarField::Internal G
    (
	this->Gname(),
	nut.v()*(dev(twoSymm(tradU().v())) && tgradU().v())
    );
    tgradU.clear();

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
	fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
	C1_*alpha()*rho()*G*epsilon_()/k_()
      - fvm::Sp(C2_*f2()*alpha()*rho()*epsilon_/k_, epsilon_)
      + epsilonSource()
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
	fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvM;;laplacian(alpha*rho*DkEff(), k_)
     ==
	alpha()*rho()*G
      - fvm::Sp(alpha()*rho()*epsilon_()/k_(), k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
