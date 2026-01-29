/*---------------------------------------------------------------------------*\
Copyright (C) 2015 Cyrille Bonamy, Julien Chauchat, Tian-Jian Hsu
                   and contributors

License
    This file is part of SedFOAM.

    SedFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SedFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with SedFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "twophasekEpsilonVeg2.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void twophasekEpsilonVeg2<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_.min(nutMax_);
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    this->nud_ = Clambda_*sqr(kd_)/epsilond_;
    this->nud_.min(nutMax_);
    this->nud_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nud_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
twophasekEpsilonVeg2<BasicTurbulenceModel>::twophasekEpsilonVeg2
(
    const alphaField& beta,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& betaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        beta,
        rho,
        U,
        betaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    writeTke_
    (
        Switch::getOrAddToDict
        (
            "writeTke",
            this->coeffDict_,
            false
        )
    ),
    C3ep_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C3ep",
            this->coeffDict_,
            1.2
        )
    ),
    C4ep_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C4ep",
            this->coeffDict_,
            1.0
        )
    ),
    KE2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "KE2",
            this->coeffDict_,
            1.0
        )
    ),
    KE4_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "KE4",
            this->coeffDict_,
            1.0
        )
    ),
    KE6_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "KE6",
            this->coeffDict_,
            0.2
        )
    ),
    KE7_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "KE7",
            this->coeffDict_,
            0.15
        )
    ),
    KE8_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "KE8",
            this->coeffDict_,
            1
        )
    ),
    Cmu_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    Clambda_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Clambda",
            this->coeffDict_,
            0.01
        )
    ),
    C1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
        )
    ),
    nutMax_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "nutMax",
            this->coeffDict_,
            1e-1
        )
    ),
    alphak_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphak",
            this->coeffDict_,
            1.0
        )
    ),
    alphaEps_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaEps",
            this->coeffDict_,
            1.3
        )
    ),
    tmfexp_(U.db().lookupObject<volScalarField> ("tmfexp")),
    ESD3_(U.db().lookupObject<volScalarField> ("ESD3")),
    ESD4_(U.db().lookupObject<volScalarField> ("ESD4")),
    ESD5_(U.db().lookupObject<volScalarField> ("ESD5")),
    ESD6_(U.db().lookupObject<volScalarField> ("ESD6")),
    ESD7_(U.db().lookupObject<volScalarField> ("ESD7")),
    ESD8_(U.db().lookupObject<volScalarField> ("ESD8")),
    ESD_(U.db().lookupObject<volScalarField> ("ESD")),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    kt_
    (
        IOobject
        (
            IOobject::groupName("kt", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    kd_
    (
        IOobject
        (
            IOobject::groupName("kd", U.group()),
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
            IOobject::groupName("epsilon", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilont_
    (
        IOobject
        (
            IOobject::groupName("epsilont", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilonw_
    (
        IOobject
        (
            IOobject::groupName("epsilonw", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilond_
    (
        IOobject
        (
            IOobject::groupName("epsilond", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    nud_
    (
        IOobject
        (
            IOobject::groupName("nud", U.group()),
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
bool twophasekEpsilonVeg2<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict_);
        Clambda_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict_);
        C2_.readIfPresent(this->coeffDict_);
        alphak_.readIfPresent(this->coeffDict_);
        alphaEps_.readIfPresent(this->coeffDict_);
        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void twophasekEpsilonVeg2<BasicTurbulenceModel>::correct()
{
    if (not this->turbulence_)
    {
        return;
    }

    // Local references
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    const surfaceScalarField& phi = this->phi_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    volTensorField GradU(fvc::grad(U));
    volSymmTensorField Sij(symm(GradU));

    //volScalarField::Internal G
    volScalarField G
    (
        this->GName(),
        nut*(GradU && dev(twoSymm(GradU)))
    );

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsilontEqn
    (
        fvm::ddt(epsilont_)
      + fvm::div(phi, epsilont_)
      - fvm::Sp(fvc::div(phi), epsilont_)
      - fvm::laplacian(DepsilonEff(), epsilont_)
     ==
      - fvm::SuSp(-C1_*G/kt_, epsilont_)
      + fvm::Sp(-C2_*epsilont_/kt_, epsilont_)
      + fvm::Sp(C3ep_*ESD_, epsilont_)
      + ESD2()*fvm::Sp(C3ep_*KE2_, epsilont_)
      + fvm::Sp(C4ep_*KE4_*ESD5_*nut/kt_, epsilont_)
    );
    if (writeTke_)
    {
        #include "writeTKEBudget_kEpsilonVeg2.H"
    }
    epsilontEqn.ref().relax();
    fvOptions.constrain(epsilontEqn.ref());
    epsilontEqn.ref().boundaryManipulate(epsilont_.boundaryFieldRef());
    solve(epsilontEqn);
    fvOptions.correct(epsilont_);
    bound(epsilont_, this->epsilonMin_);

    epsilonw_ = 0.55*ESD7_*sqrt(kt_)*kt_;
    
    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> ktEqn
    (
        fvm::ddt(kt_)
      + fvm::div(phi, kt_)
      - fvm::Sp(fvc::div(phi), kt_)
      - fvm::laplacian(DkEff(), kt_)
     ==
      - fvm::SuSp(-G/kt_, kt_)
      - fvm::Sp(epsilont_/kt_, kt_)
      + fvm::Sp(ESD_, kt_)
      + fvm::Sp(KE4_*ESD4_*nut/kt_, kt_)
      + ESD2()*fvm::Sp(KE2_, kt_)
      // Transfer from dispersive to turbulent fluctuating energy
      + fvm::Sp(7*nut*pow(ESD7_, 2)*kd_/kt_, kt_)
      // Dissipation in the wake
      + fvm::Sp(-0.55*ESD7_*sqrt(kt_), kt_)
      //Spectral shortcut
      //- fvc::Sp(KE8_*ESD8_, kt_)
    );

    ktEqn.ref().relax();
    fvOptions.constrain(ktEqn.ref());
    solve(ktEqn);
    fvOptions.correct(kt_);
    bound(kt_, this->kMin_);

    // Wake dissipation
    //epsilond_ = KE7_*ESD7_*kd_*sqrt(kd_);
    epsilond_ = 0.4*ESD7_*kd_*sqrt(kd_);
    // Wake turbulent kinetic energy equation
    tmp<fvScalarMatrix> kdEqn
    (
    fvm::ddt(kd_)
      + fvm::div(phi, kd_)
      - fvm::Sp(fvc::div(phi), kd_)
      - fvm::laplacian(DkdEff(), kd_)
      ==
      //turbulent production by vegetation drag force
      - fvm::SuSp(-KE6_*ESD6_/kd_, kd_)
      // Transfer from dispersive to turbulent fluctuating energy
      + fvm::Sp(-7*nut*pow(ESD7_, 2), kd_)
      // Dissipation in the wake
      + fvm::Sp(-0.4*ESD7_*sqrt(kd_), kd_)
      //Spectral shortcut
      //+ fvc::Sp(KE8_*ESD8_, kt_)
    );

    kdEqn.ref().relax();
    fvOptions.constrain(kdEqn.ref());
    solve(kdEqn);
    fvOptions.correct(kd_);
    bound(kd_, this->kMin_);

    k_ = kt_;
    epsilon_ = epsilont_ + epsilonw_;

    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
