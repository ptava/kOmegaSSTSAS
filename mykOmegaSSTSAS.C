/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
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

#include "mykOmegaSSTSAS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void mykOmegaSSTSAS<BasicTurbulenceModel>::makeDumpField
(
    autoPtr<volScalarField>& fld,
    const word& name
) const
{
    const bool needCreate = !fld.valid() || (fld->size() != this->mesh_.nCells());

    if (!needCreate) return;

    fld.reset
    (
        new volScalarField
        (
            IOobject
            (
                name,
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar(dimensionSet(0, 1, 0, 0, 0), Zero)
        )
    );
}

template<class BasicTurbulenceModel>
template<Foam::autoPtr<Foam::volScalarField>
         mykOmegaSSTSAS<BasicTurbulenceModel>::* MemberPtr>
inline void mykOmegaSSTSAS<BasicTurbulenceModel>::setDumpField
(
    const scalarField& vals,
    const word& name
) const
{
    // Sanity check
    if (vals.size() != this->mesh_.nCells())
    {
        FatalErrorInFunction
            << "Size mismatch for field '" << name << "': vals.size() = "
            << vals.size() << " vs nCells() = " << this->mesh_.nCells()
            << exit(FatalError);
    }

    // Resolve member pointer to actual autoPtr
    autoPtr<volScalarField>& fld = const_cast<autoPtr<volScalarField>&>(
        this->*MemberPtr
    );

    // Ensure field exists and matches current mesh
    makeDumpField(fld, name);

    // Assign internal values
    fld->primitiveFieldRef() = vals;
}

template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> mykOmegaSSTSAS<BasicTurbulenceModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    volScalarField::Internal L
    (
        sqrt(this->k_())/(pow025(this->betaStar_)*this->omega_())
    );

    // force update of new private members C1 and C2
    volScalarField::Internal C1 =
    (
        kappa_*sqrt(S2)
        /(
            mag(fvc::laplacian(this->U_))()()
            + dimensionedScalar
            (
                "ROOTVSMALL",
                dimensionSet(0, -1, -1, 0, 0),
                ROOTVSMALL
            )
        )
    );

    volScalarField::Internal C2 =
        Cs_*sqrt(kappa_*zeta2_/(beta/this->betaStar_ - gamma))*delta()();

    volScalarField::Internal Lvk = max(C1,C2);

    this->template setDumpField<&mykOmegaSSTSAS<BasicTurbulenceModel>::C1_>
    (
        C1.field(), "C1"
    );
    this->template setDumpField<&mykOmegaSSTSAS<BasicTurbulenceModel>::C2_>
    (
        C2.field(), "C2"
    );
    this->template setDumpField<&mykOmegaSSTSAS<BasicTurbulenceModel>::Lvk_>
    (
        Lvk.field(), "Lvk"
    );

    return fvm::Su
    (
        this->alpha_()*this->rho_()
       *min
        (
            max
            (
                zeta2_*kappa_*S2*sqr(L/Lvk)
              - (2*C_/sigmaPhi_)*this->k_()
               *max
                (
                    magSqr(fvc::grad(this->omega_)()())/sqr(this->omega_()),
                    magSqr(fvc::grad(this->k_)()())/sqr(this->k_())
                ),
                dimensionedScalar(dimensionSet(0, 0, -2, 0, 0), Zero)
            ),
            // Limit SAS production of omega for numerical stability,
            // particularly during start-up
            this->omega_()/(0.1*this->omega_.time().deltaT())
        ),
        this->omega_
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
mykOmegaSSTSAS<BasicTurbulenceModel>::mykOmegaSSTSAS
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kOmegaSST<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    Cs_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.11
        )
    ),
    kappa_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    zeta2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "zeta2",
            this->coeffDict_,
            3.51
        )
    ),
    sigmaPhi_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaPhi",
            this->coeffDict_,
            2.0/3.0
        )
    ),
    C_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C",
            this->coeffDict_,
            2
        )
    ),

    delta_
    (
        LESdelta::New
        (
            IOobject::groupName("delta", alphaRhoPhi.group()),
            *this,
            this->coeffDict_
        )
    )
{
    if (type == typeName)
    {
        this->correctNut();
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool mykOmegaSSTSAS<BasicTurbulenceModel>::read()
{
    if (kOmegaSST<BasicTurbulenceModel>::read())
    {
        Cs_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());
        sigmaPhi_.readIfPresent(this->coeffDict());
        zeta2_.readIfPresent(this->coeffDict());
        C_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
