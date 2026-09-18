/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2025 OpenFOAM Foundation
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

#include "coefficientMulticomponentMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::coefficientMulticomponentMixture<ThermoType>::
coefficientMulticomponentMixture
(
    const dictionary& dict
)
:
    multicomponentMixture<ThermoType>(dict),
    mixture_("mixture", this->specieThermos()[0]),
    scratch_("mixture", this->specieThermos()[0])
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const typename
Foam::coefficientMulticomponentMixture<ThermoType>::thermoMixtureType&
Foam::coefficientMulticomponentMixture<ThermoType>::thermoMixture
(
    const scalarFieldListSlice& Y
) const
{
    const boolList& active(this->speciesActive());

    // Scale in place rather than forming Y[i]*thermo temporaries. Each such
    // temporary copy-constructs the specie, name string included, three
    // times on its way through the equation of state, the thermo and the
    // transport wrapper; assignment does not copy the name (specie::operator=)
    // and operator*= scales only Y. s*Y and Y*s are the same IEEE product,
    // so the mixture is bit-identical. Profiled at 6 % of a production step
    // on a 33-species case (rde_engine benchmarks/openfoam-hotpath).
    mixture_ = this->specieThermos()[0];
    mixture_ *= Y[0];

    for (label i=1; i<Y.size(); i++)
    {
        if (active[i])
        {
            scratch_ = this->specieThermos()[i];
            scratch_ *= Y[i];
            mixture_ += scratch_;
        }
    }

    return mixture_;
}


template<class ThermoType>
const typename
Foam::coefficientMulticomponentMixture<ThermoType>::transportMixtureType&
Foam::coefficientMulticomponentMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice& Y
) const
{
    return thermoMixture(Y);
}


template<class ThermoType>
const typename
Foam::coefficientMulticomponentMixture<ThermoType>::transportMixtureType&
Foam::coefficientMulticomponentMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice&,
    const thermoMixtureType& mixture
) const
{
    return mixture;
}


// ************************************************************************* //
