/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2022 Skandinaviska Enskilda Bankan AB (publ)

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file svivolsurface.hpp
    \brief SVI volatility (smile) surface
*/

#ifndef quantlib_svi_vol_surface_hpp
#define quantlib_svi_vol_surface_hpp

#include <ql/experimental/volatility/sviinterpolatedsmilesection.hpp>
#include <ql/experimental/volatility/svismilesection.hpp>
#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>

namespace QuantLib {


    //! SVI volatility (smile) surface
    /*! blah blah

        \test None yet
     */
    class SviVolSurface : public BlackVarianceTermStructure {
      public:
        //! \name Constructors
        //@{
        /*! @param todo TO DO
        */
        SviVolSurface(const std::vector<SviSmileSection>& smileSections);
        //@}

        // Inherited via BlackVarianceTermStructure
        //! Inspectors
        //@{
        Date maxDate() const override;
        Real minStrike() const override;
        Real maxStrike() const override;
        //@}

        Real blackVarianceImpl(Time t, Real strike) const override;

      private:
        
    };


    inline Date SviVolSurface::maxDate() const { return Date(); }

    inline Real SviVolSurface::minStrike() const { return 0.0; }

    inline Real SviVolSurface::maxStrike() const { return QL_MAX_REAL; }
}
#endif
