/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2025 AcadiaSoft Inc.
 Copyright (C) 2010 Dimitri Reiswich

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

/*! \file bachelierdeltacalculator.hpp
    \brief Bachelier formula delta calculator class, following the design / interface of
   BlackDeltaCalculator in ql/experimental/fx/blackcalculator.hpp
*/

#ifndef quantlib_bachelier_delta_calculator_hpp
#define quantlib_bachelier_delta_calculator_hpp

#include <ql/pricingengines/bacheliercalculator.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/math/solvers1d/brent.hpp>
#include <ql/quotes/deltavolquote.hpp>

namespace QuantLib {

    //! Bachelier delta calculator class
    /*! Class includes many operations needed for converting 
        strikes into deltas and vice versa for the different delta
        types in a normal Bachelier model.
    */
    class BachelierDeltaCalculator {
      public:
        // A parsimonious constructor is chosen, which for example
        // doesn't need a strike. The reason for this is, that we'd
        // like this class to calculate deltas for different strikes
        // many times, e.g. in a numerical routine, which will be the
        // case in the smile setup procedure.
        BachelierDeltaCalculator(Option::Type ot,
                             DeltaVolQuote::DeltaType dt,
                             Real spot,
                             DiscountFactor dDiscount,   // domestic discount
                             DiscountFactor fDiscount,   // foreign discount
                             Real stdDev);

        // Give strike, receive delta according to specified type
        Real deltaFromStrike(Real strike) const;
        // Give delta according to specified type, receive strike
        Real strikeFromDelta(Real delta) const;

        Real moneyness(Real strike) const; // (F-K)/stdDev
        Real cumD(Real strike) const;    // N(d) or N(-d)
        Real nD(Real strike) const;      // n(d)


        void setDeltaType(DeltaVolQuote::DeltaType dt);
        void setOptionType(Option::Type ot);

        // The following function can be calculated without an explicit strike
        Real atmStrike(DeltaVolQuote::AtmType atmT) const;

      private:
        // alternative delta type
        Real strikeFromDelta(Real delta, DeltaVolQuote::DeltaType dt) const;


        DeltaVolQuote::DeltaType dt_;
        Option::Type ot_;
        DiscountFactor dDiscount_, fDiscount_;

        Real stdDev_, spot_, forward_;
        Integer phi_;
        Real fExpPos_,fExpNeg_;
    };


    class BachelierDeltaPremiumAdjustedSolverClass {
      public:
        BachelierDeltaPremiumAdjustedSolverClass(
                        Option::Type ot,
                        DeltaVolQuote::DeltaType dt,
                        Real spot,
                        DiscountFactor dDiscount,   // domestic discount
                        DiscountFactor fDiscount,   // foreign  discount
                        Real stdDev,
                        Real delta);

        Real operator()(Real strike) const;

      private:
        BachelierDeltaCalculator bdc_;
        Real delta_;
    };
  }


#endif
