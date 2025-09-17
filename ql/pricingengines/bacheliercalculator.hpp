/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2025 AcadiaSoft Inc.

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

/*! \file bacheliercalculator.hpp
    \brief Bachelier-formula calculator class
*/

#ifndef quantlib_bacheliercalculator_hpp
#define quantlib_bacheliercalculator_hpp

#include <ql/instruments/payoffs.hpp>

namespace QuantLib {

    //! Bachelier 1976 calculator class
    /*! \bug When the variance is null, division by zero occur during
             the calculation of delta, delta forward, gamma, gamma
             forward, rho, dividend rho, vega, and strike sensitivity.
    */
    class BachelierCalculator {
      private:
        class Calculator;
      public:
        BachelierCalculator(const ext::shared_ptr<StrikedTypePayoff>& payoff,
                        Real forward,
                        Real stdDev,
                        Real discount = 1.0);
        BachelierCalculator(Option::Type optionType,
                        Real strike,
                        Real forward,
                        Real stdDev,
                        Real discount = 1.0);
        virtual ~BachelierCalculator() = default;

        Real value() const;

        /*! Sensitivity to change in the underlying forward price. */
        Real deltaForward() const;
        /*! Sensitivity to change in the underlying spot price. */
        virtual Real delta(Real spot) const;

        /*! Sensitivity in percent to a percent change in the
            underlying forward price. */
        Real elasticityForward() const;
        /*! Sensitivity in percent to a percent change in the
            underlying spot price. */
        virtual Real elasticity(Real spot) const;

        /*! Second order derivative with respect to change in the
            underlying forward price. */
        Real gammaForward() const;
        /*! Second order derivative with respect to change in the
            underlying spot price. */
        virtual Real gamma(Real spot) const;

        /*! Sensitivity to time to maturity. */
        virtual Real theta(Real spot,
                           Time maturity) const;
        /*! Sensitivity to time to maturity per day,
            assuming 365 day per year. */
        virtual Real thetaPerDay(Real spot,
                                 Time maturity) const;

        /*! Sensitivity to volatility. */
        Real vega(Time maturity) const;

        /*! Sensitivity to discounting rate. */
        Real rho(Time maturity) const;

        /*! Sensitivity to dividend/growth rate. */
        Real dividendRho(Time maturity) const;

        /*! Probability of being in the money in the bond martingale
            measure, i.e. N(d2).
            It is a risk-neutral probability, not the real world one.
        */
        Real itmCashProbability() const;

        /*! Probability of being in the money in the asset martingale
            measure, i.e. N(d1).
            It is a risk-neutral probability, not the real world one.
        */
        Real itmAssetProbability() const;

        /*! Sensitivity to strike. */
        Real strikeSensitivity() const;

        /*! gamma w.r.t. strike. */
        Real strikeGamma() const;

        Real alpha() const;
        Real beta() const;
      protected:
        void initialize(const ext::shared_ptr<StrikedTypePayoff>& p);
        Real strike_, forward_, stdDev_, discount_, variance_;
        Real d_;
        Real alpha_, beta_, DalphaDd_, D2alphaD2d_, DbetaDd_, D2betaD2d_;
        Real n_d_, cum_d_;
        Real y_, DyDs_, DyDstrike_, DyDforward_;
        Real x_, DxDs_, DxDstrike_, DxDforward_;
    };

    // inline
    inline Real BachelierCalculator::thetaPerDay(Real spot,
                                             Time maturity) const {
        return theta(spot, maturity)/365.0;
    }

    inline Real BachelierCalculator::itmCashProbability() const {
        return cum_d_;
    }

    inline Real BachelierCalculator::itmAssetProbability() const {
        return cum_d_;
    }

    inline Real BachelierCalculator::alpha() const {
        return alpha_;
    }

    inline Real BachelierCalculator::beta() const {
        return beta_;
    }

}

#endif
