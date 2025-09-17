/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003, 2004, 2005, 2006 Ferdinando Ametrano
 Copyright (C) 2006 StatPro Italia srl

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

/*! \file diffusioncalculator.hpp
    \brief Black/Bachelier-formula calculator class
*/

#ifndef quantlib_diffusioncalculator_hpp
#define quantlib_diffusioncalculator_hpp

#include <ql/instruments/payoffs.hpp>
#include <ql/termstructures/volatility/volatilitytype.hpp>
namespace QuantLib {
    //! Calculator class supporting shifted lognormal (Black) and normal (Bachelier) models
    /*! \bug When the variance is null, division by zero occur during
             the calculation of delta, delta forward, gamma, gamma
             forward, rho, dividend rho, vega, and strike sensitivity.
    */
    class DiffusionCalculator {
      public:
        // Implementation of the calculator shifted lognormal (black) or normal (bachelier)
        class Impl {
          public:
            virtual Real value() const;

            /*! Sensitivity to change in the underlying forward price. */
            virtual Real deltaForward() const;
            /*! Sensitivity to change in the underlying spot price. */
            virtual Real delta(Real spot) const;

            /*! Sensitivity in percent to a percent change in the
                underlying forward price. */
            virtual Real elasticityForward() const;
            /*! Sensitivity in percent to a percent change in the
                underlying spot price. */
            virtual Real elasticity(Real spot) const;

            /*! Second order derivative with respect to change in the
                underlying forward price. */
            virtual Real gammaForward() const;
            /*! Second order derivative with respect to change in the
                underlying spot price. */
            virtual Real gamma(Real spot) const;

            /*! Sensitivity to time to maturity. */
            virtual Real theta(Real spot, Time maturity) const;
            /*! Sensitivity to time to maturity per day,
                assuming 365 day per year. */
            virtual Real thetaPerDay(Real spot, Time maturity) const;

            /*! Sensitivity to volatility. */
            virtual Real vega(Time maturity) const;

            /*! Sensitivity to discounting rate. */
            virtual Real rho(Time maturity) const;

            /*! Sensitivity to dividend/growth rate. */
            virtual Real dividendRho(Time maturity) const;

            /*! Probability of being in the money in the bond martingale
                measure, i.e. N(d2).
                It is a risk-neutral probability, not the real world one.
            */
            virtual Real itmCashProbability() const;

            /*! Probability of being in the money in the asset martingale
                measure, i.e. N(d1).
                It is a risk-neutral probability, not the real world one.
            */
            virtual Real itmAssetProbability() const;

            /*! Sensitivity to strike. */
            virtual Real strikeSensitivity() const;

            /*! gamma w.r.t. strike. */
            virtual Real strikeGamma() const;

            virtual Real alpha() const;
            virtual Real beta() const;
        };


        DiffusionCalculator(const ext::shared_ptr<StrikedTypePayoff>& payoff,
                        Real forward,
                        Real stdDev,
                        Real discount = 1.0,
                        VolatilityType vType = ShiftedLognormal,
                        Real displacement = 0.0);
        DiffusionCalculator(Option::Type optionType,
                        Real strike,
                        Real forward,
                        Real stdDev,
                        Real discount = 1.0,
                        VolatilityType vType = ShiftedLognormal,
                        Real displacement = 0.0);
        virtual ~DiffusionCalculator() = default;

        Real value() const { return impl_->value(); };

        /*! Sensitivity to change in the underlying forward price. */
        Real deltaForward() const { return impl_->deltaForward(); };
        /*! Sensitivity to change in the underlying spot price. */
        virtual Real delta(Real spot) const { return impl_->delta(spot); };

        /*! Sensitivity in percent to a percent change in the
            underlying forward price. */
        Real elasticityForward() const { return impl_->elasticityForward(); };
        /*! Sensitivity in percent to a percent change in the
            underlying spot price. */
        virtual Real elasticity(Real spot) const { return impl_->elasticity(spot); };

        /*! Second order derivative with respect to change in the
            underlying forward price. */
        Real gammaForward() const { return impl_->gammaForward(); };
        /*! Second order derivative with respect to change in the
            underlying spot price. */
        virtual Real gamma(Real spot) const { return impl_->gamma(spot); };

        /*! Sensitivity to time to maturity. */
        virtual Real theta(Real spot, Time maturity) const { return impl_->theta(spot, maturity); };
        /*! Sensitivity to time to maturity per day,
            assuming 365 day per year. */
        virtual Real thetaPerDay(Real spot, Time maturity) const {
            return impl_->thetaPerDay(spot, maturity);
        };

        /*! Sensitivity to volatility. */
        Real vega(Time maturity) const { return impl_->vega(maturity); };

        /*! Sensitivity to discounting rate. */
        Real rho(Time maturity) const { return impl_->rho(maturity); };

        /*! Sensitivity to dividend/growth rate. */
        Real dividendRho(Time maturity) const { return impl_->dividendRho(maturity); };

        /*! Probability of being in the money in the bond martingale
            measure, i.e. N(d2) in Black model .
            It is a risk-neutral probability, not the real world one.
        */
        Real itmCashProbability() const { return impl_->itmCashProbability(); };

        /*! Probability of being in the money in the asset martingale
            measure, i.e. N(d1) in Black model.
            It is a risk-neutral probability, not the real world one.
        */
        Real itmAssetProbability() const { return impl_->itmAssetProbability(); };

        /*! Sensitivity to strike. */
        Real strikeSensitivity() const { return impl_->strikeSensitivity(); };

        /*! gamma w.r.t. strike. */
        Real strikeGamma() const { return impl_->strikeGamma(); };

        Real alpha() const { return impl_->alpha(); };
        Real beta() const const { return impl_->beta(); };

      private:
        std::unique_ptr<Impl> impl_;
    };
}

#endif
