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
#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
namespace QuantLib {

    // TODO move into own hpp/cpp
    enum class DiffusionModelType {
      Black,
      Bachelier,
      AsInputVolatilityType
    };

    // Return the implied volatility and its type (shifted lognormal or normal) and a displacement
    // given the model type requested (Black, Bachelier, or AsInputVolatilityType)
    // and the input vol termstructure with a given volType and (shifted lognormal or normal) and
    // displacement
    std::tuple<double, VolatilityType, double>
    convertInputVariance(DiffusionModelType outputModelType,
                                      double displacement,
                                      QuantLib::ext::shared_ptr<BlackVolTermStructure> volTS,
                                      double forward,
                                      double strike,
                                      double t);

    std::tuple<double, VolatilityType, double>
    convertInputVolatility(DiffusionModelType outputModelType,
                                 double displacement,
                                 QuantLib::ext::shared_ptr<BlackVolTermStructure> volTS,
                                 double forward,
                                 double strike,
                                 double t);

    double convertNormalToShiftedLogNormalVol(
        double forward, double strike, double ttm, double nVol, double displacement);

    double convertShiftedLognormalToNormalVol(
        double forward, double strike, double ttm, double slnVol, double displacement);


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

            virtual ~Impl() = default;
            virtual Real value() const = 0;

            /*! Sensitivity to change in the underlying forward price. */
            virtual Real deltaForward() const = 0;
            /*! Sensitivity to change in the underlying spot price. */
            virtual Real delta(Real spot) const = 0;

            /*! Sensitivity in percent to a percent change in the
                underlying forward price. */
            virtual Real elasticityForward() const = 0;
            /*! Sensitivity in percent to a percent change in the
                underlying spot price. */
            virtual Real elasticity(Real spot) const = 0;

            /*! Second order derivative with respect to change in the
                underlying forward price. */
            virtual Real gammaForward() const = 0;
            /*! Second order derivative with respect to change in the
                underlying spot price. */
            virtual Real gamma(Real spot) const = 0;

            /*! Sensitivity to time to maturity. */
            virtual Real theta(Real spot, Time maturity) const = 0;
            /*! Sensitivity to time to maturity per day,
                assuming 365 day per year. */
            virtual Real thetaPerDay(Real spot, Time maturity) const = 0;

            /*! Sensitivity to volatility. */
            virtual Real vega(Time maturity) const = 0;

            /*! Sensitivity to discounting rate. */
            virtual Real rho(Time maturity) const = 0;

            /*! Sensitivity to dividend/growth rate. */
            virtual Real dividendRho(Time maturity) const = 0;

            /*! Probability of being in the money in the bond martingale
                measure, i.e. N(d2).
                It is a risk-neutral probability, not the real world one.
            */
            virtual Real itmCashProbability() const = 0;

            /*! Probability of being in the money in the asset martingale
                measure, i.e. N(d1).
                It is a risk-neutral probability, not the real world one.
            */
            virtual Real itmAssetProbability() const = 0;

            /*! Sensitivity to strike. */
            virtual Real strikeSensitivity() const = 0;

            /*! gamma w.r.t. strike. */
            virtual Real strikeGamma() const = 0;

            virtual Real alpha() const = 0;
            virtual Real beta() const = 0;
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
        Real beta() const { return impl_->beta(); };

      private:
        std::unique_ptr<Impl> impl_;
    };
}

#endif
