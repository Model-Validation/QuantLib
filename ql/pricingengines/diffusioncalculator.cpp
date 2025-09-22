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

#include <ql/math/comparison.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/pricingengines/bacheliercalculator.hpp>
#include <ql/pricingengines/blackcalculator.hpp>
#include <ql/pricingengines/blackformula.hpp>
#include <ql/pricingengines/diffusioncalculator.hpp>
namespace QuantLib {

    double convertEuropeanImpliedNormalVolToShiftedLogNormalVol(
        double forward, double strike, double ttm, double nVol, double displacement) {
        auto optionType = forward > strike ? Option::Type::Put : Option::Type::Call;
        double premium = bachelierBlackFormula(optionType, strike, forward, nVol * ttm);
        double slnVol =
            blackFormulaImpliedStdDev(optionType, strike, forward, premium, 1.0, displacement);
        return slnVol / sqrt(ttm);
    }

    double convertEuropeanImpliedShiftedLognormalVolToNormalVol(
        double forward, double strike, double ttm, double slnVol, double displacement) {
        auto optionType = forward > strike ? Option::Type::Put : Option::Type::Call;
        double price = blackFormula(optionType, strike, forward, slnVol, 1.0, displacement);
        double nVol = bachelierBlackFormulaImpliedVol(optionType, strike, forward, ttm, price);
        return nVol;
    }


        // Return the implied volatility and its type (shifted lognormal or normal) and a displacement
    // given the model type requested (Black, Bachelier, or AsInputVolatilityType)
    // and the input vol termstructure with a given volType and (shifted lognormal or normal) and
    // displacement
    std::tuple<double, VolatilityType, double>
    getImpliedVarianceFromModelType(DiffusionModelType outputModelType,
                                      double displacement,
                                      QuantLib::ext::shared_ptr<BlackVolTermStructure> volTS,
                                      double forward,
                                      double strike,
                                      double t) {

        Real variance = volTS->blackVariance(t, strike);
        VolatilityType volType = volTS->volType();
        Real inputDisplacement = volTS->shift();
        Real outputDisplacement = displacement;
        if (outputModelType == DiffusionModelType::Bachelier &&
            volType == VolatilityType::ShiftedLognormal) {
            double slnVol = volTS->blackVol(t, strike);

            double nVol = convertEuropeanImpliedShiftedLognormalVolToNormalVol(
                forward, strike, t, slnVol, inputDisplacement);
            volType = VolatilityType::Normal;
            outputDisplacement = 0;
            variance = nVol * nVol * t;
        } else if (outputModelType == DiffusionModelType::Black &&
                   volType == VolatilityType::Normal) {

            double nVol = volTS->blackVol(t, strike);
            double slnVol = convertEuropeanImpliedNormalVolToShiftedLogNormalVol(
                forward, strike, t, nVol, outputDisplacement);
            volType = VolatilityType::ShiftedLognormal;
            variance = slnVol * slnVol * t;
        }
        return std::make_tuple(variance, volType, outputDisplacement);
    }

    DiffusionCalculator::DiffusionCalculator(const ext::shared_ptr<StrikedTypePayoff>& p,
                                     Real forward,
                                     Real stdDev,
                                     Real discount,
                                     VolatilityType vType,
                                     Real displacement)
    {
        switch(vType){
            case ShiftedLognormal:
                impl_ = std::make_unique<BlackCalculator>(p, forward, stdDev, discount, displacement);
                break;
            case Normal:
                impl_ = std::make_unique<BachelierCalculator>(p, forward, stdDev, discount);
                break;
            default:
                QL_FAIL("unknown volatility type");
        }
    }

    DiffusionCalculator::DiffusionCalculator(Option::Type optionType,
                                     Real strike,
                                     Real forward,
                                     Real stdDev,
                                     Real discount, 
                                     VolatilityType vType,
                                     Real displacement)
    {
        switch(vType){
            case ShiftedLognormal:
                impl_ = std::make_unique<BlackCalculator>(optionType, strike, forward, stdDev, discount, displacement);
                break;
            case Normal:
                impl_ = std::make_unique<BachelierCalculator>(optionType, strike, forward, stdDev, discount);
                break;
            default:
                QL_FAIL("unknown volatility type");
        }
    }
}
