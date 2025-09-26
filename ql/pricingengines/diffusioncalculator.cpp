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

    double convertNormalToShiftedLogNormalVol(
        double forward, double strike, double ttm, double nVol, double displacement) {
        QL_REQUIRE(!close_enough(forward, 0.0), "forward is 0.0, can not convert normal vol to shifted lognormal vol");
        auto optionType = forward > strike ? Option::Type::Put : Option::Type::Call;
        double premium = bachelierBlackFormula(optionType, strike, forward, nVol * std::sqrt(ttm));
        double slnVol =
            blackFormulaImpliedStdDev(optionType, strike, forward, premium, 1.0, displacement, nVol / forward );
        return slnVol / sqrt(ttm);
    }

    double convertShiftedLognormalToNormalVol(
        double forward, double strike, double ttm, double slnVol, double displacement) {
        auto optionType = forward > strike ? Option::Type::Put : Option::Type::Call;
        double price = blackFormula(optionType, strike, forward, slnVol * std::sqrt(ttm), 1.0, displacement);
        double nVol = bachelierBlackFormulaImpliedVol(optionType, strike, forward, ttm, price);
        return nVol;
    }

    double convertShiftedLognormalToShiftedLognormalVol(
        double forward, double strike, double ttm, double slnVol, double oldDisplacement,
        double newDisplacement) {
        if (close_enough(oldDisplacement, newDisplacement))
            return slnVol;
        auto optionType = forward > strike ? Option::Type::Put : Option::Type::Call;
        double price = blackFormula(optionType, strike, forward, slnVol, 1.0, oldDisplacement);
        double newSlnVol =
            blackFormulaImpliedStdDev(optionType, strike, forward, price, 1.0, newDisplacement);
        return newSlnVol / sqrt(ttm);
    }

    std::tuple<double, VolatilityType, double>
    convertInputVolatility(DiffusionModelType outputModelType,
                           double outputDisplacement,
                           VolatilityType inputVolType,
                           double inputDisplacement,
                           double inputVol,
                           double forward,
                           double strike,
                           double t) {
        if (outputModelType == DiffusionModelType::AsInputVolatilityType){
            return std::make_tuple(inputVol, inputVolType, inputDisplacement);
        }

        auto targetVolType = outputModelType == DiffusionModelType::Black ?
                                 VolatilityType::ShiftedLognormal :
                                 VolatilityType::Normal;
        auto targetDisplacement =
            outputModelType == DiffusionModelType::Black ? outputDisplacement : 0.0;


        if (inputVolType == VolatilityType::ShiftedLognormal &&
            targetVolType == VolatilityType::Normal) {
            double slnVol = inputVol;
            double nVol =
                convertShiftedLognormalToNormalVol(forward, strike, t, slnVol, inputDisplacement);
            return std::make_tuple(nVol, VolatilityType::Normal, targetDisplacement);
        }

        if (inputVolType == VolatilityType::Normal &&
            targetVolType == VolatilityType::ShiftedLognormal) {
            double nVol = inputVol;
            double slnVol =
                convertNormalToShiftedLogNormalVol(forward, strike, t, nVol, targetDisplacement);
            return std::make_tuple(slnVol, VolatilityType::ShiftedLognormal, targetDisplacement);
        }

        if (inputVolType == VolatilityType::ShiftedLognormal &&
            targetVolType == VolatilityType::ShiftedLognormal &&
            !close_enough(inputDisplacement, targetDisplacement)) {
            // need to convert the vol to the new displacement
            double slnVol = inputVol;
            double newSlnVol = convertShiftedLognormalToShiftedLognormalVol(
                forward, strike, t, slnVol, inputDisplacement, targetDisplacement);
            return std::make_tuple(newSlnVol, VolatilityType::ShiftedLognormal,
                                   targetDisplacement);
        }

        return std::make_tuple(inputVol, targetVolType, targetDisplacement);
    }

    // Convert the given implied volatility to the implied volatility of the target model type
    // (Black, Bachelier, or AsInputVolatilityType i.e no conversion)
    std::tuple<double, VolatilityType, double>
    convertInputVolatility(DiffusionModelType outputModelType,
                                       double displacement,
                                       QuantLib::ext::shared_ptr<BlackVolTermStructure> volTS,
                                       double forward,
                                       double strike,
                                       double t) {
        return convertInputVolatility(outputModelType, displacement, volTS->volType(),
                                      volTS->shift(), volTS->blackVol(t, strike), forward, strike,
                                      t);
    }

    std::tuple<double, VolatilityType, double>
    convertInputVariance(DiffusionModelType outputModelType,
                         double displacement,
                         QuantLib::ext::shared_ptr<BlackVolTermStructure> volTS,
                         double forward,
                         double strike,
                         double t) {

        auto [vol, volType, outputDisplacement] =
            convertInputVolatility(outputModelType, displacement, volTS, forward, strike, t);
        return std::make_tuple(vol * vol * t, volType, outputDisplacement);
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
