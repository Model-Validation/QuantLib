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

#include <ql/pricingengines/diffusioncalculator.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/math/comparison.hpp>
#include <ql/pricingengines/bacheliercalculator.hpp>
#include <ql/pricingengines/blackcalculator.hpp>
namespace QuantLib {

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
