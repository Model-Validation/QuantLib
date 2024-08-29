/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2002, 2003 Ferdinando Ametrano
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006 StatPro Italia srl
 Copyright (C) 2024 SEB AB STh

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

/*! \file b_spline_curve_print_information_hpp
    \brief Example with rate time curves
*/

#ifndef b_spline_curve_print_information_hpp
#define b_spline_curve_print_information_hpp

#include <ql/compounding.hpp>
#include <ql/handle.hpp>
#include <ql/termstructures/yield/piecewiseyieldcurve.hpp>
#include <ql/termstructures/yield/ratehelpers.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/time/date.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/types.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <iomanip>
#include <ios>
#include <iostream>
#include <vector>

namespace RateTimePrintInfo {
    using namespace QuantLib;

    template <typename T, typename U>
    void printInformation(std::vector<ext::shared_ptr<RateHelper>>& rateHelpers,
                          ext::shared_ptr<PiecewiseYieldCurve<T, U>>& piecewiseCurve,
                          Date& settlementDate,
                          Handle<YieldTermStructure>& curve,
                          DayCounter& dayCounter) {
        std::cout << "Rate helper data: " << std::endl << std::endl;
        for (const auto& helper : rateHelpers) {
            std::cout << io::iso_date(helper->pillarDate()) << ":\t" << std::setprecision(3)
                      << io::rate(helper->quote()->value()) << std::endl;
        }

        // Print curve data
        const auto& times = piecewiseCurve->times();
        const auto& dates = piecewiseCurve->dates();
        const auto& data = piecewiseCurve->data();

        std::cout << std::endl << "Curve raw rates and raw rate / time:" << std::endl << std::endl;
        for (size_t i = 0; i < data.size(); ++i) {
            Rate r;
            if (times[i] == 0.0) {
                r = data[i + 1] / times[i + 1];
            } else {
                r = data[i] / times[i];
            };
            std::cout << io::iso_date(dates[i]) << ":\t" << std::setprecision(3)
                      << io::rate(data[i]) << "\t" << std::setw(10) << io::rate(r) << std::endl;
        }

        // Print out some values from the zero curve
        std::cout << std::endl << "Zero rates from curve:" << std::endl << std::endl << std::endl;
        for (Size i = 0; i < 360; i += 30) {
            Date date = settlementDate + i;
            Rate zeroRate = curve->zeroRate(date, dayCounter, Continuous);
            std::cout << io::iso_date(date) << ":\t" << std::setprecision(3) << io::rate(zeroRate)
                      << std::endl;
        }

        // Print quote errors
        std::cout << std::endl << "Quote errors:" << std::endl << std::endl;
        for (const auto& helper : rateHelpers) {
            try {
                helper->setTermStructure(&*curve.currentLink());
                Date pillarDate = helper->pillarDate();
                Real quoteError = helper->quoteError();
                std::cout << io::iso_date(pillarDate) << ":\t" << std::scientific
                          << std::setprecision(2) << std::showpos << quoteError << std::endl;
            } catch (...) {
                continue;
            }
        }
    }
}

#endif // b_spline_curve_print_information_hpp
