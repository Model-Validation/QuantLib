/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*!
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

#if !defined(BOOST_ALL_NO_LIB) && defined(BOOST_MSVC)
#endif
#include "printinformation.hpp"
#include <ql/compounding.hpp>
#include <ql/currencies/europe.hpp>
#include <ql/handle.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/math/interpolations/linearinterpolation.hpp>
#include <ql/math/interpolations/ratetimeinterpolation.hpp>
#include <ql/quote.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/yield/bootstraptraits.hpp>
#include <ql/termstructures/yield/oisratehelper.hpp>
#include <ql/termstructures/yield/piecewiseyieldcurve.hpp>
#include <ql/termstructures/yield/ratehelpers.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/time/calendar.hpp>
#include <ql/time/calendars/sweden.hpp>
#include <ql/time/date.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/frequency.hpp>
#include <ql/time/period.hpp>
#include <ql/time/timeunit.hpp>
#include <ql/types.hpp>
#include <exception>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

using namespace QuantLib;
using namespace RateTimePrintInfo;

// Function to create rate helpers
std::vector<ext::shared_ptr<RateHelper>>
createRateHelpersOld(const ext::shared_ptr<OvernightIndex>& overnightIndex,
                     const Calendar& calendar) {
    // Define the necessary variables
    std::vector<ext::shared_ptr<RateHelper>> rateHelpers;

    // DepositRateHelper
    Rate depositRate = 0.03746;
    rateHelpers.push_back(ext::make_shared<DepositRateHelper>(depositRate, overnightIndex));

    // DatedOISRateHelper
    std::vector<std::tuple<Rate, Date, Date>> oisData = {
        {0.0375, Date(2, July, 2024), Date(3, July, 2024)},
        {0.0371, Date(20, August, 2024), Date(21, August, 2024)},
        {0.0354, Date(1, October, 2024), Date(2, October, 2024)},
        {0.0339, Date(12, November, 2024), Date(13, November, 2024)},
        {0.0325, Date(7, January, 2025), Date(8, January, 2025)},
        {0.0311, Date(4, February, 2025), Date(5, February, 2025)},
        {0.0298, Date(1, April, 2025), Date(2, April, 2025)},
        {0.0286, Date(13, May, 2025), Date(14, May, 2025)},
        {0.0276, Date(1, July, 2025), Date(2, July, 2025)},
        {0.0266, Date(26, August, 2025), Date(27, August, 2025)},
        {0.0257, Date(30, September, 2025), Date(1, October, 2025)},
        {0.0250, Date(11, November, 2025), Date(12, November, 2025)},
        {0.0244, Date(5, January, 2026), Date(7, January, 2026)},
        {0.0239, Date(10, February, 2026), Date(11, February, 2026)},
        {0.0235, Date(7, April, 2026), Date(8, April, 2026)},
        {0.0232, Date(12, May, 2026), Date(13, May, 2026)}};

    for (const auto& ois : oisData) {
        Rate rate = std::get<0>(ois);
        Date startDate = std::get<1>(ois);
        Date endDate = std::get<2>(ois);
        rateHelpers.push_back(ext::make_shared<DatedOISRateHelper>(
            startDate, endDate, Handle<Quote>(ext::make_shared<SimpleQuote>(rate)),
            overnightIndex));
    }

    return rateHelpers;
}

// Existing createRateHelpers function
static std::vector<ext::shared_ptr<RateHelper>>
createRateHelpers(const ext::shared_ptr<OvernightIndex>& overnightIndex, const Calendar& calendar) {
    // Define the necessary variables
    std::vector<ext::shared_ptr<RateHelper>> rateHelpers;

    // DepositRateHelper
    Rate depositRate = 0.03746;
    rateHelpers.push_back(ext::make_shared<DepositRateHelper>(depositRate, overnightIndex));

    // DatedOISRateHelper
    std::vector<std::tuple<Rate, Date, Date>> oisData = {
        {0.0375, Date(2, July, 2024), Date(3, July, 2024)},
        {0.0371, Date(20, August, 2024), Date(21, August, 2024)},
        {0.0354, Date(1, October, 2024), Date(2, October, 2024)},
        {0.0339, Date(12, November, 2024), Date(13, November, 2024)},
        {0.0325, Date(7, January, 2025), Date(8, January, 2025)},
        {0.0311, Date(4, February, 2025), Date(5, February, 2025)},
        {0.0298, Date(1, April, 2025), Date(2, April, 2025)},
        {0.0286, Date(13, May, 2025), Date(14, May, 2025)},
        {0.0276, Date(1, July, 2025), Date(2, July, 2025)},
        {0.0266, Date(26, August, 2025), Date(27, August, 2025)},
        {0.0257, Date(30, September, 2025), Date(1, October, 2025)},
        {0.0250, Date(11, November, 2025), Date(12, November, 2025)},
        {0.0244, Date(5, January, 2026), Date(7, January, 2026)},
        {0.0239, Date(10, February, 2026), Date(11, February, 2026)},
        {0.0235, Date(7, April, 2026), Date(8, April, 2026)},
        {0.0232, Date(12, May, 2026), Date(13, May, 2026)}};

    for (const auto& ois : oisData) {
        Rate rate = std::get<0>(ois);
        Date startDate = std::get<1>(ois);
        Date endDate = std::get<2>(ois);
        rateHelpers.push_back(ext::make_shared<DatedOISRateHelper>(
            startDate, endDate, Handle<Quote>(ext::make_shared<SimpleQuote>(rate)),
            overnightIndex));
    }

    // Additional OISRateHelpers
    std::vector<std::pair<Period, Rate>> oisRateData = {
        {2 * Years, 0.02870671},  {3 * Years, 0.02661334},  {4 * Years, 0.02546776},
        {5 * Years, 0.02477894},  {6 * Years, 0.02451402},  {7 * Years, 0.02444801},
        {8 * Years, 0.02449672},  {9 * Years, 0.02459572},  {10 * Years, 0.02467563},
        {12 * Years, 0.02483502}, {15 * Years, 0.02489491}, {20 * Years, 0.02429415},
        {25 * Years, 0.02333397}, {30 * Years, 0.02235453}};

    for (const auto& oisRate : oisRateData) {
        Period tenor = oisRate.first;
        Rate rate = oisRate.second;
        rateHelpers.push_back(ext::make_shared<OISRateHelper>(
            2, tenor, Handle<Quote>(ext::make_shared<SimpleQuote>(rate)), overnightIndex));
    }

    return rateHelpers;
}


int main() {
    /* ### Setup phase ### */
    Calendar calendar = Sweden();

    // Define the SEK-STINATN index
    ext::shared_ptr<OvernightIndex> overnightIndex =
        ext::make_shared<OvernightIndex>("SEK-STINATN", 1, SEKCurrency(), calendar, Actual360());

    // Create rate helpers
    std::vector<ext::shared_ptr<RateHelper>> rateHelpers =
        createRateHelpers(overnightIndex, calendar);

    // Slice the vector to include only the short end
    std::vector<ext::shared_ptr<RateHelper>> shortEndHelpers(rateHelpers.begin(),
                                                             rateHelpers.begin() + 16);

    // Bootstrap the piecewise linear zero curve
    DayCounter dayCounter = Actual360();

    Date settlementDate = calendar.advance(Date::todaysDate(), 2, Days);

    /* ### Test PiecewiseYieldCurve<RateTime, Linear> implementation ### */
    std::cout << std::endl
              << "### Test PiecewiseYieldCurve<RateTime, Linear> implementation ###" << std::endl
              << std::endl;

    try {
        ext::shared_ptr<PiecewiseYieldCurve<RateTime, Linear>> piecewiseLinearRateTimeCurve =
            ext::make_shared<PiecewiseYieldCurve<RateTime, Linear>>(settlementDate, shortEndHelpers,
                                                                    dayCounter);

        Handle<YieldTermStructure> curveLinearRateTime(piecewiseLinearRateTimeCurve);

        // Initiate lazy bootstrap
        curveLinearRateTime->zeroRate(settlementDate, dayCounter, Continuous, NoFrequency);

        // Print the rate helper quotes for verification
        printInformation(shortEndHelpers, piecewiseLinearRateTimeCurve,
                                          settlementDate,
                         curveLinearRateTime, dayCounter);


    } catch (std::exception& e) {
        std::cerr << "Piecewise RateTime, Linear curve error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown piecewise RateTime, Linear zero curve error" << std::endl;
        return 1;
    }

    /* ### Test PiecewiseYieldCurve<ZeroYield, RateTimeLinear> implementation ### */
    std::cout << std::endl
              << "### Test PiecewiseYieldCurve<ZeroYield, RateTimeLinear> implementation ###" << std::endl
              << std::endl;

    try {
        ext::shared_ptr<PiecewiseYieldCurve<ZeroYield, RateTimeLinear>>
            piecewiseRateTimeLinearCurve =
                ext::make_shared<PiecewiseYieldCurve<ZeroYield, RateTimeLinear>>(
                    settlementDate, shortEndHelpers, dayCounter);

        Handle<YieldTermStructure> curveRateTimeLinear(piecewiseRateTimeLinearCurve);

        // Initiate lazy bootstrap
        curveRateTimeLinear->zeroRate(settlementDate, dayCounter, Continuous, NoFrequency);

        // Print the rate helper quotes for verification
        printInformation(shortEndHelpers, piecewiseRateTimeLinearCurve, settlementDate,
                         curveRateTimeLinear, dayCounter);

    } catch (std::exception& e) {
        std::cerr << "Piecewise ZeroYield, RateTimeLinear curve error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown ZeroYield, RateTimeLinear curve error" << std::endl;
        return 1;
    }

    /* ### Test PiecewiseYieldCurve<ZeroYield, MixedRateTimeLinearParabolic> implementation ### */
    std::cout << std::endl
              << "### Test PiecewiseYieldCurve<ZeroYield, MixedRateTimeLinearParabolic> implementation ###"
              << std::endl
              << std::endl;

    try {
        ext::shared_ptr<PiecewiseYieldCurve<ZeroYield, MixedRateTimeLinearParabolic>>
            piecewiseRateTimeLinearCurve =
                ext::make_shared<PiecewiseYieldCurve<ZeroYield, MixedRateTimeLinearParabolic>>(
                    settlementDate, rateHelpers, dayCounter, MixedRateTimeLinearParabolic(16));

        Handle<YieldTermStructure> curveMixedRateTimeLinearParabolic(piecewiseRateTimeLinearCurve);

        // Initiate lazy bootstrap
        curveMixedRateTimeLinearParabolic->zeroRate(settlementDate, dayCounter, Continuous,
                                                    NoFrequency);

        // Print the rate helper quotes for verification
        printInformation(rateHelpers, piecewiseRateTimeLinearCurve, settlementDate,
                         curveMixedRateTimeLinearParabolic, dayCounter);

    } catch (std::exception& e) {
        std::cerr << "Piecewise ZeroYield, MixedRateTimeLinearParabolic curve error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown piecewise ZeroYield, MixedRateTimeLinearParabolic curve error" << std::endl;
        return 1;
    }
    return 0;
}
