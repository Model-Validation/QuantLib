/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*!
 Copyright (C) 2021 Skandinaviska Enskilda Banken AB (publ)

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

#include <ql/qldefines.hpp>
#ifdef BOOST_MSVC
#    include <ql/auto_link.hpp>
#endif

#include <ql/indexes/ibor/all.hpp>
#include <ql/instruments/vanillaswap.hpp>
#include <ql/math/interpolations/linearinterpolation.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>
#include <ql/settings.hpp>
#include <ql/termstructures/yield/forwardcurve.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/date.hpp>
#include <ql/time/daycounters/all.hpp>
#include <ql/cashflows/floatingratecoupon.hpp>

#include <iostream>
#include <iomanip>

using namespace QuantLib;

int main()
{
    std::cout << "Beginning vanilla swap pricer...\n";

    try {
        // First, set the evaluation date and general settings
        Date today(10, June, 2021);
        Settings::instance().evaluationDate() = today;

        // Second, create the market data
        std::vector<Date> fwdDatesSek{
            Date(14, June, 2021), Date(16, June, 2021), Date(15, Sep, 2021), Date(15, Dec, 2021),
            Date(16, Mar, 2022),  Date(15, June, 2022), Date(21, Sep, 2022), Date(21, Dec, 2022),
            Date(15, Mar, 2023),  Date(21, June, 2023), Date(20, Sep, 2023), Date(20, Dec, 2023),
            Date(20, Mar, 2024),  Date(17, Mar, 2025),  Date(3, Nov, 2025),  Date(2, Nov, 2026),
            Date(1, Nov, 2027),   Date(30, Oct, 2028),  Date(29, Oct, 2029), Date(4, Nov, 2030),
            Date(3, May, 2032),   Date(30, Oct, 2034),  Date(1, Nov, 2038),  Date(2, Nov, 2043),
            Date(2, Nov, 2048)};

        std::vector<Rate> fwdRatesSek{-0.038, -0.036, -0.029, -0.096, 0.001,  0.024,  0.051,
                                      0.049,  0.156,  0.209,  0.266,  0.266,  0.399,  0.5783,
                                      0.7554, 0.9498, 1.0863, 1.2242, 1.3294, 1.4232, 1.5225,
                                      1.6955, 1.5774, 1.2658, 1.0346};
        std::transform(fwdRatesSek.begin(), fwdRatesSek.end(), fwdRatesSek.begin(),
                       [](Rate f) { return f / 100.0; });

        Handle<YieldTermStructure> fwdCurve(ext::make_shared<InterpolatedForwardCurve<Linear> >(
            fwdDatesSek, fwdRatesSek, Actual360(), TARGET()));

        Handle<YieldTermStructure> discCurve(
            ext::make_shared<FlatForward>(today, 0.01, Actual360()));

        // Third, create some helper variables
        std::vector<Date> fixedScheduleDates{Date(7, July, 2021), Date(7, July, 2022),
                                             Date(7, July, 2023), Date(8, July, 2024),
                                             Date(7, July, 2025), Date(7, July, 2026)};

        std::vector<Date> floatScheduleDates{
            Date(7, July, 2021), Date(7, Oct, 2021), Date(7, Jan, 2022), Date(7, Apr, 2022),
            Date(7, July, 2022), Date(7, Oct, 2022), Date(9, Jan, 2023), Date(11, Apr, 2023),
            Date(7, July, 2023), Date(9, Oct, 2023), Date(8, Jan, 2024), Date(8, Apr, 2024),
            Date(8, July, 2024), Date(7, Oct, 2024), Date(7, Jan, 2025), Date(7, Apr, 2025),
            Date(7, July, 2025), Date(7, Oct, 2025), Date(7, Jan, 2026), Date(7, Apr, 2026),
            Date(7, July, 2026)};
        
        // Fourth, create the vanilla swap instrument
        Swap::Type type = Swap::Type::Receiver;
        Real nominal = 1e6;
        Schedule fixedSchedule(fixedScheduleDates);
        Rate fixedRate = 0.0183;
        DayCounter fixedDayCount = Thirty360(Thirty360::BondBasis);
        Schedule floatSchedule(floatScheduleDates);
        ext::shared_ptr<IborIndex> iborIndex = ext::make_shared<IborIndex>(SEKLibor(Period(0, Months), fwdCurve));
        Spread spread = 0.0;
        DayCounter floatingDayCount = Actual360();

        VanillaSwap swap(type, nominal, fixedSchedule, fixedRate, fixedDayCount, floatSchedule, iborIndex,
                         spread, floatingDayCount);
        
        // Fifth, create the vanilla swap pricing engine
        swap.setPricingEngine(ext::shared_ptr<PricingEngine>(new DiscountingSwapEngine(discCurve)));

        swap.NPV();

        // Sixth, output the interesting stuff
        std::cout << std::setw(3) << "i, " << std::setw(8) << "date, " << std::setw(10)
                  << "nominal, " << std::setw(15) << "accrualStartDate, " << std::setw(15)
                  << "accrualEndDate, " << std::setw(8) << "rate, " << std::setw(8) << "amount"
                  << std::endl;

        for (size_t i = 0; i < swap.floatingLeg().size(); ++i) {
            std::cout << std::setprecision(10) << std::fixed;
            ext::shared_ptr<FloatingRateCoupon> coupon =
                ext::dynamic_pointer_cast<FloatingRateCoupon>(swap.floatingLeg().at(i));
            std::cout << i << ", " << coupon->rate() << std::endl;
            //swap.floatingLeg().at(i);
        }
        std::cout << "End table.\n" << std::endl;

        std::cout << "Used forward rate (Date(7, July, 2021), Date(7, Oct, 2021)):\n"
                  << fwdCurve->forwardRate(Date(7, July, 2021), Date(7, Oct, 2021),
                                           floatingDayCount, Compounding::Simple) *
                         100
                  << std::endl;
        std::cout << "Expected forward rate (Date(7, July, 2021), Date(7, July 2021)):\n"
                  << fwdCurve->forwardRate(Date(7, July, 2021), Date(7, July, 2021),
                                           floatingDayCount, Compounding::Simple) *
                         100
                  << std::endl;

    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "unknown error" << std::endl;
        return 1;
    }
}

