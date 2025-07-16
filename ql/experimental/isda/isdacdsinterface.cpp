/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2024 QuantLib Contributors

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.
*/

#define QL_ENABLE_ISDA_CDS
#ifdef QL_ENABLE_ISDA_CDS

#    include <ql/errors.hpp>
#    include <ql/experimental/isda/isdacdsinterface.hpp>
#    include <ql/math/interpolations/linearinterpolation.hpp>
#    include <ql/time/calendars/weekendsonly.hpp>
#    include <ql/time/daycounters/actual360.hpp>
#    include <ql/time/daycounters/actual365fixed.hpp>

// Include ISDA CDS Standard Model headers
extern "C" {
#    include "isda/bastypes.h"
#    include "isda/cdate.h"
#    include "isda/cds.h"
#    include "isda/convert.h"
#    include "isda/ldate.h"
#    include "isda/macros.h"
#    include "isda/stub.h"
#    include "isda/tcurve.h"
}

namespace QuantLib {


    TCurve* IsdaCdsInterface::convertYieldCurve(const Handle<YieldTermStructure>& yieldCurve,
                                                const Date& baseDate) {

        QL_REQUIRE(!yieldCurve.empty(), "Yield curve handle is empty");

        // Build curve dates (quarterly out to 30 years)
        std::vector<Date> dates = buildCurveDates(yieldCurve, baseDate);
        std::vector<Real> discountFactors = extractDiscountFactors(yieldCurve, dates);

        // Convert to ISDA format
        std::vector<long> isdaDates;
        std::vector<double> isdaRates;

        for (Size i = 0; i < dates.size(); ++i) {
            isdaDates.push_back(dateToTDate(dates[i]));

            // Convert discount factor to zero rate
            Real yearFraction = Actual365Fixed().yearFraction(baseDate, dates[i]);
            if (yearFraction > 0.0) {
                Real zeroRate = -std::log(discountFactors[i]) / yearFraction;
                isdaRates.push_back(zeroRate);
            } else {
                isdaRates.push_back(0.0);
            }
        }

        // Create ISDA curve
        TCurve* curve = JpmcdsMakeTCurve(dateToTDate(baseDate), isdaDates.data(), isdaRates.data(),
                                         static_cast<int>(isdaDates.size()),
                                         1.0,            // Annual compounding
                                         JPMCDS_ACT_365F // Actual/365 Fixed day count
        );

        QL_REQUIRE(curve != nullptr, "Failed to create ISDA curve");
        return curve;
    }

    TCurve* IsdaCdsInterface::convertDefaultCurve(
        const Handle<DefaultProbabilityTermStructure>& defaultCurve,
        const Date& baseDate,
        Real recoveryRate) {

        QL_REQUIRE(!defaultCurve.empty(), "Default curve handle is empty");

        // Build curve dates
        std::vector<Date> dates = buildCurveDates(Handle<YieldTermStructure>(), baseDate);
        std::vector<Real> hazardRates = extractHazardRates(defaultCurve, dates);

        // Convert to ISDA format
        std::vector<long> isdaDates;
        std::vector<double> isdaRates;

        for (Size i = 0; i < dates.size(); ++i) {
            isdaDates.push_back(dateToTDate(dates[i]));
            isdaRates.push_back(hazardRates[i]);
        }

        // Create ISDA curve
        TCurve* curve = JpmcdsMakeTCurve(dateToTDate(baseDate), isdaDates.data(), isdaRates.data(),
                                         static_cast<int>(isdaDates.size()),
                                         1.0,            // Annual compounding
                                         JPMCDS_ACT_365F // Actual/365 Fixed day count
        );

        QL_REQUIRE(curve != nullptr, "Failed to create ISDA spread curve");
        return curve;
    }

    Real IsdaCdsInterface::priceCds(const CreditDefaultSwap& cds,
                                    const Handle<YieldTermStructure>& yieldCurve,
                                    const Handle<DefaultProbabilityTermStructure>& defaultCurve,
                                    Real recoveryRate,
                                    const Date& today) {

        // Convert curves to ISDA format
        std::unique_ptr<TCurve, decltype(&freeTCurve)> discCurve(
            convertYieldCurve(yieldCurve, today), &freeTCurve);
        std::unique_ptr<TCurve, decltype(&freeTCurve)> spreadCurve(
            convertDefaultCurve(defaultCurve, today, recoveryRate), &freeTCurve);

        // Extract CDS parameters
        Date startDate = cds.protectionStartDate();
        Date endDate = cds.maturity();
        Real notional = cds.notional();
        Real spread = cds.runningSpread();

        // Calculate step-in date (typically T+1)
        const Date stepInDate = WeekendsOnly().advance(today, 1, Days);

        // Price the CDS using ISDA methodology
        double price = 0.0;
        int result = JpmcdsCdsPrice(dateToTDate(today),      // today
                                    dateToTDate(today),      // valueDate
                                    dateToTDate(stepInDate), // stepinDate
                                    dateToTDate(startDate),  // startDate
                                    dateToTDate(endDate),    // endDate
                                    spread,                  // couponRate
                                    TRUE,                    // payAccOnDefault
                                    nullptr,                 // couponInterval (use default 3M)
                                    nullptr,                 // stubType (use default)
                                    JPMCDS_ACT_360,          // paymentDcc
                                    JPMCDS_BAD_DAY_MODIFIED, // badDayConv
                                    nullptr,                 // calendar
                                    discCurve.get(),         // discCurve
                                    spreadCurve.get(),       // spreadCurve
                                    recoveryRate,            // recoveryRate
                                    TRUE,                    // isPriceClean
                                    &price                   // output price
        );

        QL_REQUIRE(result == SUCCESS, "ISDA CDS pricing failed");

        return price * notional;
    }

    Real IsdaCdsInterface::calculateParSpread(
        const Schedule& schedule,
        const Handle<YieldTermStructure>& yieldCurve,
        const Handle<DefaultProbabilityTermStructure>& defaultCurve,
        Real recoveryRate,
        const Date& today,
        const Date& protectionStart) {

        // Convert curves to ISDA format
        std::unique_ptr<TCurve, decltype(&freeTCurve)> discCurve(
            convertYieldCurve(yieldCurve, today), &freeTCurve);
        std::unique_ptr<TCurve, decltype(&freeTCurve)> spreadCurve(
            convertDefaultCurve(defaultCurve, today, recoveryRate), &freeTCurve);

        // Calculate step-in date (typically T+1)
        Date stepinDate = WeekendsOnly().advance(today, 1, Days);

        // Single maturity for par spread calculation
        Date endDate = schedule.endDate();
        double parSpread = 0.0;

        long endTDate = dateToTDate(endDate);
        int result = JpmcdsCdsParSpreads(dateToTDate(today),           // today
                                         dateToTDate(stepinDate),      // stepinDate
                                         dateToTDate(protectionStart), // startDate
                                         1,                            // nbEndDates
                                         &endTDate,                    // endDates
                                         TRUE,                         // payAccOnDefault
                                         nullptr,                 // couponInterval (use default 3M)
                                         nullptr,                 // stubType (use default)
                                         JPMCDS_ACT_360,          // paymentDcc
                                         JPMCDS_BAD_DAY_MODIFIED, // badDayConv
                                         nullptr,                 // calendar
                                         discCurve.get(),         // discCurve
                                         spreadCurve.get(),       // spreadCurve
                                         recoveryRate,            // recoveryRate
                                         &parSpread               // output par spread
        );

        QL_REQUIRE(result == SUCCESS, "ISDA par spread calculation failed");

        return parSpread;
    }

    TCurve* IsdaCdsInterface::bootstrapDefaultCurve(const Handle<YieldTermStructure>& yieldCurve,
                                                    const std::vector<Real>& cdsQuotes,
                                                    const std::vector<Period>& tenors,
                                                    Real recoveryRate,
                                                    const Date& today) {

        QL_REQUIRE(cdsQuotes.size() == tenors.size(), "CDS quotes and tenors must have same size");

        // Convert curves to ISDA format
        std::unique_ptr<TCurve, decltype(&freeTCurve)> discCurve(
            convertYieldCurve(yieldCurve, today), &freeTCurve);

        // Build end dates from tenors
        std::vector<Date> endDates;
        std::vector<long> isdaEndDates;
        std::vector<double> isdaCouponRates;

        for (Size i = 0; i < tenors.size(); ++i) {
            Date endDate = WeekendsOnly().advance(today, tenors[i]);
            endDates.push_back(endDate);
            isdaEndDates.push_back(dateToTDate(endDate));
            isdaCouponRates.push_back(cdsQuotes[i]);
        }

        // Calculate step-in date and protection start
        Date stepinDate = WeekendsOnly().advance(today, 1, Days);
        Date protectionStart = today;

        // Bootstrap the curve
        TCurve* curve = JpmcdsCleanSpreadCurve(dateToTDate(today),               // today
                                               discCurve.get(),                  // discCurve
                                               dateToTDate(protectionStart),     // startDate
                                               dateToTDate(stepinDate),          // stepinDate
                                               dateToTDate(today),               // cashSettleDate
                                               static_cast<long>(tenors.size()), // nbDate
                                               isdaEndDates.data(),              // endDates
                                               isdaCouponRates.data(),           // couponRates
                                               nullptr,        // includes (all included)
                                               recoveryRate,   // recoveryRate
                                               TRUE,           // payAccOnDefault
                                               nullptr,        // couponInterval (use default 3M)
                                               JPMCDS_ACT_360, // paymentDcc
                                               nullptr,        // stubType (use default)
                                               JPMCDS_BAD_DAY_MODIFIED, // badDayConv
                                               nullptr                  // calendar
        );

        QL_REQUIRE(curve != nullptr, "Failed to bootstrap ISDA spread curve");
        return curve;
    }

    void IsdaCdsInterface::freeTCurve(TCurve* curve) {
        if (curve) {
            JpmcdsFreeTCurve(curve);
        }
    }

    long IsdaCdsInterface::dateToTDate(const Date& date) {
        return static_cast<long>(date.serialNumber() + 25569); // Excel to ISDA date conversion
    }

    Date IsdaCdsInterface::tDateToDate(const long tDate) {
        return Date(static_cast<Date::serial_type>(tDate - 25569));
    }

    std::vector<Date> IsdaCdsInterface::buildCurveDates(const Handle<YieldTermStructure>& curve,
                                                        const Date& baseDate) {

        std::vector<Date> dates;

        // Build standard curve grid: 1W, 2W, 1M, 2M, 3M, 6M, 1Y, 2Y, ..., 30Y
        WeekendsOnly calendar;

        dates.push_back(calendar.advance(baseDate, 1, Weeks));
        dates.push_back(calendar.advance(baseDate, 2, Weeks));
        dates.push_back(calendar.advance(baseDate, 1, Months));
        dates.push_back(calendar.advance(baseDate, 2, Months));
        dates.push_back(calendar.advance(baseDate, 3, Months));
        dates.push_back(calendar.advance(baseDate, 6, Months));

        // Annual points
        for (int i = 1; i <= 30; ++i) {
            dates.push_back(calendar.advance(baseDate, i, Years));
        }

        return dates;
    }

    std::vector<Real>
    IsdaCdsInterface::extractDiscountFactors(const Handle<YieldTermStructure>& curve,
                                             const std::vector<Date>& dates) {

        std::vector<Real> discountFactors;
        discountFactors.reserve(dates.size());

        std::transform(dates.begin(), dates.end(), std::back_inserter(discountFactors),
                       [&curve](const Date& date) { return curve->discount(date); });

        return discountFactors;
    }

    std::vector<Real>
    IsdaCdsInterface::extractHazardRates(const Handle<DefaultProbabilityTermStructure>& curve,
                                         const std::vector<Date>& dates) {

        std::vector<Real> hazardRates;
        hazardRates.reserve(dates.size());

        std::transform(dates.begin(), dates.end(), std::back_inserter(hazardRates),
                       [&curve](const Date& date) { return curve->hazardRate(date); });

        return hazardRates;
    }

}

#endif // QL_ENABLE_ISDA_CDS