/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2024 QuantLib Contributors

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.
*/

/*! \file isdacdsinterface.hpp
    \brief Interface to ISDA CDS Standard Model C library
*/

#ifndef quantlib_isda_cds_interface_hpp
#define quantlib_isda_cds_interface_hpp

#ifdef QL_ENABLE_ISDA_CDS

#    include <ql/default.hpp>
#    include <ql/instruments/creditdefaultswap.hpp>
#    include <ql/termstructures/defaulttermstructure.hpp>
#    include <ql/termstructures/yieldtermstructure.hpp>
#    include <ql/time/date.hpp>
#    include <ql/time/schedule.hpp>

// Forward declare ISDA types to avoid including C headers in this header
struct _TCurve; // NOLINT(bugprone-reserved-identifier, clang-diagnostic-reserved-identifier)
typedef struct _TCurve TCurve;

namespace QuantLib {

    //! Interface to ISDA CDS Standard Model C library
    /*! This class provides a bridge between QuantLib objects and the
        ISDA CDS Standard Model C library, allowing for ISDA-compliant
        CDS pricing and analytics.

        \ingroup experimental
    */
    class IsdaCdsInterface {
      public:
        //! Default constructor
        IsdaCdsInterface() = default;

        //! Destructor
        ~IsdaCdsInterface() = default;

        //! Convert QuantLib yield curve to ISDA TCurve format
        /*! \param yieldCurve QuantLib yield term structure
            \param baseDate Base date for the curve
            \return ISDA TCurve pointer (caller must free with freeTCurve)
        */
        static TCurve* convertYieldCurve(const Handle<YieldTermStructure>& yieldCurve,
                                         const Date& baseDate);

        //! Convert QuantLib default probability curve to ISDA spread curve format
        /*! \param defaultCurve QuantLib default probability term structure
            \param baseDate Base date for the curve
            \param recoveryRate Recovery rate assumption
            \return ISDA TCurve pointer (caller must free with freeTCurve)
        */
        static TCurve*
        convertDefaultCurve(const Handle<DefaultProbabilityTermStructure>& defaultCurve,
                            const Date& baseDate,
                            Real recoveryRate = 0.4);

        //! Price a CDS using ISDA standard model
        /*! \param cds QuantLib CDS instrument
            \param yieldCurve Discount curve
            \param defaultCurve Default probability curve
            \param recoveryRate Recovery rate
            \param today Valuation date
            \return CDS price using ISDA methodology
        */
        static Real priceCds(const CreditDefaultSwap& cds,
                             const Handle<YieldTermStructure>& yieldCurve,
                             const Handle<DefaultProbabilityTermStructure>& defaultCurve,
                             Real recoveryRate,
                             const Date& today);

        //! Calculate CDS par spread using ISDA standard model
        /*! \param schedule CDS payment schedule
            \param yieldCurve Discount curve
            \param defaultCurve Default probability curve
            \param recoveryRate Recovery rate
            \param today Valuation date
            \param protectionStart Protection start date
            \return Par spread using ISDA methodology
        */
        static Real calculateParSpread(const Schedule& schedule,
                                       const Handle<YieldTermStructure>& yieldCurve,
                                       const Handle<DefaultProbabilityTermStructure>& defaultCurve,
                                       Real recoveryRate,
                                       const Date& today,
                                       const Date& protectionStart);

        //! Bootstrap default curve from CDS quotes using ISDA methodology
        /*! \param yieldCurve Discount curve
            \param cdsQuotes Vector of CDS spread quotes
            \param tenors Vector of CDS tenors
            \param recoveryRate Recovery rate
            \param today Valuation date
            \return ISDA bootstrapped spread curve
        */
        static TCurve* bootstrapDefaultCurve(const Handle<YieldTermStructure>& yieldCurve,
                                             const std::vector<Real>& cdsQuotes,
                                             const std::vector<Period>& tenors,
                                             Real recoveryRate,
                                             const Date& today);

        //! Free an ISDA TCurve
        /*! \param curve ISDA TCurve to free
         */
        static void freeTCurve(TCurve* curve);

        //! Convert QuantLib date to ISDA TDate
        /*! \param date QuantLib date
            \return ISDA TDate
        */
        static long dateToTDate(const Date& date);

        //! Convert ISDA TDate to QuantLib date
        /*! \param tDate ISDA TDate
            \return QuantLib Date
        */
        static Date tDateToDate(long tDate);

      private:
        // Helper functions for curve conversion
        static std::vector<Date> buildCurveDates(const Handle<YieldTermStructure>& curve,
                                                 const Date& baseDate);

        static std::vector<Real> extractDiscountFactors(const Handle<YieldTermStructure>& curve,
                                                        const std::vector<Date>& dates);

        static std::vector<Real>
        extractHazardRates(const Handle<DefaultProbabilityTermStructure>& curve,
                           const std::vector<Date>& dates);
    };

}

#endif // QL_ENABLE_ISDA_CDS

#endif // quantlib_isda_cds_interface_hpp