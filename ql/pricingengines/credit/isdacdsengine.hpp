/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2014 Jose Aparicio
 Copyright (C) 2014 Peter Caspers

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <https://www.quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#ifndef quantlib_isda_cds_engine_hpp
#define quantlib_isda_cds_engine_hpp

#include <ql/instruments/creditdefaultswap.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/termstructures/defaulttermstructure.hpp>
#include <ql/optional.hpp>

namespace QuantLib {

    /*! References:

        [1] The Pricing and Risk Management of Credit Default Swaps, with a
            Focus on the ISDA Model,
            OpenGamma Quantitative Research, Version as of 15-Oct-2013

        [2] ISDA CDS Standard Model Proposed Numerical Fix \ Thursday,
            November 15, 2012, Markit

        [3] Markit Interest Rate Curve XML Specifications,
            Version 1.16, Tuesday, 15 October 2013

    */
    //! Isda CDS engine base
    //! \ingroup engines
    class IsdaCdsEngineBase {
    public:
        enum NumericalFix {
            None,  // as in [1] footnote 26 (i.e. 10^{-50} is added to
                   // denominators $f_i+h_i$$)
            Taylor // as in [2] i.e. for $f_i+h_i < 10^{-4}$ a Taylor expansion
                   // is used to avoid zero denominators
        };

        enum AccrualBias {
            HalfDayBias, // as in [1] formula (50), second (error) term is
                         // included
            NoBias // as in [1], but second term in formula (50) is not included
        };

        enum ForwardsInCouponPeriod {
            Flat, // as in [1], formula (52), second (error) term is included
            Piecewise // as in [1], but second term in formula (52) is not
                      // included
        };
        IsdaCdsEngineBase(const Handle<YieldTermStructure>& discountCurve,
                        const Handle<DefaultProbabilityTermStructure>& probability,
                        ext::optional<bool> includeSettlementDateFlows,
                        NumericalFix numericalFix,
                        AccrualBias accrualBias,
                        ForwardsInCouponPeriod forwardsInCouponPeriod)
            : discountCurve_(discountCurve), probability_(probability),includeSettlementDateFlows_(includeSettlementDateFlows),
            numericalFix_(numericalFix), accrualBias_(accrualBias), forwardsInCouponPeriod_(forwardsInCouponPeriod) {}
        virtual ~IsdaCdsEngineBase() {}

    protected:
        virtual Real survivalProbability(const Date& d) const = 0;
        virtual Real defaultProbability(const Date& d1, const Date& d2) const = 0;
        virtual Real expectedLoss(const Date& defaultDate, const Date& d1, const Date& d2, const Real notional) const = 0;
        virtual Real claimLoss(const Date& defaultDate, const Real notional) const = 0;
        void calculate(const Date& refDate, const CreditDefaultSwap::arguments& arguments,
                       CreditDefaultSwap::results& results) const;

        Handle<YieldTermStructure> discountCurve_;
        Handle<DefaultProbabilityTermStructure> probability_;
        ext::optional<bool> includeSettlementDateFlows_;
        NumericalFix numericalFix_;
        AccrualBias accrualBias_;
        ForwardsInCouponPeriod forwardsInCouponPeriod_;

    };

    //! Isda CDS engine
    //! \ingroup engines
    class IsdaCdsEngine : public CreditDefaultSwap::engine, public IsdaCdsEngineBase {
      public:
        /*! According to [1] the settings for the flags
            AccrualBias / ForwardsInCouponPeriod corresponding
            to the standard model implementation C code are

            prior 1.8.2    HalfDayBias / Flat
            1.8.2          NoBias / Flat

            The theoretical correct setting would be NoBias / Piecewise

            Todo: Clarify in which version of the standard model
            implementation C code the numerical problem of zero denominators
            is solved and how exactly.
        */

        /*! Constructor where the client code is responsible for providing a
            default curve and an interest rate curve compliant with the ISDA
            specifications.

            To be precisely consistent with the ISDA specification
                bool IborCoupon::Settings::usingAtParCoupons();
            must be true. This is not checked in order not to
            kill the engine completely in this case.

            Furthermore, the ibor index in the swap rate helpers should not
            provide the evaluation date's fixing.
        */

        IsdaCdsEngine(const Handle<DefaultProbabilityTermStructure>& probability,
                      Real recoveryRate,
                      const Handle<YieldTermStructure>& discountCurve,
                      const ext::optional<bool> includeSettlementDateFlows = ext::nullopt,
                      const NumericalFix numericalFix = Taylor,
                      const AccrualBias accrualBias = HalfDayBias,
                      const ForwardsInCouponPeriod forwardsInCouponPeriod = Piecewise);

        // Handle<YieldTermStructure> isdaRateCurve() const { return discountCurve_; }
        // Handle<DefaultProbabilityTermStructure> isdaCreditCurve() const { return probability_; }

        void calculate() const override;

      protected:
        virtual Real survivalProbability(const Date& d) const override;
        virtual Real defaultProbability(const Date& d1, const Date& d2) const override;
        virtual Real expectedLoss(const Date& defaultDate, const Date& d1, const Date& d2, const Real notional) const override;
        virtual Real claimLoss(const Date& defaultDate, const Real notional) const override;
        
        // mutable Handle<DefaultProbabilityTermStructure> probability_;
        mutable Real recoveryRate_;

    };
}

#endif
