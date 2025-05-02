/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2015 Peter Caspers

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

/*! \file swaptionhelper.hpp
    \brief Swaption calibration helper
*/

#ifndef quantlib_swaption_calibration_helper_hpp
#define quantlib_swaption_calibration_helper_hpp

#include <ql/models/calibrationhelper.hpp>
#include <ql/instruments/swaption.hpp>
#include <ql/termstructures/volatility/volatilitytype.hpp>
#include <ql/cashflows/rateaveraging.hpp>

namespace QuantLib {

    //! calibration helper for interest-rate swaptions
    /*! \warning passing an overnight index to the constructor will
                 result in an overnight-indexed swap being built, but
                 model-based engines will treat it as a vanilla swap.
                 This is at best a decent proxy, at worst simply wrong.
                 Use with caution.
    */
    class SwaptionHelper : public BlackCalibrationHelper {
      public:
        SwaptionHelper(const Period& maturity,
                       const Period& length,
                       const Handle<Quote>& volatility,
                       ext::shared_ptr<IborIndex> index,
                       const Period& fixedLegTenor,
                       DayCounter fixedLegDayCounter,
                       DayCounter floatingLegDayCounter,
                       Handle<YieldTermStructure> termStructure,
                       CalibrationErrorType errorType = RelativePriceError,
                       Real strike = Null<Real>(),
                       Real nominal = 1.0,
                       VolatilityType type = ShiftedLognormal,
                       Real shift = 0.0,
                       Natural settlementDays = Null<Size>(),
                       RateAveraging::Type averagingMethod = RateAveraging::Compound);

        SwaptionHelper(const Date& exerciseDate,
                       const Period& length,
                       const Handle<Quote>& volatility,
                       ext::shared_ptr<IborIndex> index,
                       const Period& fixedLegTenor,
                       DayCounter fixedLegDayCounter,
                       DayCounter floatingLegDayCounter,
                       Handle<YieldTermStructure> termStructure,
                       CalibrationErrorType errorType = RelativePriceError,
                       Real strike = Null<Real>(),
                       Real nominal = 1.0,
                       VolatilityType type = ShiftedLognormal,
                       Real shift = 0.0,
                       Natural settlementDays = Null<Size>(),
                       RateAveraging::Type averagingMethod = RateAveraging::Compound);

        SwaptionHelper(const Date& exerciseDate,
                       const Date& endDate,
                       const Handle<Quote>& volatility,
                       ext::shared_ptr<IborIndex> index,
                       const Period& fixedLegTenor,
                       DayCounter fixedLegDayCounter,
                       DayCounter floatingLegDayCounter,
                       Handle<YieldTermStructure> termStructure,
                       CalibrationErrorType errorType = RelativePriceError,
                       Real strike = Null<Real>(),
                       Real nominal = 1.0,
                       VolatilityType type = ShiftedLognormal,
                       Real shift = 0.0,
                       Natural settlementDays = Null<Size>(),
                       RateAveraging::Type averagingMethod = RateAveraging::Compound);

        void addTimesTo(std::list<Time>& times) const override;
        Real modelValue() const override;
        Real blackPrice(Volatility volatility) const override;

        // populated in call to blackPrice():
        Real timeToExpiry() const { return timeToExpiry_; }
        Real swapLength() const { return swapLength_; }
        Real strike() const { return strike_; }
        Real atmForward() const { return atmForward_; }
        Real annuity() const { return annuity_; };
        Real vega() const { return vega_; }
        Real stdDev() const { return stdDev_; }

        const ext::shared_ptr<FixedVsFloatingSwap>& underlying() const {
            calculate();
            return swap_;
        }
        /*! \deprecated Use the SwaptionHelper::underlying method instead.
                        Deprecated in version 1.34.
        */
        [[deprecated("Use the SwaptionHelper::underlying method instead")]]
        ext::shared_ptr<VanillaSwap> underlyingSwap() const {
            calculate();
            auto vanilla = ext::dynamic_pointer_cast<VanillaSwap>(swap_);
            QL_REQUIRE(vanilla, "underlying is not a vanilla swap");
            return vanilla;
        }
        ext::shared_ptr<Swaption> swaption() const { calculate(); return swaption_; }

      private:
        void performCalculations() const override;
        ext::shared_ptr<FixedVsFloatingSwap> makeSwap(Schedule fixedSchedule,
                                                      Schedule floatSchedule,
                                                      Rate exerciseRate,
                                                      Swap::Type type) const;
        mutable Date exerciseDate_, endDate_;
        const Period maturity_, length_, fixedLegTenor_;
        const ext::shared_ptr<IborIndex> index_;
        const Handle<YieldTermStructure> termStructure_;
        const DayCounter fixedLegDayCounter_, floatingLegDayCounter_;
        const Real strike_, nominal_;
        const Natural settlementDays_;
        const RateAveraging::Type averagingMethod_;
        mutable Rate exerciseRate_;
        mutable ext::shared_ptr<FixedVsFloatingSwap> swap_;
        mutable ext::shared_ptr<Swaption> swaption_;
        mutable Real timeToExpiry_ = Null<Real>();
        mutable Real swapLength_ = Null<Real>();
        mutable Real atmForward_ = Null<Real>();
        mutable Real annuity_ = Null<Real>();
        mutable Real vega_ = Null<Real>();
        mutable Real stdDev_ = Null<Real>();
    };
}

#endif
