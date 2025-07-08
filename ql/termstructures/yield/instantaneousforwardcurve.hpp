/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2005, 2006, 2007, 2008, 2009 StatPro Italia srl
 Copyright (C) 2009, 2015 Ferdinando Ametrano
 Copyright (C) 2015 Paolo Mazzocchi
  Copyright (C) 2024 SEB AB Sverrir Thorvaldsson

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

/*! \file termforwardcurve.hpp
    \brief interpolated term-forward-rate structure
*/

#ifndef quantlib_instantaneous_forward_curve_hpp
#define quantlib_instantaneous_forward_curve_hpp

#include <ql/termstructures/yield/forwardstructure.hpp>
#include <ql/termstructures/yield/forwardcurve.hpp>

// TODO This is a bad and misleading name, should be InterpolatedInstantaneousForwardCurve
namespace QuantLib {
    template <class Interpolator>
    class InterpolatedInstantaneousForwardCurve : public InterpolatedForwardCurve<Interpolator> {
    public:
        InterpolatedInstantaneousForwardCurve(const std::vector<Date>& dates,
                                              const std::vector<Rate>& forwards,
                                              const DayCounter& dayCounter,
                                              const Calendar& calendar = Calendar(),
                                              const Interpolator& interpolator = {})
        : InterpolatedForwardCurve<Interpolator>(
              dates, forwards, dayCounter, calendar, interpolator) {}

        [[nodiscard]] InterestRate forwardRate(const Date& d1,
                                               const Date& d2,
                                               const DayCounter& resultDayCounter,
                                               Compounding comp,
                                               Frequency freq,
                                               bool extrapolate) const override;

        [[nodiscard]] InterestRate forwardRate(
            Time t1, Time t2, Compounding comp, Frequency freq, bool extrapolate) const override;

        [[nodiscard]] InterestRate zeroRate(const Date& d,
                                            const DayCounter& resultDayCounter,
                                            Compounding comp,
                                            Frequency freq = Annual,
                                            bool extrapolate = false) const override;

        [[nodiscard]] InterestRate zeroRate(Time t,
                                            Compounding comp,
                                            Frequency freq = Annual,
                                            bool extrapolate = false) const override;


    protected:
        InterpolatedInstantaneousForwardCurve(const DayCounter& dayCounter,
                                     const Interpolator& interpolator = {})
        : InterpolatedForwardCurve<Interpolator>(dayCounter, interpolator) {}

        InterpolatedInstantaneousForwardCurve(const Date& referenceDate,
                                              const DayCounter& dayCounter,
                                              const std::vector<Handle<Quote>>& jumps = {},
                                              const std::vector<Date>& jumpDates = {},
                                              const Interpolator& interpolator = {})
        : InterpolatedForwardCurve<Interpolator>(referenceDate, dayCounter, jumps, jumpDates, interpolator) {}

        InterpolatedInstantaneousForwardCurve(Natural settlementDays,
                                     const Calendar& calendar,
                                              const DayCounter& dayCounter,
                                              const std::vector<Handle<Quote>>& jumps = {},
                                              const std::vector<Date>& jumpDates = {},
                                     const Interpolator& interpolator = {})
        : InterpolatedForwardCurve<Interpolator>(
              settlementDays, calendar, dayCounter, jumps, jumpDates, interpolator) {}

        [[nodiscard]] DiscountFactor discountImpl(Time t) const override;
    };

    #ifndef __DOXYGEN__
    // template definitions

    template <class T>
    DiscountFactor InterpolatedInstantaneousForwardCurve<T>::discountImpl(Time t) const
    {
        return static_cast<DiscountFactor>(std::exp(-this->zeroYieldImpl(t) * t));
    }
    #endif

    template <class T>
    InterestRate
    InterpolatedInstantaneousForwardCurve<T>::zeroRate(const Date& d,
                                                        const DayCounter& resultDayCounter,
                                                        Compounding comp,
                                                        Frequency freq,
                                                        bool extrapolate) const {
        const Time t = ForwardRateStructure::timeFromReference(d);
        return InterestRate(this->zeroYieldImpl(t), this->dayCounter(), Continuous, NoFrequency)
            .equivalentRate(resultDayCounter, comp, freq, this->referenceDate(), d);
    }

    template <class T>
    InterestRate InterpolatedInstantaneousForwardCurve<T>::zeroRate(Time t,
                                                        Compounding comp,
                                                        Frequency freq,
                                                        bool extrapolate) const {
        return InterestRate(this->zeroYieldImpl(t), this->dayCounter(), Continuous, NoFrequency)
            .equivalentRate(comp, freq, t);
    }

    template <class T>
    InterestRate
    InterpolatedInstantaneousForwardCurve<T>::forwardRate(const Date& d1,
                                                           const Date& d2,
                                                           const DayCounter& resultDayCounter,
                                                           Compounding comp,
                                                           Frequency freq,
                                                           bool extrapolate) const {
        //QL_REQUIRE(d1 <= d2, "d2 (" << d2 << ") < d1 (" << d1 << ")");
        //ForwardRateStructure::checkRange(d1, extrapolate);
        const Time t1 = ForwardRateStructure::timeFromReference(d1);
        const Time t2 = ForwardRateStructure::timeFromReference(d2);
        if (d1 == d2) {
            const Real r = this->forwardImpl(t1);
            return InterestRate(r, this->dayCounter(), Continuous, NoFrequency);
        }
        const Real r = (t2 * this->zeroYieldImpl(t2) - t1 * this->zeroYieldImpl(t1)) / (t2 - t1);
        return InterestRate(r, this->dayCounter(), Continuous, NoFrequency)
            .equivalentRate(comp, freq, t2 - t1);
    }

    template <class T>
    InterestRate InterpolatedInstantaneousForwardCurve<T>::forwardRate(
        const Time t1, const Time t2, Compounding comp, Frequency freq, bool extrapolate) const {
        QL_REQUIRE(t2 >= t1, "t2 (" << t2 << ") < t1 (" << t1 << ")");
        ForwardRateStructure::checkRange(t1, extrapolate);
        if (t2 == t1) {
            const Real r = this->forwardImpl(t1);
            return InterestRate(r, this->dayCounter(), Continuous, NoFrequency);
        }
        const Real r = (t2 * this->zeroYieldImpl(t2) - t1 * this->zeroYieldImpl(t1)) / (t2-t1);
        return InterestRate(r, this->dayCounter(), Continuous, NoFrequency).equivalentRate(comp, freq, t2 - t1);
    }
}

#endif
