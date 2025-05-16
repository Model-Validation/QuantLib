/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004 StatPro Italia srl
 Copyright (C) 2009 Ferdinando Ametrano

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

#include <ql/termstructures/yield/forwardstructure.hpp>

namespace QuantLib {

    ForwardRateStructure::ForwardRateStructure(const DayCounter& dayCounter)
    : YieldTermStructure(dayCounter) {}

    ForwardRateStructure::ForwardRateStructure(const Date& referenceDate,
                                               const Calendar& cal,
                                               const DayCounter& dayCounter,
                                               const std::vector<Handle<Quote>>& jumps,
                                               const std::vector<Date>& jumpDates)
    : YieldTermStructure(referenceDate, cal, dayCounter, jumps, jumpDates) {}

    ForwardRateStructure::ForwardRateStructure(Natural settlementDays,
                                               const Calendar& calendar,
                                               const DayCounter& dayCounter,
                                               const std::vector<Handle<Quote>>& jumps,
                                               const std::vector<Date>& jumpDates)
    : YieldTermStructure(settlementDays, calendar, dayCounter, jumps, jumpDates) {}

    Rate ForwardRateStructure::zeroYieldImpl(Time t) const {
        if (t == 0.0)
            return forwardImpl(0.0);
        // implement smarter integration if plan to use the following code
        Rate sum = 0.5*forwardImpl(0.0);
        Size N = 1000;
        Time dt = t/static_cast<Time>(N);
        for (Time i=dt; i<t; i+=dt)
            sum += forwardImpl(i);
        sum += 0.5*forwardImpl(t);
        return static_cast<Rate>(sum * dt / t);
    }

}
