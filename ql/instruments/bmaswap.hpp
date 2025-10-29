/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2007 Roland Lichters
 Copyright (C) 2007 StatPro Italia srl

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

/*! \file bmaswap.hpp
    \brief swap paying Libor against BMA coupons
*/

#ifndef quantlib_bma_swap_hpp
#define quantlib_bma_swap_hpp

#include <ql/instruments/swap.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/indexes/bmaindex.hpp>

namespace QuantLib {

    //! swap paying Libor against BMA coupons
    class BMASwap : public Swap {
      public:
        /*! In this constructor, the type (Payer or Receiver) refers
            to the BMA leg.
        */
        BMASwap(Type type,
                Real nominal,
                // index leg
                Schedule indexSchedule,
                Rate indexFraction,
                Rate indexSpread,
                const ext::shared_ptr<IborIndex>& index,
                const DayCounter& indexDayCount,
                // BMA leg
                Schedule bmaSchedule,
                const ext::shared_ptr<BMAIndex>& bmaIndex,
                const DayCounter& bmaDayCount,
                // Payment adjustments
                Calendar indexPaymentCalendar = Calendar(),
                BusinessDayConvention indexPaymentConvention = Following,
                Natural indexPaymentLag = 0,
                Calendar bmaPaymentCalendar = Calendar(),
                BusinessDayConvention bmaPaymentConvention = Following,
                Natural bmaPaymentLag = 0,
                // overnight conventions
                Natural overnightLoackoutDays = 0,
                bool telescopicValueDates = false);

        //! \name Inspectors
        //@{
        Real indexFraction() const;
        Spread indexSpread() const;
        Real nominal() const;
        //! "Payer" or "Receiver" refers to the BMA leg
        Type type() const;
        const Leg& bmaLeg() const;
        const Leg& indexLeg() const;
        //@}

        //! \name Results
        //@{
        Real indexLegBPS() const;
        Real indexLegNPV() const;
        Rate fairIndexFraction() const;
        Spread fairIndexSpread() const;

        Real bmaLegBPS() const;
        Real bmaLegNPV() const;
        //@}

      private:
        Type type_;
        Real nominal_;
        Rate indexFraction_;
        Rate indexSpread_;
        Calendar indexPaymentCalendar_;
        BusinessDayConvention indexPaymentConvention_;
        Natural indexPaymentLag_;
        Calendar bmaPaymentCalendar_;
        BusinessDayConvention bmaPaymentConvention_;
        Natural bmaPaymentLag_;
        Natural overnightLockoutDays_;
        bool telescopicValueDates_;
    };

}

#endif
