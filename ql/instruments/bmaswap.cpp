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
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include <ql/cashflows/averagebmacoupon.hpp>
#include <ql/cashflows/iborcoupon.hpp>
#include <ql/cashflows/overnightindexedcoupon.hpp>
#include <ql/instruments/bmaswap.hpp>

namespace QuantLib {

    BMASwap::BMASwap(Type type,
                     Real nominal,
                     // Index leg
                     Schedule indexSchedule,
                     Real indexFraction,
                     Spread indexSpread,
                     const ext::shared_ptr<IborIndex>& index,
                     const DayCounter& indexDayCount,
                     // BMA leg
                     Schedule bmaSchedule,
                     const ext::shared_ptr<BMAIndex>& bmaIndex,
                     const DayCounter& bmaDayCount,
                     // Payment adjustments
                     Calendar indexPaymentCalendar,
                     BusinessDayConvention indexPaymentConvention,
                     Natural indexPaymentLag,
                     Calendar bmaPaymentCalendar,
                     BusinessDayConvention bmaPaymentConvention,
                     Natural bmaPaymentLag,
                     // overnight conventions
                     Natural overnightLockoutDays,
                     bool telescopicValueDates)
    : Swap(2), type_(type), nominal_(nominal), indexFraction_(indexFraction),
      indexSpread_(indexSpread), indexPaymentCalendar_(indexPaymentCalendar),
      indexPaymentConvention_(indexPaymentConvention), indexPaymentLag_(indexPaymentLag),
      bmaPaymentCalendar_(bmaPaymentCalendar), bmaPaymentConvention_(bmaPaymentConvention),
      bmaPaymentLag_(bmaPaymentLag), overnightLockoutDays_(overnightLockoutDays) {

        auto ois = QuantLib::ext::dynamic_pointer_cast<OvernightIndex>(index);

        Calendar effectiveIndexPaymentCalendar =
            indexPaymentCalendar_.empty() ?
                (indexSchedule.calendar().empty() ? NullCalendar() : indexSchedule.calendar()) :
                indexPaymentCalendar_;
        Calendar effectiveBmaPaymentCalendar =
            bmaPaymentCalendar_.empty() ?
                (bmaSchedule.calendar().empty() ? NullCalendar() : bmaSchedule.calendar()) :
                bmaPaymentCalendar_;

        if (ois) {
            legs_[0] = OvernightLeg(std::move(indexSchedule), ois)
                           .withNotionals(nominal)
                           .withPaymentDayCounter(indexDayCount)
                           .withPaymentAdjustment(indexPaymentConvention_)
                           .withPaymentLag(indexPaymentLag_)
                           .withPaymentCalendar(effectiveIndexPaymentCalendar)
                           .withTelescopicValueDates(telescopicValueDates)
                           .withGearings(indexFraction)
                           .withSpreads(indexSpread)
                           .withLockoutDays(overnightLockoutDays);

        } else {
            legs_[0] = IborLeg(std::move(indexSchedule), index)
                           .withNotionals(nominal)
                           .withPaymentDayCounter(indexDayCount)
                           .withPaymentAdjustment(indexPaymentConvention_)
                           .withPaymentLag(indexPaymentLag_)
                           .withPaymentCalendar(effectiveIndexPaymentCalendar)
                           .withFixingDays(index->fixingDays())
                           .withGearings(indexFraction)
                           .withSpreads(indexSpread);
        }


        legs_[1] = AverageBMALeg(std::move(bmaSchedule), bmaIndex)
                       .withNotionals(nominal)
                       .withPaymentDayCounter(bmaDayCount)
                       .withPaymentAdjustment(bmaPaymentConvention_)
                       .withPaymentLag(bmaPaymentLag_)
                       .withPaymentCalendar(effectiveIndexPaymentCalendar);


        for (Size j = 0; j < 2; ++j) {
            for (auto& i : legs_[j])
                registerWith(i);
        }

        switch (type_) {
            case Payer:
                payer_[0] = +1.0;
                payer_[1] = -1.0;
                break;
            case Receiver:
                payer_[0] = -1.0;
                payer_[1] = +1.0;
                break;
            default:
                QL_FAIL("Unknown BMA-swap type");
        }
    }

    Real BMASwap::indexFraction() const {
        return indexFraction_;
    }

    Spread BMASwap::indexSpread() const {
        return indexSpread_;
    }

    Real BMASwap::nominal() const {
        return nominal_;
    }

    Swap::Type BMASwap::type() const {
        return type_;
    }

    const Leg& BMASwap::indexLeg() const {
        return legs_[0];
    }

    const Leg& BMASwap::bmaLeg() const {
        return legs_[1];
    }


    Real BMASwap::indexLegBPS() const {
        calculate();
        QL_REQUIRE(legBPS_[0] != Null<Real>(), "result not available");
        return legBPS_[0];
    }

    Real BMASwap::indexLegNPV() const {
        calculate();
        QL_REQUIRE(legNPV_[0] != Null<Real>(), "result not available");
        return legNPV_[0];
    }

    Real BMASwap::fairIndexFraction() const {
        static Spread basisPoint = 1.0e-4;

        Real spreadNPV = (indexSpread_ / basisPoint) * indexLegBPS();
        Real pureIndexNPV = indexLegNPV() - spreadNPV;
        QL_REQUIRE(pureIndexNPV != 0.0, "result not available (null index NPV)");
        return -indexFraction_ * (bmaLegNPV() + spreadNPV) / pureIndexNPV;
    }

    Spread BMASwap::fairIndexSpread() const {
        static Spread basisPoint = 1.0e-4;

        return indexSpread_ - NPV() / (indexLegBPS() / basisPoint);
    }

    Real BMASwap::bmaLegBPS() const {
        calculate();
        QL_REQUIRE(legBPS_[1] != Null<Real>(), "result not available");
        return legBPS_[1];
    }

    Real BMASwap::bmaLegNPV() const {
        calculate();
        QL_REQUIRE(legNPV_[1] != Null<Real>(), "result not available");
        return legNPV_[1];
    }

}
