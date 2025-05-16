/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2004, 2009 Ferdinando Ametrano
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006 StatPro Italia srl

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

/*! \file yieldtermstructure.hpp
    \brief Interest-rate term structure
*/

#ifndef quantlib_yield_term_structure_hpp
#define quantlib_yield_term_structure_hpp

#include "ql/indexes/interestrateindex.hpp"
#include <ql/termstructure.hpp>
#include <ql/interestrate.hpp>
#include <ql/quote.hpp>
#include <vector>

namespace QuantLib {

    //! Interest-rate term structure
    /*! This abstract class defines the interface of concrete
        interest rate structures which will be derived from this one.

        \ingroup yieldtermstructures

        \test observability against evaluation date changes is checked.
    */
    class YieldTermStructure : public TermStructure {
      public:
        /*! \name Constructors
            See the TermStructure documentation for issues regarding
            constructors.
        */
        //@{
        explicit YieldTermStructure(const DayCounter& dc = DayCounter());
        YieldTermStructure(const Date& referenceDate,
                           const Calendar& cal = Calendar(),
                           const DayCounter& dc = DayCounter(),
                           std::vector<Handle<Quote> > jumps = {},
                           const std::vector<Date>& jumpDates = {});
        YieldTermStructure(Natural settlementDays,
                           const Calendar& cal,
                           const DayCounter& dc = DayCounter(),
                           std::vector<Handle<Quote> > jumps = {},
                           const std::vector<Date>& jumpDates = {});
        YieldTermStructure(const Date& referenceDate,
                           const ext::shared_ptr<InterestRateIndex>& index,
                           std::vector<Handle<Quote>> jumps = {},
                           const std::vector<Date>& jumpDates = {});
        YieldTermStructure(const ext::shared_ptr<InterestRateIndex>& index,
                           std::vector<Handle<Quote>> jumps = {},
                           const std::vector<Date>& jumpDates = {});

        //@}

        /*! \name Discount factors

            These methods return the discount factor from a given date or time
            to the reference date.  In the latter case, the time is calculated
            as a fraction of year from the reference date.
        */
        //@{
        DiscountFactor discount(const Date& d,
                                bool extrapolate = false) const;
        /*! The same day-counting rule used by the term structure
            should be used for calculating the passed time t.
        */
        DiscountFactor discount(Time t,
                                bool extrapolate = false) const;
        //@}

        /*! \name Zero-yield rates

            These methods return the implied zero-yield rate for a
            given date or time.  In the former case, the time is
            calculated as a fraction of year from the reference date.
        */
        //@{
        /*! The resulting interest rate has the required daycounting
            rule.
        */
        virtual InterestRate zeroRate(const Date& d,
                              const DayCounter& resultDayCounter,
                              Compounding comp,
                              Frequency freq = Annual,
                              bool extrapolate = false) const;

        /*! The resulting interest rate has the same day-counting rule
            used by the term structure. The same rule should be used
            for calculating the passed time t.
        */
        virtual InterestRate zeroRate(Time t,
                              Compounding comp,
                              Frequency freq = Annual,
                              bool extrapolate = false) const;
        //@}

        /*! \name Forward rates

            These methods returns the forward interest rate between two dates
            or times.  In the former case, times are calculated as fractions
            of year from the reference date.

            If both dates (times) are equal the instantaneous forward rate is
            returned.
        */
        //@{
        /*! The resulting interest rate has the required day-counting
            rule.
        */
        virtual InterestRate forwardRate(const Date& d1,
                                 const Date& d2,
                                 const DayCounter& resultDayCounter,
                                 Compounding comp,
                                 Frequency freq = NoFrequency,
                                 bool extrapolate = false) const;
        /*! The resulting interest rate has the required day-counting
            rule.
            \warning dates are not adjusted for holidays
        */
        virtual InterestRate forwardRate(const Date& d,
                                 const Period& p,
                                 const DayCounter& resultDayCounter,
                                 Compounding comp,
                                 Frequency freq = NoFrequency,
                                 bool extrapolate = false) const;

        /*! The resulting interest rate has the same day-counting rule
            used by the term structure. The same rule should be used
            for calculating the passed times t1 and t2.
        */
        virtual InterestRate forwardRate(Time t1,
                                 Time t2,
                                 Compounding comp,
                                 Frequency freq = Annual,
                                 bool extrapolate = false) const;
        //@}
        virtual InterestRate termForwardRate(Time t, bool extrapolate = false) const;

        virtual InterestRate termForwardRate(Date date, bool extrapolate = false) const;

        virtual InterestRate termForwardRate(Date date,
                                             const DayCounter& resultDayCounter,
                                             Compounding resultCompounding = Simple,
                                             Frequency resultFrequency = NoFrequency,
                                             bool extrapolate = false) const;
        //! \name Jump inspectors
        //@{
        const std::vector<Date>& jumpDates() const;
        const std::vector<Time>& jumpTimes() const;
        //@}

        //! \name Observer interface
        //@{
        void update() override;
        //@}

        bool supportsDiscount() const;
        bool supportsZero() const;
        bool isTermForward() const;

        ext::shared_ptr<InterestRateIndex> index() const;
        Period tenor() const;

    protected:
        /*! \name Calculations

            This method must be implemented in derived classes to
            perform the actual calculations. When it is called,
            range check has already been performed; therefore, it
            must assume that extrapolation is required.
        */
        //@{
        //! discount factor calculation
        virtual DiscountFactor discountImpl(Time) const {
            QL_FAIL("Discount factor not implemented for YieldTermStructure");
        }

        virtual Rate forwardImpl(Time) const {
            QL_FAIL("Forward rate not implemented for YieldTermStructure");
        }
        //@{
        /*! Returns the zero yield rate for the given date calculating it
            from the instantaneous forward rate \f$ f(t) \f$ as
            \f[
            z(t) = \int_0^t f(\tau) d\tau
            \f]

            \warning This default implementation uses an highly inefficient
                     and possibly wildly inaccurate numerical integration.
                     Derived classes should override it if a more efficient
                     implementation is available.
        */
        virtual Rate zeroYieldImpl(Time) const {
            QL_FAIL("Zero rate not implemented for YieldTermStructure");
        }

        //@}
        //@}
        bool supportsDiscount_;
        bool supportsZero_;
        bool isTermForward_;

      private:
        // methods
        void setJumps(const Date& referenceDate);
        // data members
        std::vector<Handle<Quote> > jumps_;
        std::vector<Date> jumpDates_;
        std::vector<Time> jumpTimes_;
        Size nJumps_ = 0;
        Date latestReference_;
        ext::shared_ptr<InterestRateIndex> index_;
    };

    // inline definitions

    inline
    DiscountFactor YieldTermStructure::discount(const Date& d,
                                                bool extrapolate) const {
        return discount(timeFromReference(d), extrapolate);
    }

    inline
    InterestRate YieldTermStructure::forwardRate(const Date& d,
                                                 const Period& p,
                                                 const DayCounter& dayCounter,
                                                 Compounding comp,
                                                 Frequency freq,
                                                 bool extrapolate) const {
        return forwardRate(d, index_->advance(d),
            dayCounter, comp, freq, extrapolate);
    }

    inline InterestRate YieldTermStructure::termForwardRate(Time t, bool extrapolate) const {
        checkRange(t, extrapolate);
        if (isTermForward()) {
            return {this->forwardImpl(t), this->dayCounter(), Simple, NoFrequency};
        } else {
            // TODO Why is this needed, please document
            Integer days = static_cast<Integer>(std::round(t / this->dayCounter().yearFraction(Date(367), Date(368))));
            Date date = this->referenceDate() + days;
            return termForwardRate(date, extrapolate);
        }
    }

    inline InterestRate YieldTermStructure::termForwardRate(Date date, bool extrapolate) const {
        if (isTermForward()) {
            const Time t = timeFromReference(date);
            return termForwardRate(t, extrapolate);
        } else {
            return forwardRate(date, this->tenor(), this->dayCounter(), Simple, NoFrequency, extrapolate);
        }
    }

    inline InterestRate YieldTermStructure::termForwardRate(Date date,
                                                            const DayCounter& resultDayCounter,
                                                            Compounding resultCompounding,
                                                            Frequency resultFrequency,
                                                            bool extrapolate) const {
        checkRange(date, extrapolate);
        if (isTermForward()) {
            const Time t = timeFromReference(date);
            return termForwardRate(t, extrapolate).equivalentRate(resultDayCounter, resultCompounding, resultFrequency, date,
                                index_->advance(date));
        } else {
            return forwardRate(date, this->tenor(), resultDayCounter, resultCompounding,
                               resultFrequency, extrapolate);
        }
    }

    inline const std::vector<Date>& YieldTermStructure::jumpDates() const {
        return this->jumpDates_;
    }

    inline const std::vector<Time>& YieldTermStructure::jumpTimes() const {
        return this->jumpTimes_;
    }

}

#endif
