/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
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

#ifndef quantlib_term_forward_curve_hpp
#define quantlib_term_forward_curve_hpp

#include <ql/termstructures/yield/termforwardratestructure.hpp>

namespace QuantLib {

    //! YieldTermStructure based on interpolation of termed forward rates, cannot calculate discount or zero rates
    /*! \ingroup yieldtermstructures */
    template <class Interpolator>
    class InterpolatedTermForwardCurve : public TermForwardRateStructure, public InterpolatedCurve<Interpolator> {
      public:
        // constructors
        InterpolatedTermForwardCurve(const std::vector<Date>& dates,
                                     const std::vector<Rate>& forwards,
                                     const ext::shared_ptr<IborIndex>& index,
                                     Compounding comp = Simple,
                                     Frequency freq = NoFrequency,
                                     const Interpolator& interpolator = {})
        : TermForwardRateStructure(dates.at(0), index),
          InterpolatedCurve<Interpolator>(std::vector<Time>(), forwards, interpolator), comp_(comp),
          freq_(freq),
          dates_(dates) {
            initialize();
        }

        // Add dates_ and initialize

        InterpolatedTermForwardCurve(const std::vector<Date>& dates,
                                     const std::vector<Rate>& forwards,
                                     const ext::shared_ptr<IborIndex>& index,
                                     const Interpolator& interpolator)
        : TermForwardRateStructure(dates.at(0), index),
          InterpolatedCurve<Interpolator>(std::vector<Time>(), forwards, interpolator),
          comp_(Simple),
          freq_(NoFrequency), dates_(dates) {
            initialize();
        }

        // InterpolatedTermForwardCurve(const std::vector<Date>& dates,
        //                              const std::vector<Rate>& forwards,
        //                              const DayCounter& dayCounter,
        //                              const Compounding comp = Continuous,
        //                              const Frequency freq = NoFrequency,
        //                              const Interpolator& interpolator = {})
        //: InterpolatedForwardCurve<Interpolator>(
        //      dates, forwards, dayCounter, interpolator),
        //  comp_(comp), freq_(freq) {}

        InterestRate forwardRate(const Date& d1,
                                 const Date& d2,
                                 const DayCounter& resultDayCounter,
                                 Compounding comp,
                                 Frequency freq,
                                 bool extrapolate) const override;

        InterestRate forwardRate(
            Time t1, Time t2, Compounding comp, Frequency freq, bool extrapolate) const override;

        InterestRate forwardRate(const Date& d,
                                 const Period& p,
                                 const DayCounter& resultDayCounter,
                                 Compounding resultCompounding = Simple,
                                 Frequency resultFrequency = NoFrequency,
                                 bool extrapolate = false) const override;

        InterestRate termForwardRate(Time t, bool extrapolate = false) const override;

        InterestRate termForwardRate(Date date, bool extrapolate = false) const override;

        InterestRate termForwardRate(Date date,
                                     const DayCounter& resultDayCounter,
                                     Compounding resultCompounding = Simple,
                                     Frequency resultFrequency = NoFrequency,
                                     bool extrapolate = false) const override;
        //! \name TermStructure interface
        //@{
        Date maxDate() const override;
        //@}
        //! \name other inspectors
        //@{
        const std::vector<Time>& times() const;
        const std::vector<Date>& dates() const;
        const std::vector<Real>& data() const;
        const std::vector<Rate>& forwards() const;
        std::vector<std::pair<Date, Real>> nodes() const;
        //@}

      protected:
        explicit InterpolatedTermForwardCurve(const ext::shared_ptr<IborIndex>&,
                                              const Interpolator& interpolator = {});

        InterpolatedTermForwardCurve(const Date& referenceDate,
                                     const ext::shared_ptr<IborIndex>& index,
                                     const Interpolator& interpolator);

        Compounding comp_;
        Frequency freq_;

        Rate forwardImpl(Time t) const override;

        mutable std::vector<Date> dates_;

      private:
        void initialize();
    };

    // inline definitions

    template <class T>
    inline Date InterpolatedTermForwardCurve<T>::maxDate() const {
        if (this->maxDate_ != Date())
            return this->maxDate_;
        return this->dates_.back();
    }

    template <class T>
    inline const std::vector<Time>& InterpolatedTermForwardCurve<T>::times() const {
        return this->times_;
    }

    template <class T>
    inline const std::vector<Date>& InterpolatedTermForwardCurve<T>::dates() const {
        return this->dates_;
    }

    template <class T>
    inline const std::vector<Real>& InterpolatedTermForwardCurve<T>::data() const {
        return this->data_;
    }

    template <class T>
    inline const std::vector<Rate>& InterpolatedTermForwardCurve<T>::forwards() const {
        return this->data_;
    }

    template <class T>
    inline std::vector<std::pair<Date, Real>> InterpolatedTermForwardCurve<T>::nodes() const {
        std::vector<std::pair<Date, Real>> results(this->dates_.size());
        for (Size i = 0; i < this->dates_.size(); ++i)
            results[i] = std::make_pair(this->dates_[i], this->data_[i]);
        return results;
    }

    #ifndef __DOXYGEN__

    #endif
    template <class T>
    Rate InterpolatedTermForwardCurve<T>::forwardImpl(Time t) const {
        if (t <= this->times_.back())
            return this->interpolation_(t, true);

        // flat fwd extrapolation
        return this->data_.back();
    }

    template <class T>
    InterestRate InterpolatedTermForwardCurve<T>::forwardRate(const Date& d1,
                                                              const Date& d2,
                                                              const DayCounter& resultDayCounter,
                                                              Compounding comp,
                                                              Frequency freq,
                                                              bool extrapolate) const {
        return termForwardRate(d1, extrapolate).equivalentRate(resultDayCounter, comp, freq, d1, d2);
    }

    template <class T>
    InterestRate InterpolatedTermForwardCurve<T>::forwardRate(const Time t1,
                                                              const Time t2,
                                                              Compounding comp,
                                                              Frequency freq,
                                                              bool extrapolate) const {
        QL_FAIL("Cannot get forward rate for a termed forward curve based on time only, pass dates so that period length can be inferred.");
    }

    template <class Interpolator>
    InterestRate InterpolatedTermForwardCurve<Interpolator>::forwardRate(const Date& d,
        const Period& p,
        const DayCounter& resultDayCounter,
        Compounding resultCompounding,
        Frequency resultFrequency,
        bool extrapolate) const {
        QL_REQUIRE(p == this->index()->tenor(), "Period " << p
                                                          << " not compatible with curve period "
                                                          << this->index()->tenor() << ".");
        return termForwardRate(d, resultDayCounter, resultCompounding, resultFrequency,
                               extrapolate);
    }

    template <class Interpolator>
    InterestRate InterpolatedTermForwardCurve<Interpolator>::termForwardRate(
        Time t, bool extrapolate) const {
        TermStructure::checkRange(t, extrapolate);
        return {this->forwardImpl(t), this->index()->dayCounter(), comp_, freq_};
    }

    template <class Interpolator>
    InterestRate
    InterpolatedTermForwardCurve<Interpolator>::termForwardRate(Date date,
                                                                bool extrapolate) const {
        const Time t = TermStructure::timeFromReference(date);
        return termForwardRate(t, extrapolate);
    }

    template <class Interpolator>
    InterestRate
    InterpolatedTermForwardCurve<Interpolator>::termForwardRate(Date date,
                                                                const DayCounter& resultDayCounter,
                                                                Compounding resultCompounding,
                                                                Frequency resultFrequency,
                                                                bool extrapolate) const {
        TermStructure::checkRange(date, extrapolate);
        const Time t = this->timeFromReference(date);
        return InterestRate(this->forwardImpl(t), this->index()->dayCounter(), comp_, freq_)
            .equivalentRate(resultDayCounter, resultCompounding, resultFrequency, date, this->index()->maturityDate(date));
    }

    template <class Iterator>
    inline InterpolatedTermForwardCurve<Iterator>::InterpolatedTermForwardCurve(
        const ext::shared_ptr<IborIndex>& index, const Iterator& interpolator)
    : TermForwardRateStructure(index), InterpolatedCurve<Iterator>(interpolator), comp_(Simple),
      freq_(NoFrequency) {}

    template <class T>
    InterpolatedTermForwardCurve<T>::InterpolatedTermForwardCurve(
        const Date& referenceDate, const ext::shared_ptr<IborIndex>& index, const T& interpolator)
    : TermForwardRateStructure(referenceDate, index), InterpolatedCurve<T>(interpolator), comp_(Simple),
      freq_(NoFrequency) {}


    template <class T>
    void InterpolatedTermForwardCurve<T>::initialize() {
        QL_REQUIRE(dates_.size() >= T::requiredPoints, "not enough input dates given");
        QL_REQUIRE(this->data_.size() == dates_.size(), "dates/data count mismatch");
        YieldTermStructure::supportsDiscount_ = false;
        YieldTermStructure::supportsZero_ = false;
        YieldTermStructure::isTermForward_ = true;

        this->setupTimes(dates_, dates_[0], dayCounter());
        this->setupInterpolation();
        this->interpolation_.update();
    }


}

#endif
