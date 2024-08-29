#ifndef quantlib_rt_curve_hpp
#define quantlib_rt_curve_hpp

#include <ql/math/comparison.hpp>
#include <ql/math/interpolations/backwardflatinterpolation.hpp>
#include <ql/termstructures/interpolatedcurve.hpp>
#include <ql/termstructures/yield/rtstructure.hpp>
#include <utility>

namespace QuantLib {

    //! YieldTermStructure based on interpolation of RT values
    /*! \ingroup yieldtermstructures */
    template <class Interpolator>
    class InterpolatedRateTimeCurve : public RateTimeStructure,
                                      public InterpolatedCurve<Interpolator> {
      public:
        // constructor
        InterpolatedRateTimeCurve(const std::vector<Date>& dates,
                                  const std::vector<Rate>& rate_times,
                                  const DayCounter& dayCounter,
                                  const Calendar& cal = Calendar(),
                                  const std::vector<Handle<Quote>>& jumps = {},
                                  const std::vector<Date>& jumpDates = {},
                                  const Interpolator& interpolator = {});
        InterpolatedRateTimeCurve(const std::vector<Date>& dates,
                                  const std::vector<Rate>& rate_times,
                                  const DayCounter& dayCounter,
                                  const Calendar& calendar,
                                  const Interpolator& interpolator);
        InterpolatedRateTimeCurve(const std::vector<Date>& dates,
                                  const std::vector<Rate>& rate_times,
                                  const DayCounter& dayCounter,
                                  const Interpolator& interpolator);
        //! \name TermStructure interface
        //@{
        Date maxDate() const override;
        //@}
        //! \name other inspectors
        //@{
        const std::vector<Time>& times() const;
        const std::vector<Date>& dates() const;
        const std::vector<Real>& data() const;
        const std::vector<Rate>& rate_times() const;
        std::vector<std::pair<Date, Real>> nodes() const;
        //@}

      protected:
        explicit InterpolatedRateTimeCurve(const DayCounter&,
                                           const Interpolator& interpolator = {});
        InterpolatedRateTimeCurve(const Date& referenceDate,
                                  const DayCounter&,
                                  const std::vector<Handle<Quote>>& jumps = {},
                                  const std::vector<Date>& jumpDates = {},
                                  const Interpolator& interpolator = {});
        InterpolatedRateTimeCurve(Natural settlementDays,
                                  const Calendar&,
                                  const DayCounter&,
                                  const std::vector<Handle<Quote>>& jumps = {},
                                  const std::vector<Date>& jumpDates = {},
                                  const Interpolator& interpolator = {});

        //! \name RateTimeStructure implementation
        //@{
        Rate zeroYieldImpl(Time t) const override;
        //@}
        mutable std::vector<Date> dates_;

      private:
        void initialize();
    };

    //! Term structure based on flat interpolation of rate*time values
    /*! \ingroup yieldtermstructures */

    typedef InterpolatedRateTimeCurve<BackwardFlat> RateTimeCurve;


    // inline definitions

    template <class T>
    inline Date InterpolatedRateTimeCurve<T>::maxDate() const {
        if (this->maxDate_ != Date())
            return this->maxDate_;
        return dates_.back();
    }

    template <class T>
    inline const std::vector<Time>& InterpolatedRateTimeCurve<T>::times() const {
        return this->times_;
    }

    template <class T>
    inline const std::vector<Date>& InterpolatedRateTimeCurve<T>::dates() const {
        return dates_;
    }

    template <class T>
    inline const std::vector<Real>& InterpolatedRateTimeCurve<T>::data() const {
        return this->data_;
    }

    template <class T>
    inline const std::vector<Rate>& InterpolatedRateTimeCurve<T>::rate_times() const {
        return this->data_;
    }

    template <class T>
    inline std::vector<std::pair<Date, Real>> InterpolatedRateTimeCurve<T>::nodes() const {
        std::vector<std::pair<Date, Real>> results(dates_.size());
        for (Size i = 0; i < dates_.size(); ++i)
            results[i] = std::make_pair(dates_[i], this->data_[i]);
        return results;
    }

#ifndef __DOXYGEN__

    // template definitions


    template <class T>
    Rate InterpolatedRateTimeCurve<T>::zeroYieldImpl(Time t) const {
        if (t == 0.0)
            return Rate(this->interpolation_.derivative(0.0, true));

        if (t <= this->times_.back())
            return Rate(this->interpolation_(t, true) / t);

        // flat r extrapolation
        return Rate(this->data_.back() / this->times_.back());
    }

    template <class T>
    InterpolatedRateTimeCurve<T>::InterpolatedRateTimeCurve(const DayCounter& dayCounter,
                                                            const T& interpolator)
    : RateTimeStructure(dayCounter), InterpolatedCurve<T>(interpolator) {}

    template <class T>
    InterpolatedRateTimeCurve<T>::InterpolatedRateTimeCurve(const Date& referenceDate,
                                                            const DayCounter& dayCounter,
                                                            const std::vector<Handle<Quote>>& jumps,
                                                            const std::vector<Date>& jumpDates,
                                                            const T& interpolator)
    : RateTimeStructure(referenceDate, Calendar(), dayCounter, jumps, jumpDates),
      InterpolatedCurve<T>(interpolator) {}

    template <class T>
    InterpolatedRateTimeCurve<T>::InterpolatedRateTimeCurve(Natural settlementDays,
                                                            const Calendar& calendar,
                                                            const DayCounter& dayCounter,
                                                            const std::vector<Handle<Quote>>& jumps,
                                                            const std::vector<Date>& jumpDates,
                                                            const T& interpolator)
    : RateTimeStructure(settlementDays, calendar, dayCounter, jumps, jumpDates),
      InterpolatedCurve<T>(interpolator) {}

    template <class T>
    InterpolatedRateTimeCurve<T>::InterpolatedRateTimeCurve(const std::vector<Date>& dates,
                                                            const std::vector<Rate>& rate_times,
                                                            const DayCounter& dayCounter,
                                                            const Calendar& calendar,
                                                            const std::vector<Handle<Quote>>& jumps,
                                                            const std::vector<Date>& jumpDates,
                                                            const T& interpolator)
    : RateTimeStructure(dates.at(0), calendar, dayCounter, jumps, jumpDates),
      InterpolatedCurve<T>(std::vector<Time>(), rate_times, interpolator), dates_(dates) {
        initialize();
    }

    template <class T>
    InterpolatedRateTimeCurve<T>::InterpolatedRateTimeCurve(const std::vector<Date>& dates,
                                                            const std::vector<Rate>& rate_times,
                                                            const DayCounter& dayCounter,
                                                            const Calendar& calendar,
                                                            const T& interpolator)
    : RateTimeStructure(dates.at(0), calendar, dayCounter),
      InterpolatedCurve<T>(std::vector<Time>(), rate_times, interpolator), dates_(dates) {
        initialize();
    }

    template <class T>
    InterpolatedRateTimeCurve<T>::InterpolatedRateTimeCurve(const std::vector<Date>& dates,
                                                            const std::vector<Rate>& rate_times,
                                                            const DayCounter& dayCounter,
                                                            const T& interpolator)
    : RateTimeStructure(dates.at(0), Calendar(), dayCounter),
      InterpolatedCurve<T>(std::vector<Time>(), rate_times, interpolator), dates_(dates) {
        initialize();
    }

#endif

    template <class T>
    void InterpolatedRateTimeCurve<T>::initialize() {
        QL_REQUIRE(dates_.size() >= T::requiredPoints, "not enough input dates given");
        QL_REQUIRE(this->data_.size() == dates_.size(), "dates/data count mismatch");

        this->setupTimes(dates_, dates_[0], dayCounter());
        this->setupInterpolation();
        this->interpolation_.update();
    }

}

#endif
