/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Ferdinando Ametrano

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

/*! \file spreadedinterpolation.hpp
    \brief spreaded interpolation
*/

#ifndef quantlib_spreaded_interpolation_hpp
#define quantlib_spreaded_interpolation_hpp

#include <ql/termstructures/yield/bootstraptraits.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/daycounters/business252.hpp>
#include <ql/math/interpolation.hpp>
#include <ql/time/date.hpp>
#include <optional>

namespace QuantLib {
    template <class T>
    Date lowerDate(const Real t, const Date& refDate) {
        // Implementation for generic case
        // ...

        // Specialization for Actual360
        if constexpr (std::is_same_v<T, Actual360>) {
            return refDate + static_cast<Date::serial_type>(std::round(t * 360.0));
        }

        // Specialization for Actual365Fixed
        if constexpr (std::is_same_v<T, Actual365Fixed>) {
            return refDate + static_cast<Date::serial_type>(std::round(t * 365.0));
        }

        // Default implementation
        QL_FAIL("lowerDate not implemented\n");
    }

    namespace detail {

        template <class Traits, class I1, class I2, class Interpolator>
        class SpreadedInterpolationImpl;

    };

    //! spreaded interpolation between discrete points
    /*! \ingroup interpolations
        \warning See the Interpolation class for information about the
                 required lifetime of the underlying data.
    */
    template <class Traits>
    class SpreadedInterpolation : public Interpolation {
      public:
        /*! \pre the \f$ x \f$ values must be sorted. */
        SpreadedInterpolation() = default;

        template <class I1, class I2, class Interpolator>
        SpreadedInterpolation(const I1& xBegin,
                              const I1& xEnd,
                              const I2& yBegin,
                              Interpolator factory,
                              ext::shared_ptr<YieldTermStructure> baseCurve,
                              const std::optional<ext::shared_ptr<InterestRateIndex>>& index = std::nullopt) {
            impl_ = ext::shared_ptr<Interpolation::Impl>(
                new detail::SpreadedInterpolationImpl<Traits, I1, I2, Interpolator>(
                    xBegin, xEnd, yBegin, factory, baseCurve, index
                )
            );
            impl_->update();
        }
    };

    namespace detail {
        template <class Traits, class I1, class I2, class Interpolator>
        class SpreadedInterpolationImpl final : public Interpolation::templateImpl<I1, I2> {
          public:
            SpreadedInterpolationImpl(
                const I1& xBegin,
                const I1& xEnd,
                const I2& yBegin,
                const Interpolator& factory = Interpolator(),
                const ext::shared_ptr<YieldTermStructure>& baseCurve = ext::make_shared<YieldTermStructure>(),
                const std::optional<ext::shared_ptr<InterestRateIndex>>& index = std::nullopt)
            : Interpolation::templateImpl<I1, I2>(xBegin, xEnd, yBegin),
              factory_(std::move(factory)), baseCurve_(std::move(baseCurve)),
              spreadY_(xEnd - xBegin), index_(index) {}

            void update() override {
                // Calculate the spreaded values
                for (size_t i = 0; i < spreadY_.size(); ++i) {
                    Real xValue = this->xBegin_[i];
                    Real yValue = this->yBegin_[i];
                    Rate baseValue = Traits::rate(baseCurve_, xValue, index_);
                    spreadY_[i] = yValue - baseValue;
                }

                spread_ = factory_.interpolate(this->xBegin_, this->xEnd_, spreadY_.begin());
            }
            Real value(Real x) const override {
                return spread_(x, true) + Traits::rate(baseCurve_, x, index_);
            }
            Real primitive(Real x) const override {
                QL_FAIL("Spreaded spline primitive not implemented");
            }
            Real derivative(Real x) const override {
                QL_FAIL("Spreaded spline derivative not implemented");
            }
            Real secondDerivative(Real x) const override {
                QL_FAIL("Spreaded spline second derivative not implemented");
            }
            Real xMin() const override { return std::max(spread_.xMin(), 0.0); }
            Real xMax() const override { return std::min(spread_.xMax(), baseCurve_->maxTime()); }


          private:
            Interpolator factory_;
            Interpolation spread_;
            ext::shared_ptr<YieldTermStructure> baseCurve_;

            std::vector<Real> spreadY_;
            std::optional<ext::shared_ptr<InterestRateIndex>> index_;
        };
    };

    //! %SpreadedInterpolationModel interpolation factory and traits
    /*! \ingroup interpolations */
    template <class Traits, class Interpolator>
    class SpreadedInterpolationModel {
      public:
        SpreadedInterpolationModel(
            Interpolator factory,
            const ext::shared_ptr<YieldTermStructure>& baseCurve,
            const std::optional<ext::shared_ptr<InterestRateIndex>>& index = std::nullopt)
        : factory_(std::move(factory)), baseCurve_(baseCurve), index_(std::move(index)) {}

        template <class I1, class I2>
        SpreadedInterpolation<Traits>
        interpolate(const I1& xBegin, const I1& xEnd, const I2& yBegin) const {
            return {xBegin, xEnd, yBegin, factory_, baseCurve_, index_};
        }

        [[nodiscard]] ext::shared_ptr< SpreadedInterpolation<Traits> >
        interpolate(const std::vector<Real>& x, const std::vector<Real>& y) const {
            return ext::make_shared< SpreadedInterpolation<Traits> >(x.begin(), x.end(), y.begin(),
                                                                   factory_, baseCurve_, index_);
        }

        static constexpr bool global = true;
        static constexpr Size requiredPoints = 2;

      protected:
        Interpolator factory_;
        ext::shared_ptr<YieldTermStructure> baseCurve_;
        std::optional<ext::shared_ptr<InterestRateIndex>> index_;
    };
}

#endif
