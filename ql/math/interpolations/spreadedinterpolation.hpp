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

/*! \file mixedinterpolation.hpp
    \brief mixed interpolation between discrete points
*/

#ifndef quantlib_spreaded_interpolation_hpp
#define quantlib_spreaded_interpolation_hpp

#include <ql/math/interpolation.hpp>
#include <ql/termstructures/interpolatedcurve.hpp>

namespace QuantLib {

    namespace detail {

        template <class I1, class I2, class Interpolator>
        class SpreadedInterpolationImpl;

    };

    //! spreaded interpolation between discrete points
    /*! \ingroup interpolations
        \warning See the Interpolation class for information about the
                 required lifetime of the underlying data.
    */
    class SpreadedInterpolation : public Interpolation {
      public:
        /*! \pre the \f$ x \f$ values must be sorted. */
        SpreadedInterpolation() = default;

        template <class I1, class I2, class Interpolator>
        SpreadedInterpolation(const I1& xBegin,
                              const I1& xEnd,
                              const I2& yBegin,
                              Interpolator& factory,
                              ext::shared_ptr<Interpolation>& baseCurve) {
            impl_ = ext::shared_ptr<Interpolation::Impl>(
                new detail::SpreadedInterpolationImpl<I1, I2, Interpolator>(
                    xBegin, xEnd, yBegin, factory, baseCurve
                )
            );
            impl_->update();
        }
    };

    namespace detail {
        template <class I1, class I2, class Interpolator>
        class SpreadedInterpolationImpl : public Interpolation::templateImpl<I1, I2> {
          public:
            SpreadedInterpolationImpl(
                const I1& xBegin,
                const I1& xEnd,
                const I2& yBegin,
                const Interpolator& factory = Interpolator(),
                const ext::shared_ptr<Interpolation>& baseCurve = ext::make_shared<Interpolation>())
            : Interpolation::templateImpl<I1, I2>(xBegin, xEnd, yBegin),
              factory_(std::move(factory)), baseCurve_(std::move(baseCurve)),
              spreadY_(xEnd - xBegin) {}

            void update() {
                // Calculate the spreaded values
                for (size_t i = 0; i < spreadY_.size(); ++i) {
                    Real xValue = xBegin_[i];
                    Real yValue = yBegin_[i];
                    Real baseValue = baseCurve_->operator()(xValue);
                    spreadY_[i] = yValue - baseValue;
                }

                spread_ = factory_.interpolate(xBegin_, xEnd_, spreadY_.begin());
            }
            Real value(Real x) const { return spread_(x, true) + baseCurve_->operator()(x, true); }
            Real primitive(Real x) const {
                return spread_.primitive(x, true) + baseCurve_->primitive(x, true);
            }
            Real derivative(Real x) const {
                return spread_.derivative(x, true) + baseCurve_->derivative(x, true);
            }
            Real secondDerivative(Real x) const {
                return spread_.secondDerivative(x, true) + baseCurve_->secondDerivative(x, true);
            }
            Real xMin() const override { return std::max(spread_.xMin(), baseCurve_->xMin()); }
            Real xMax() const override { return std::min(spread_.xMax(), baseCurve_->xMax()); }


          private:
            Interpolator factory_;
            Interpolation spread_;
            ext::shared_ptr<Interpolation> baseCurve_;

            std::vector<Real> spreadY_;
        };
    };

    //! %SpreadedInterpolationModel interpolation factory and traits
    /*! \ingroup interpolations */
    template <class Interpolator>
    class SpreadedInterpolationModel {
      public:
        SpreadedInterpolationModel(Interpolator& factory, ext::shared_ptr<Interpolation>& baseCurve)
        : factory_(std::move(factory)), baseCurve_(*baseCurve) {}

        template <class I1, class I2> 
        Interpolation interpolate(const I1& xBegin, const I1& xEnd, const I2& yBegin) const {
            return SpreadedInterpolation(xBegin, xEnd, yBegin, factory_, ext::make_shared<Interpolation>(baseCurve_));
        }
        static const bool global = true;
        static const Size requiredPoints = 2;

      protected:
        Interpolator factory_;
        Interpolation baseCurve_;
    };
}

#endif
