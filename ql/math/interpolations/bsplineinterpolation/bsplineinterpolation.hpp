/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2024 SEB AB STh

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


#ifndef quantlib_b_spline_interpolation_hpp
#define quantlib_b_spline_interpolation_hpp

#include "bsplineevaluator.hpp"
#include <vector>
#include <ql/errors.hpp>
#include <ql/math/interpolation.hpp>
#include <ql/math/interpolations/bsplineinterpolation/splinestructure.hpp>
#include <ql/shared_ptr.hpp>
#include <ql/types.hpp>

/*! \file bsplineinterpolation.hpp
    \brief B-spline interpolation based on linear equality and inequality constraints on coefficients,
    and a quadratic objective. This allows separation spline knots and interpolation nodes, underdetermined constraints 
*/

namespace QuantLib {
    namespace detail {
        template <class I1, class I2>
        class BSplineInterpolationImpl : public Interpolation::templateImpl<I1, I2> {
        
            public:
                BSplineInterpolationImpl(const I1& xBegin,
                                         const I1& xEnd,
                                         const I2& yBegin,
                                       const ext::shared_ptr<BSplineStructure>& splineStructure)
                : Interpolation::templateImpl<I1, I2>(xBegin, xEnd, yBegin),
                  xSize_(Size(xEnd - xBegin)),
                splineStructure_(ext::make_shared<BSplineStructure>(*splineStructure)) {
                    splineEvaluator_ = BSplineEvaluator(splineStructure_->knots(), splineStructure_->degree());
                }

                void update() override {
                    coeffs_ = splineStructure_->interpolate(std::vector<double>(xBegin_, xEnd_),
                        std::vector<double>(yBegin_, yBegin_+xSize_));
                }

                Real extrapolate(Real x) const { return 0.0; }

                Real value(Real x) const override {
                    const Real xMin = splineStructure_->range().first;
                    const Real xMax = splineStructure_->range().second;
                    if (x < xMin || x > xMax) {
                        return this->extrapolate(x);
                    }

                    return splineEvaluator_.value(coeffs_, x);
                }

                Real primitive(Real) const override {
                    QL_FAIL("Primitive calculation not implemented "
                            "for kernel interpolation");
                }

                Real derivative(Real x) const override {
                    const Real xMin = this->xBegin_[0];
                    const Real xMax = this->xBegin_[xSize_ - 1];
                    if (x < xMin || x > xMax) {
                        return this->extrapolate(x);
                    }
                    QL_FAIL("First derivative calculation not implemented "
                            "for kernel interpolation");
                    splineStructure_->evaluateAll(x);
                }

                Real secondDerivative(Real) const override {
                    QL_FAIL("Second derivative calculation not implemented "
                            "for kernel interpolation");
                }
                Real xMin() const override { return splineStructure_->range().first; }
                Real xMax() const override { return splineStructure_->range().second; }

             private:
                Size xSize_;
                ext::shared_ptr<BSplineStructure> splineStructure_;
                BSplineEvaluator splineEvaluator_;
                Eigen::VectorXd coeffs_;

        };

    }; // end namespace detail
    //! B-spline interpolation between discrete points
    /*!
        \ingroup interpolations
        \warning See the Interpolation class for information about the
                 required lifetime of the underlying data.
    */
    class BSplineInterpolation : public Interpolation {
      public:
        /*! \pre the \f$ x \f$ values must be sorted.

        */
        template <class I1, class I2>
        BSplineInterpolation(const I1& xBegin,
                             const I1& xEnd,
                             const I2& yBegin,
                             const ext::shared_ptr<BSplineStructure>& splineStructure) {
            impl_ = ext::make_shared<detail::BSplineInterpolationImpl<I1, I2>>(xBegin, xEnd, yBegin,
                                                                               splineStructure);
            impl_->update();
        }

        // Non-template constructor
        BSplineInterpolation(const std::vector<double>& x,
                             const std::vector<double>& y,
                             const ext::shared_ptr<BSplineStructure>& splineStructure) {
            using ConstIterator = std::vector<double>::const_iterator;

            impl_ =
                ext::make_shared<detail::BSplineInterpolationImpl<ConstIterator, ConstIterator> >(
                x.begin(), x.end(), y.begin(), splineStructure);
            impl_->update();
        }
    };
        

    //! %BSplineModel interpolation factory and traits (unfortunately BSpline is taken)
    /*! \ingroup interpolations */
    class BSplineModel {
      public:
        BSplineModel(ext::shared_ptr<BSplineStructure>& splineStructure)
        : splineStructure_(std::move(splineStructure)) {}
        template <class I1, class I2>
        Interpolation interpolate(const I1& xBegin, const I1& xEnd, const I2& yBegin) const {
            return BSplineInterpolation(xBegin, xEnd, yBegin, splineStructure_);
        }
        static const bool global = true;
        static const Size requiredPoints = 2;

      private:
        ext::shared_ptr<BSplineStructure> splineStructure_;
    };

}

#endif
