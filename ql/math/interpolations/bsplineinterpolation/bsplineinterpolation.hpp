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

#ifndef b_spline_interpolation_hpp
#define b_spline_interpolation_hpp

#include "bsplinestructure.hpp"
#include <ql/shared_ptr.hpp>
#include <ql/types.hpp>
#include <ql/math/interpolation.hpp>
#include <ql/math/interpolations/bsplineinterpolation/splinesegment.hpp>
#include <vector>

/*! \file bsplineinterpolation.hpp
    \brief B-spline interpolation based on linear equality and inequality constraints on
   coefficients, and a quadratic objective. This allows separation spline knots and interpolation
   nodes, under-determined constraints
*/

namespace QuantLib {
    namespace detail {

        /*!
         * \brief Implementation of B-spline interpolation.
         * \tparam I1 Iterator type for x values.
         * \tparam I2 Iterator type for y values.
         */
        template <class I1, class I2>
        class BSplineInterpolationImpl final : public Interpolation::templateImpl<I1, I2> {

          public:
            /*!
             * \brief Constructor for BSplineInterpolationImpl.
             * \param xBegin Iterator to the beginning of x values.
             * \param xEnd Iterator to the end of x values.
             * \param yBegin Iterator to the beginning of y values.
             * \param splineStructure Shared pointer to the B-spline structure.
             */
            BSplineInterpolationImpl(const I1& xBegin,
                                        const I1& xEnd,
                                        const I2& yBegin,
                                        const ext::shared_ptr<BSplineStructure>& splineStructure)
                : Interpolation::templateImpl<I1, I2>(xBegin, xEnd, yBegin),
                    xSize_(static_cast<Size>(std::distance(xBegin, xEnd))),
              splineStructure_(ext::make_shared<BSplineStructure>(*splineStructure)) {
            }

            /*!
             * \brief Update the interpolation coefficients.
             */
            void update() override {
                coefficients_ = splineStructure_->interpolate(
                    std::vector<double>(this->xBegin_, this->xEnd_),
                    std::vector<double>(this->yBegin_, this->yBegin_ + xSize_));
            }

            /*!
             * \brief Get the interpolated value at a given point.
             * \param x The point at which to get the value.
             * \return The interpolated value.
             */
            Real value(Real x) const override {
                // For B-splines, we allow natural polynomial extrapolation
                // The splineStructure_->value() method naturally extends the polynomial
                // beyond the domain bounds, providing continuous extrapolation
                return splineStructure_->value(coefficients_, x);
            }

            /*!
             * \brief Get the primitive of the interpolation.
             * \param x The point at which to get the primitive.
             * \return The primitive value.
             */
            Real primitive(Real x) const override {
                // Natural polynomial extrapolation for primitives
                return splineStructure_->value(coefficients_, x, -1);
            }

            /*!
             * \brief Get the first derivative of the interpolation.
             * \param x The point at which to get the derivative.
             * \return The first derivative value.
             */
            Real derivative(Real x) const override {
                // Natural polynomial extrapolation for derivatives
                return splineStructure_->value(coefficients_, x, 1);
            }

            /*!
             * \brief Get the second derivative of the interpolation.
             * \param x The point at which to get the second derivative.
             * \return The second derivative value.
             */
            Real secondDerivative(Real x) const override {
                // Natural polynomial extrapolation for second derivatives
                return splineStructure_->value(coefficients_, x, 2);
            }

            /*!
             * \brief Get the minimum x value of the interpolation range.
             * \return The minimum x value.
             */
            Real xMin() const override { return splineStructure_->range().first; }

            /*!
             * \brief Get the maximum x value of the interpolation range.
             * \return The maximum x value.
             */
            Real xMax() const override { return splineStructure_->range().second; }

            [[nodiscard]] std::vector<Real> getCoefficients() const override {
                return {coefficients_.data(), coefficients_.data() + coefficients_.size()};
            }

          private:
            Size xSize_; /*!< Size of the x values. */
            ext::shared_ptr<BSplineStructure>
                splineStructure_;              /*!< Shared pointer to the B-spline structure. */
            Eigen::VectorXd coefficients_;           /*!< Interpolation coefficients. */
        };

    }; // end namespace detail

    //! B-spline interpolation between discrete points
    /*!
        \ingroup interpolations
        \warning See the Interpolation class for information about the
                 required lifetime of the underlying data.
    */
    class BSplineInterpolation final : public Interpolation {
      public:
        /*! \pre the \f$ x \f$ values must be sorted.

        */
        BSplineInterpolation() = default;

        /*!
         * \brief Constructor for BSplineInterpolation with composite spline structure.
         * \tparam I1 Iterator type for x values.
         * \tparam I2 Iterator type for y values.
         * \param xBegin Iterator to the beginning of x values.
         * \param xEnd Iterator to the end of x values.
         * \param yBegin Iterator to the beginning of y values.
         * \param splineStructure Shared pointer to the composite B-spline structure.
         */
        template <class I1, class I2>
        BSplineInterpolation(const I1& xBegin,
                             const I1& xEnd,
                             const I2& yBegin,
                             const ext::shared_ptr<BSplineStructure>& splineStructure,
                             bool enableExtrapolation = true) {
            splineStructure_ = splineStructure;

            impl_ = ext::make_shared<detail::BSplineInterpolationImpl<I1, I2>>(xBegin, xEnd, yBegin,
                                                                               splineStructure);
            impl_->update();
            
            // Enable extrapolation by default for B-splines to allow natural polynomial continuation
            if (enableExtrapolation) {
                this->enableExtrapolation();
            }
        }

        /*!
         * \brief Non-template constructor for BSplineInterpolation.
         * \param x Vector of x values.
         * \param y Vector of y values.
         * \param splineStructure Shared pointer to the B-spline structure.
         */
        BSplineInterpolation(const std::vector<double>& x,
                             const std::vector<double>& y,
                             const ext::shared_ptr<BSplineStructure>& splineStructure,
                             bool enableExtrapolation = true) {
            using ConstIterator = std::vector<double>::const_iterator;

            splineStructure_ = splineStructure;
            impl_ =
                ext::make_shared<detail::BSplineInterpolationImpl<ConstIterator, ConstIterator>>(
                    x.begin(), x.end(), y.begin(), splineStructure);
            impl_->update();
            
            // Enable extrapolation by default for B-splines
            if (enableExtrapolation) {
                this->enableExtrapolation();
            }
        }

        // ReSharper disable once CppInconsistentNaming
        [[nodiscard]] ext::shared_ptr<BSplineStructure> get_structure() const {
            return splineStructure_;
        }

        // ReSharper disable once CppInconsistentNaming
        [[nodiscard]] std::vector<Real> get_coefficients() const {
            return impl_->getCoefficients();
        }

    private:
        ext::shared_ptr<BSplineStructure> splineStructure_;

    };

    //! %BSplineModel interpolation factory and traits (unfortunately BSpline is taken)
    /*! \ingroup interpolations */
    class BSplineModel {
      public:
        /*!
         * \brief Constructor for BSplineModel.
         * \param splineStructure Shared pointer to the B-spline structure.
         */
        BSplineModel(ext::shared_ptr<BSplineStructure> splineStructure = {})
        : splineStructure_(std::move(splineStructure)) {}

        /*!
         * \brief Interpolate values using the B-spline model.
         * \tparam I1 Iterator type for x values.
         * \tparam I2 Iterator type for y values.
         * \param xBegin Iterator to the beginning of x values.
         * \param xEnd Iterator to the end of x values.
         * \param yBegin Iterator to the beginning of y values.
         * \return Interpolation object.
         */
        template <class I1, class I2>
        [[nodiscard]] BSplineInterpolation interpolate(const I1& xBegin, const I1& xEnd, const I2& yBegin, bool enableExtrapolation = true) const {
            return {xBegin, xEnd, yBegin, splineStructure_, enableExtrapolation};
        }

        [[nodiscard]] ext::shared_ptr<BSplineInterpolation> interpolate(const std::vector<Real>& x, const std::vector<Real>& y, bool enableExtrapolation = true) const {
            auto interp = ext::make_shared<BSplineInterpolation>(x.begin(), x.end(), y.begin(), splineStructure_);
            // Enable extrapolation by default for B-splines
            // This allows natural polynomial continuation beyond segment bounds
            if (enableExtrapolation) {
                interp->enableExtrapolation();
            }
            return interp;
        }

        // ReSharper disable once CppVariableCanBeMadeConstexpr
        static const bool global = true;
        // ReSharper disable once CppVariableCanBeMadeConstexpr
        static const Size requiredPoints = 2;

        /*!
         * \brief Get the start point of the spline structure.
         * \return The start point.
         */
        [[nodiscard]] Real getStartPoint() const { return splineStructure_->range().first; }

        /*!
         * \brief Get the end point of the spline structure.
         * \return The end point.
         */
        [[nodiscard]] Real getEndPoint() const { return splineStructure_->range().second; }

        // ReSharper disable once CppInconsistentNaming
        [[nodiscard]] ext::shared_ptr<BSplineStructure> get_structure() const { return splineStructure_;}

      private:
        ext::shared_ptr<BSplineStructure> splineStructure_; /*!< Shared pointer to the B-spline structure. */
    };

}

#endif // b_spline_interpolation_hpp
