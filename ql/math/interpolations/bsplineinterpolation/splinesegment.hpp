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

/*! \file splinesegment.hpp
    \brief B-spline interpolation segment. A segment has a structure for a B-spline i.e. nodes and
   degree, can have a transformation and a side, but not constraints and objective. It can therefore
   calculate a spline value given coefficients, but not solve for nor interpolate data.
*/

#ifndef spline_segment_hpp
#define spline_segment_hpp

#include "bsplineevaluator.hpp"
#include <ql/shared_ptr.hpp>
#include <ql/types.hpp>
#include <utility>
#include <vector>

/*!
* \brief Check the validity of knot indices.
* \param knotIndices Vector of knot indices.
* \param degree Degree of the spline.
* \param n Number of knots.
* \return Boolean indicating if the knot indices are valid.
*/
bool checkKnotIndices(const std::vector<Integer>& knotIndices, Integer degree, Integer n);

namespace QuantLib {

    /*!
     * \brief Class representing the structure of a B-spline.
     */
    class BSplineSegment {
      public:
        /*!
         * \brief Enum representing the side for interpolation, the 
         */
        enum class Side : std::int8_t { Left, Right, Average, Actual, Inside };

        /*!
         * \brief Enum representing the type of interpolation.
         */
        enum class InterpolationType : std::int8_t {
                Linear,
                LinearRT,
                QuadraticSplineRT,
                CubicSpline,
                CubicSplineRT,
                Hermite,
                HermiteRT,
                PiecewiseQuadraticContFwd,
                ConvexMonotone,
                Shortest,
                Closest,
                Following,
                Preceding
        };

        /*!
         * \brief Enum representing the transform for interpolation.
         */
        enum class InterpolationTransform : std::int8_t { Default, Log, Exp, RateTime };

        /*!
         * \brief Enum representing the smoothness of interpolation.
         */
        enum class InterpolationSmoothness :std::int8_t {
                Discontinuous, // Internal knots repeated k times for k-th degree spline, aka C^{-1}
                Continuous,    // aka C^0
                ContinuouslyDifferentiable,      // aka C^1
                TwiceContinuouslyDifferentiable, // aka C^2
                Hermite, // Internal knots are double, means C^{k-2} for a k-th degree spline
                Default  // Internal knots are simple, means C^{k-1} for a k-th degree spline
        };

        /*!
         * \brief Constructor for BSplineSegment.
         * \param simpleKnots Vector of simple knots.
         * \param degree Degree of the spline.
         * \param knotIndices Vector of knot indices.
         * \param smoothness Smoothness of the interpolation.
         * \param interpolationTransform Transform for the interpolation.
         * \param side Side for the interpolation.
         * \param requiredPoints Number of required points.
         * \param isGlobal Boolean indicating if the spline is global.
         */
        BSplineSegment(const std::vector<Real>& simpleKnots,
                         Integer degree,
            const std::vector<Integer>& knotIndices,
                         InterpolationSmoothness smoothness = InterpolationSmoothness::Default,
                         InterpolationTransform interpolationTransform = InterpolationTransform::Default,
                         Side side = Side::Right,
                         Size requiredPoints = 1,
                         bool isGlobal = true);

        /*!
         * \brief Default constructor for BSplineSegment.
         */
        BSplineSegment() = default;

        /*!
         * \brief Copy constructor for BSplineSegment.
         * \param other The other BSplineSegment to copy.
         */
        BSplineSegment(const BSplineSegment& other);

        /*!
         * \brief Assignment operator for BSplineSegment.
         * \param other The other BSplineSegment to assign.
         * \return A reference to this BSplineSegment.
         */
        BSplineSegment& operator=(const BSplineSegment& other);

        ~BSplineSegment() = default;

        /*!
         * \brief Move constructor for BSplineSegment.
         * \param other The other BSplineSegment to move.
         */
        BSplineSegment(BSplineSegment&& other) noexcept = default;

        /*!
         * \brief Move assignment operator for BSplineSegment.
         * \param other The other BSplineSegment to move.
         * \return A reference to this BSplineSegment.
         */
        BSplineSegment& operator=(BSplineSegment&& other) noexcept = default;

        /*!
         * \brief Get the range of the spline.
         * \return A pair representing the range of the spline.
         */
        std::pair<Real, Real> range() const;

        /*!
         * \brief Get the knots of the spline.
         * \return A vector of knots.
         */
        const std::vector<Real>& knots() const;

        /*!
         * \brief Get the degree of the spline.
         * \return The degree of the spline.
         */
        Size degree() const;

        /*!
         * \brief Get the number of variables for the spline segment, i.e. the number of spline basis functions.
        * \return The number of variables.
        */
        Size getNumVariables() const { return nKnots_ - degree_ - static_cast<Size>(1); }

        /*!
         * \brief Evaluate all basis functions at a given value.
         * \param x The value at which to evaluate the basis functions.
         * \param degree The degree of the spline.
         * \param side The side for evaluation.
         * \return A vector of evaluated basis functions.
         */
        Eigen::VectorXd
        evaluateAll(Real x, Size degree = static_cast<Size>(-1), Side side = Side::Right) const;

        /*!
         * \brief Construct and return a matrix representing the single derivative operation.
         * \param degree The degree of the spline.
         * \param differenceOperator Boolean indicating if the difference operator should be used.
         * \return A sparse matrix representing the single derivative operation.
         */
        Eigen::SparseMatrix<Real> singleDerivativeMatrix(Size degree = static_cast<Size>(-1),
                                                           bool differenceOperator = false) const;

        /*!
         * \brief Construct and return a matrix representing the antiderivative operation.
         * \param degree The degree of the spline.
         * \param t0 The value at the left end point or any given point of the antiderivative
         * spline.
         * \param differenceOperator Boolean indicating if the difference operator should be
         * used.
         * \return A matrix representing the antiderivative operation.
         */
        Eigen::MatrixXd singleAntiDerivativeMatrix(Size degree = static_cast<Size>(-1),
                                   Real t0 = std::numeric_limits<Real>::quiet_NaN(),
                                   bool differenceOperator = false) const;

        // Transformation functions
        std::function<Real(Real, Real)> transform;
        std::function<Real(Real, Real, Natural)> transformDerivative;
        std::function<Real(Real, Real)> inverse;
        Real x0_ = 0.0;

      private:
        std::vector<Real> simpleKnots_;
        Integer degree_;
        std::vector<Integer> knotIndices_;

        InterpolationSmoothness interpolationSmoothness_;
        InterpolationTransform interpolationTransform_;
        Side side_;

        Size requiredPoints_;
        bool isGlobal_;

        Real startPoint_;
        Real endPoint_;

        std::vector<Real> knots_;
        Integer nSimpleKnots_;
        Integer nKnots_;

        BSplineEvaluator spline_;

        /*!
         * \brief Set default knot indices.
         */
        void defaultKnotIndices();

    public:
        /*!
         * \brief Calculate the derivative of the spline at a given point.
         * \param x The point at which to evaluate the derivative.
         * \param nu The order of the derivative.
         * \param degree The degree of the spline.
         * \param x0 Reference point for constant of integration.
         * \param side The side for evaluation.
         * \return A vector representing the derivative at the point.
         */
        Eigen::VectorXd derivative(Real x,
                                    Integer nu = 1,
                                    Size degree = static_cast<Size>(-1),
                                    Real x0 = 0.0,
                                    Side side = Side::Right) const;

        Real value(const Eigen::VectorXd& coefficients, Real t, Side side = Side::Right);
    };

}

#endif // spline_segment_hpp
