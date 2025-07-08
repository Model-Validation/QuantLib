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


namespace QuantLib {

    /*!
     * \brief Class representing the structure of a B-spline.
     */
    class BSplineSegment {
      public:
        /*!
         * \brief Enum representing the side for interpolation, the
         */
        enum class Side : std::int8_t { Left, Right, Average, Actual, Inside, None };

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
        enum class InterpolationTransform : std::int8_t {
            Default,
            Log,
            Exp,
            RateTime,
            RateTimeAnnualToContinuous,
            ContinuousToAnnual,
            ContinuousToSimple
        };

        /*!
         * \brief Enum representing the smoothness of interpolation.
         */
        enum class InterpolationSmoothness : std::int8_t {
            Discontinuous, // Internal knots repeated k times for k-th degree spline, aka C^{-1}
            Continuous,    // aka C^0
            ContinuouslyDifferentiable,      // aka C^1
            TwiceContinuouslyDifferentiable, // aka C^2
            Hermite, // Internal knots are double, means C^{k-2} for a k-th degree spline
            Default  // Internal knots are simple, means C^{k-1} for a k-th degree spline
        };


        // Convert Side enum to string
        static constexpr std::string_view toString(Side side) {
            switch (side) {
                case Side::Left:
                    return "Left";
                case Side::Right:
                    return "Right";
                case Side::Average:
                    return "Average";
                case Side::Actual:
                    return "Actual";
                case Side::Inside:
                    return "Inside";
                case Side::None:
                    return "None";
                //default:
                //    return "Unknown Side";
            }
            return "Unknown Side";
        }

        // Convert InterpolationType enum to string
        static constexpr std::string_view toString(InterpolationType type) {
            switch (type) {
                case InterpolationType::Linear:
                    return "Linear";
                case InterpolationType::LinearRT:
                    return "LinearRT";
                case InterpolationType::QuadraticSplineRT:
                    return "QuadraticSplineRT";
                case InterpolationType::CubicSpline:
                    return "CubicSpline";
                case InterpolationType::CubicSplineRT:
                    return "CubicSplineRT";
                case InterpolationType::Hermite:
                    return "Hermite";
                case InterpolationType::HermiteRT:
                    return "HermiteRT";
                case InterpolationType::PiecewiseQuadraticContFwd:
                    return "PiecewiseQuadraticContFwd";
                case InterpolationType::ConvexMonotone:
                    return "ConvexMonotone";
                case InterpolationType::Shortest:
                    return "Shortest";
                case InterpolationType::Closest:
                    return "Closest";
                case InterpolationType::Following:
                    return "Following";
                case InterpolationType::Preceding:
                    return "Preceding";
                //default:
                //    return "Unknown InterpolationType";
            }
            return "Unknown InterpolationType";
        }

        // Convert InterpolationTransform enum to string
        static constexpr std::string_view toString(InterpolationTransform transform) {
            switch (transform) {
                case InterpolationTransform::Default:
                    return "Default";
                case InterpolationTransform::Log:
                    return "Log";
                case InterpolationTransform::Exp:
                    return "Exp";
                case InterpolationTransform::RateTime:
                    return "RateTime";
                case InterpolationTransform::RateTimeAnnualToContinuous:
                    return "RateTimeAnnualToContinuous";
                case InterpolationTransform::ContinuousToAnnual:
                    return "ContinuousToAnnual";
                case InterpolationTransform::ContinuousToSimple:
                    return "ContinuousToSimple";
                //default:
                //    return "Unknown InterpolationTransform";
            }
            return "Unknown InterpolationTransform";
        }

        // Convert InterpolationSmoothness enum to string
        static constexpr std::string_view toString(InterpolationSmoothness smoothness) {
            switch (smoothness) {
                case InterpolationSmoothness::Discontinuous:
                    return "Discontinuous";
                case InterpolationSmoothness::Continuous:
                    return "Continuous";
                case InterpolationSmoothness::ContinuouslyDifferentiable:
                    return "ContinuouslyDifferentiable";
                case InterpolationSmoothness::TwiceContinuouslyDifferentiable:
                    return "TwiceContinuouslyDifferentiable";
                case InterpolationSmoothness::Hermite:
                    return "Hermite";
                case InterpolationSmoothness::Default:
                    return "Default";
                //default:
                //    return "Unknown InterpolationSmoothness";
            }
            return "Unknown InterpolationSmoothness";
        }


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
        BSplineSegment(
            const std::vector<Real>& simpleKnots,
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


        // Accessor for sidedness
        Side side() const { return side_; }

        // Accessor for interpolationTransform
        InterpolationTransform interpolationTransform() const {
            return interpolationTransform_;
        }

        // Accessor for interpolationSmoothness
        InterpolationSmoothness interpolationSmoothness() const {
            return interpolationSmoothness_;
        }

        // Accessor for sidedness for SWIG
        std::string side_str() const { return std::string(toString(side_)); }

        // Accessor for interpolationTransform for SWIG
        std::string interpolation_transform() const {
            return std::string(toString(interpolationTransform_));
        }

        // Accessor for interpolationSmoothness for SWIG
        std::string interpolation_smoothness() const {
            return std::string(toString(interpolationSmoothness_));
        }


        /*!
         * \brief Get the number of variables for the spline segment, i.e. the number of spline
         * basis functions. \return The number of variables.
         */
        Size getNumVariables() const { return nKnots_ - degree_ - static_cast<Size>(1); }
        // ReSharper disable once CppInconsistentNaming
        Natural get_num_variables() const { return static_cast<Natural>(nKnots_ - degree_) - 1; }

        /*!
         * \brief Evaluate all basis functions at a given value.
         * \param t The value at which to evaluate the basis functions.
         * \param degree The degree of the spline.
         * \param side The side for evaluation.
         * \return A vector of evaluated basis functions.
         */
        Eigen::VectorXd
        evaluateAll(Real t, Size degree = static_cast<Size>(-1), Side side = Side::Right) const;

        // ReSharper disable once CppInconsistentNaming
        std::vector<Real>
        evaluate_all(Real x, Integer degree = -1, Side side = Side::Right) const {
            Eigen::VectorXd eigenVector = evaluateAll(x, static_cast<Size>(degree), side);
            return {eigenVector.data(), eigenVector.data() + eigenVector.size()};
        }

        /*!
         * \brief Construct and return a matrix representing the single derivative operation.
         * \param degree The degree of the spline.
         * \param differenceOperator Boolean indicating if the difference operator should be used.
         * \return A sparse matrix representing the single derivative operation.
         */
        Eigen::SparseMatrix<Real> singleDerivativeMatrix(Size degree = static_cast<Size>(-1),
                                                         bool differenceOperator = false) const;

        // ReSharper disable once CppInconsistentNaming
        std::vector<std::vector<Real>>
        single_derivative_matrix(Integer degree = -1,
                                 bool differenceOperator = false) const {
            // Call the first method to get the sparse matrix
            Eigen::SparseMatrix<Real> sparse_matrix =
                singleDerivativeMatrix(static_cast<Size>(degree), differenceOperator);

            // Initialize the output vector
            std::vector<std::vector<Real>> result(sparse_matrix.rows(),
                                                  std::vector<Real>(sparse_matrix.cols(), 0.0));

            // Convert the sparse matrix to a dense std::vector<std::vector<Real>>
            for (int k = 0; k < sparse_matrix.outerSize(); ++k) {
                for (Eigen::SparseMatrix<Real>::InnerIterator it(sparse_matrix, k); it;
                     ++it) {
                    result[it.row()][it.col()] = it.value();
                }
            }

            return result;
        }

        /*!
         * \brief Construct and return a matrix representing the antiderivative operation.
         * \param degree The degree of the spline.
         * \param t0 The value at the left end point or any given point of the antiderivative
         * spline.
         * \param differenceOperator Boolean indicating if the difference operator should be
         * used.
         * \return A matrix representing the antiderivative operation.
         */
        Eigen::SparseMatrix<Real>
        singleAntiDerivativeMatrix(Size degree = static_cast<Size>(-1),
                                                   Real t0 = std::numeric_limits<Real>::quiet_NaN(),
                                                   bool differenceOperator = false) const;

        // ReSharper disable once CppInconsistentNaming
        std::vector<std::vector<Real>>
        single_anti_derivative_matrix(Integer degree = -1,
                                      Real t0 = std::numeric_limits<Real>::quiet_NaN(),
                                      bool differenceOperator = false) const {
            // Call the first method to get the dense Eigen::MatrixXd
            Eigen::SparseMatrix<Real> sparse_matrix =
                singleAntiDerivativeMatrix(static_cast<Size>(degree), t0, differenceOperator);

            // Initialize the output vector of vectors
            std::vector<std::vector<Real>> result(sparse_matrix.rows(),
                                                  std::vector<Real>(sparse_matrix.cols()));

            // Convert the sparse matrix to a dense std::vector<std::vector<Real>>
            for (int k = 0; k < sparse_matrix.outerSize(); ++k) {
                for (Eigen::SparseMatrix<Real>::InnerIterator it(sparse_matrix, k); it; ++it) {
                    result[it.row()][it.col()] = it.value();
                }
            }

            return result;
        }

        Eigen::SparseMatrix<Real> derivativeMatrix(Integer nu = 1,
                                                   Size degree = static_cast<Size>(-1),
                                                   bool differenceOperator = false) const;

        // ReSharper disable once CppInconsistentNaming
        std::vector<std::vector<Real>>
        derivative_matrix(Integer nu = 1,
                          Integer degree = -1,
                          bool differenceOperator = false) const {
            // Call the first method to get the sparse matrix
            Eigen::SparseMatrix<Real> sparse_matrix =
                derivativeMatrix(nu, static_cast<Size>(degree), differenceOperator);

            // Initialize the output vector
            std::vector<std::vector<Real>> result(sparse_matrix.rows(),
                                                  std::vector<Real>(sparse_matrix.cols(), 0.0));

            // Convert the sparse matrix to a dense std::vector<std::vector<Real>>
            for (int k = 0; k < sparse_matrix.outerSize(); ++k) {
                for (Eigen::SparseMatrix<Real>::InnerIterator it(sparse_matrix, k); it; ++it) {
                    result[it.row()][it.col()] = it.value();
                }
            }

            return result;
        }


        Eigen::SparseMatrix<Real> antiDerivativeMatrix(Integer nu = -1,
                                             Size degree = static_cast<Size>(-1),
                                             const Eigen::VectorXd& t0 = {},
                                             bool differenceOperator = false) const;


        Eigen::SparseMatrix<Real>
        antiDerivativeMatrix(const Integer nu = -1,
                                             const Size degree = static_cast<Size>(-1),
                                             const Real t0 = std::numeric_limits<Real>::quiet_NaN(),
                                             const bool differenceOperator = false) const {
            Eigen::VectorXd t0Vector(1);
            t0Vector << t0;
            return antiDerivativeMatrix(nu, degree, t0Vector, differenceOperator);
        }

        // ReSharper disable once CppInconsistentNaming
        std::vector<std::vector<Real>>
        anti_derivative_matrix(Integer nu = -1,
                               Integer degree = -1,
                               const std::vector<Real>& t0 = {},
                               bool differenceOperator = false) const {
            // Call the first method to get the dense Eigen::MatrixXd
            Eigen::SparseMatrix<Real> sparse_matrix = antiDerivativeMatrix(
                nu, static_cast<Size>(degree),
                Eigen::Map<const Eigen::VectorXd>(t0.data(), static_cast<Eigen::Index>(t0.size())),
                differenceOperator);

            // Initialize the output vector of vectors
            // Initialize the output vector
            std::vector<std::vector<Real>> result(sparse_matrix.rows(),
                                                  std::vector<Real>(sparse_matrix.cols(), 0.0));

            // Convert the sparse matrix to a dense std::vector<std::vector<Real>>
            for (int k = 0; k < sparse_matrix.outerSize(); ++k) {
                for (Eigen::SparseMatrix<Real>::InnerIterator it(sparse_matrix, k); it; ++it) {
                    result[it.row()][it.col()] = it.value();
                }
            }

            return result;
        }

        // ReSharper disable once CppInconsistentNaming
        std::vector<Real> get_simple_knots() const { return simpleKnots_; }
        // ReSharper disable once CppInconsistentNaming
        std::vector<Integer> get_knot_indices() const { return knotIndices_; }

        // Transformation functions
        std::function<Real(Real, Real)> transform;
        std::function<Real(Real, Real, Natural)> transformDerivative;
        std::function<Real(Real, Real)> inverse;
        Real rateTimeX0 = 0.0;

      private:
        std::vector<Real> simpleKnots_;
        Size degree_;
        std::vector<Integer> knotIndices_;

        InterpolationSmoothness interpolationSmoothness_;
        InterpolationTransform interpolationTransform_;
        Side side_;

        Size requiredPoints_;
        bool isGlobal_;

        Real startPoint_ = 0.0;
        Real endPoint_;

        std::vector<Real> knots_;
        Integer nSimpleKnots_;
        Integer nKnots_;

        BSplineEvaluator spline_;
        std::unordered_map<Integer, ext::shared_ptr<Eigen::SparseMatrix<Real>>> derivativeCache_ = {};
        std::unordered_map<Integer, ext::shared_ptr<Eigen::SparseLU<Eigen::SparseMatrix<Real>>>> antiDerivativeCache_ = {};
        std::unordered_map<Integer, ext::shared_ptr<BSplineEvaluator>> splineCache_ = {};

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
        Eigen::VectorXd derivativeFunctional(Real x,
                                   Integer nu = 1,
                                   Size degree = static_cast<Size>(-1),
                                   Real x0 = 0.0,
                                   Side side = Side::None) const;

        // ReSharper disable once CppInconsistentNaming
        std::vector<Real> derivative_functional(const Real x,
                                                const Integer nu = 1,
                                                const Size degree = static_cast<Size>(-1),
                                                const Real x0 = 0.0,
                                                const Side side = Side::None) const {
            Eigen::VectorXd eigenVector = derivativeFunctional(x, nu, degree, x0, side);
            return {eigenVector.data(), eigenVector.data() + eigenVector.size()};
        }

        Real value(const Eigen::VectorXd& coefficients, Real t, Integer nu = 0, Side side = Side::Right);

        Real value(const std::vector<Real>& coefficients, Real t, Integer nu = 0, Side side = Side::Right) {
            return value(Eigen::Map<const Eigen::VectorXd>(
                             coefficients.data(), static_cast<Eigen::Index>(coefficients.size())),
                         t, nu, side);
        }

        Eigen::VectorXd valueFunctional(Real t, Side side) const;

        // ReSharper disable once CppInconsistentNaming
        std::vector<Real> value_functional(const Real t, const Side side) const {
            Eigen::VectorXd eigenVector = valueFunctional(t, side);
            return {eigenVector.data(), eigenVector.data() + eigenVector.size()};
        }

    };

    static bool checkKnotIndices(const std::vector<Integer>& knotIndices, Size degree, Integer n);
    static std::string vectorToString(const std::vector<QuantLib::Integer>& matrix);
}

#endif // spline_segment_hpp
