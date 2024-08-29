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

/*! \file splinestructure.hpp
    \brief B-spline interpolation
*/

#ifndef spline_structure_hpp
#define spline_structure_hpp

#include "bsplineevaluator.hpp"
#include "splineconstraints.hpp"
#include <vector>
#include <utility>
#include <ql/types.hpp>
#include <ql/shared_ptr.hpp>

//using EigenMatrix = std::variant<Eigen::MatrixXd, Eigen::SparseMatrix<double>>;

namespace QuantLib {
    // Enums for Side and PrimeInterpolationType

    class BSplineStructure {
      public:
        enum class Side { Left, Right, Average, Actual, Inside };

        enum class InterpolationType {
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

        enum class InterpolationSmoothness {
            Discontinous, // Internal knots repeated k times for k-th degree spline, aka C^{-1}
            Continous,    // aka C^0
            ContinuouslyDifferentiable,      // aka C^1
            TwiceContinuouslyDifferentiable, // aka C^2
            Hermite, // Internal knots are double, means C^{k-2} for a k-th degree spline
            Default  // Internal knots are simple, means C^{k-1} for a k-th degree spline
        };

        BSplineStructure(const std::vector<double>& simpleKnots,
                         Size degree,
                         const std::vector<Integer>& knotIndices = {},
                         ext::shared_ptr<SplineConstraints>& splineConstraints =
                             ext::make_shared<SplineConstraints>(),
                         InterpolationSmoothness smoothness = InterpolationSmoothness::Default,
                         Side side = Side::Right,
                         // double knotTolerance = 1e-8,
                         // double rankTolerance = 1e-10,
                         int requiredPoints = 1,
                         bool isGlobal = true);

        BSplineStructure(const BSplineStructure& other);

        const std::pair<double, double> range() const;
        const std::vector<double>& knots() const;

        Size degree() const;

        void setConstraints(ext::shared_ptr<SplineConstraints>& splineConstraints);
        void addInterpolationNodes(const std::vector<double>& interpolationNodes);

        Eigen::VectorXd interpolate(const std::vector<double>& interpolationNodes,
                                    const std::vector<double> values);

        Eigen::VectorXd solve(const std::vector<double> parameters);

        Eigen::VectorXd getSolution();

        Eigen::VectorXd
        evaluateAll(double x, Size degree = -1, Side side = Side::Right) const;

        // Method for the derivative matrix
        Eigen::SparseMatrix<double> singleDerivativeMatrix(Size degree = -1,
                                                           bool differenceOperator = false) const;

        /**
         * @brief Constructs and returns a matrix representing the antiderivative operation for a
         * spline of the given degree, in the spline basis of the knot vector, with condition x0.
         *
         * @param degree The degree of the splines for which this matrix will apply.
         * @param t0 The value at the left end point, or at any given point, of the antiderivative
         * spline.
         * @param differenceOperator If true, returns the difference operators rather than the
         * derivative.
         * @return Eigen::MatrixXd A matrix that, when multiplied by a vector of spline
         * coefficients, yields the vector of spline coefficients for the antiderivative of the
         * spline.
         */
        Eigen::MatrixXd
        singleAntiDerivativeMatrix(Size degree = -1,
                                   double t0 = std::numeric_limits<double>::quiet_NaN(),
                                   bool differenceOperator = false) const;

        /**
         * @brief Computes the derivative matrix of order nu for spline coefficients of splines for
         * a knot vector.
         *
         * @param nu The order of the derivative, negative indicates anti-derivative.
         * @param degree The degree of the splines in the basis to compute the derivative matrix
         * for.
         * @param differenceOperator If true, returns the difference operator for coefficients.
         * @param t0 Reference point for the constant of integration.
         * @return Eigen::SparseMatrix<double> or Eigen::MatrixXd. If nu >= 0 a sparse matrix,
         * otherwise dense, that when multiplied by a vector of spline coefficients, yields the
         * vector of spline coefficients for the derivative of order nu for the spline.
         */
        // EigenMatrix derivativeMatrix(int nu = 1,
        //                  Size degree = -1,
        //                  bool differenceOperator = false,
        //                  double t0 = std::numeric_limits<double>::quiet_NaN()) const;


      private:
        Size degree_;
        std::vector<double> simpleKnots_;
        Size nSimpleKnots_;
        Size nKnots_;
        std::vector<Integer> knotIndices_;
        // std::vector<Size> structure_;
        InterpolationSmoothness interpolationSmoothness_;
        std::vector<double> knots_;
        std::vector<double> interpolationNodes_;

        BSplineEvaluator spline_;
        ext::shared_ptr<SplineConstraints> splineConstraints_;

        double startPoint_;
        double endPoint_;
        Side side_;
        int requiredPoints_;
        bool isGlobal_;
        // double knotTolerance_;
        // double rankTolerance_;

        void defaultKnotIndices();

        /**
         * Calculates the derivative of the spline at a given point.
         *
         * @param x The point at which to evaluate the derivative.
         * @param nu The order of the derivative.
         * @param degree The degree of the spline.
         * @param x0 Reference point for constant of integration.
         * @param side The side for evaluation.
         * @return A vector representing the derivative at the point.
         */
        Eigen::VectorXd derivative_(double x,
                                    int nu = 1,
                                    Size degree = -1,
                                    double x0 = 0.0,
                                    Side side = Side::Right) const;
    };

    static bool checkKnotIndices(const std::vector<Integer>& knotIndices,
                                 Integer degree,
                                 Size n);
}


#endif // spline_structure_hpp