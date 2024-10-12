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

/*! \file bsplineevaluator.hpp
    \brief B-spline evaluator
*/

#ifndef b_spline_evaluator_hpp
#define b_spline_evaluator_hpp

#include <ql/types.hpp>
#include <Eigen/Sparse>
#include <vector>

namespace QuantLib {

    /*!
     * \brief Class for evaluating B-splines.
     */
    class BSplineEvaluator {
      public:
        /*!
         * \brief Default constructor for BSplineEvaluator.
         */
        BSplineEvaluator();

        /*!
         * \brief Constructor for BSplineEvaluator with specified knots and degree.
         * \param knots Vector of knots.
         * \param degree Degree of the spline.
         */
        BSplineEvaluator(const std::vector<double>& knots, Size degree);

        /*!
         * \brief Evaluate all basis functions at a given value.
         * \param x The value at which to evaluate the basis functions.
         * \return A vector of evaluated basis functions.
         */
        Eigen::VectorXd evaluateAll(double x) const;

        /*!
         * \brief Compute the value of the spline given coefficients.
         * \param coefficients Vector of coefficients for the spline.
         * \param x The point at which to evaluate the spline.
         * \return The evaluated spline value.
         */
        double value(const Eigen::VectorXd& coefficients, double x) const;

      private:
        std::vector<double> knots_; /*!< Vector of knots. */
        Size degree_;               /*!< Degree of the spline. */
        Size numBasisFunctions_;    /*!< Number of basis functions. */
        mutable std::vector<Eigen::SparseMatrix<double>>
            Rk_matrices_;                /*!< Precomputed Rk matrices. */
        mutable Eigen::VectorXd tempB1_; /*!< Temporary vector for basis function evaluation. */
        mutable Eigen::VectorXd tempB2_; /*!< Temporary vector for basis function evaluation. */

        /*!
         * \brief Find the knot span for a given value.
         * \param x The value for which to find the knot span.
         * \return The index of the knot span.
         */
        Size findKnotSpan(double x) const;

        /*!
         * \brief Precompute the Rk matrices for the B-spline.
         */
        void precomputeRkMatrices() const;

        /*!
         * \brief Initialize temporary vectors used in evaluation.
         */
        void initializeTempVectors() const;

        /*!
         * \brief Evaluate the basis functions at a given value and knot span.
         * \param B Reference to a vector to store the evaluated basis functions.
         * \param x The value at which to evaluate the basis functions.
         * \param mu The knot span index.
         */
        void evaluate(Eigen::Ref<Eigen::VectorXd> B, double x, Size mu) const;
    };
}

#endif // b_spline_evaluator_hpp
