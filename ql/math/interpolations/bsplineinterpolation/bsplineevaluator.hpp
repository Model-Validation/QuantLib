/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2025

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.
*/

/*! \file bsplineevaluator.hpp
    \brief Single-segment B-spline evaluator for BSplineSegment internal use
*/

#ifndef quantlib_single_segment_bspline_evaluator_hpp
#define quantlib_single_segment_bspline_evaluator_hpp

#include <ql/types.hpp>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace QuantLib {

    /*!
     * \brief Single-segment B-spline basis function evaluator
     * 
     * This is a simple evaluator for individual B-spline segments,
     * used internally by BSplineSegment. This is different from
     * MultiSegmentBSplineEvaluator which handles multi-segment curves.
     */
    class BSplineEvaluator {
    public:
        /*! Default constructor */
        BSplineEvaluator();
        
        /*!
         * \brief Constructor with knots and degree
         * \param knots The knot vector
         * \param degree The degree of the B-spline
         */
        BSplineEvaluator(const std::vector<Real>& knots, Integer degree);
        
        /*!
         * \brief Evaluate all basis functions at a point
         * \param x The evaluation point
         * \return Vector of basis function values
         */
        Eigen::VectorXd evaluateAll(Real x) const;
        
        /*!
         * \brief Evaluate the B-spline value given coefficients
         * \param coefficients The B-spline coefficients
         * \param x The evaluation point
         * \return The interpolated value
         */
        Real value(const Eigen::VectorXd& coefficients, Real x) const;
        
    private:
        std::vector<Real> knots_;
        Integer degree_;
        Size numBasisFunctions_;
        mutable std::vector<Eigen::SparseMatrix<Real>> Rk_matrices_;
        mutable Eigen::VectorXd tempB1_;
        mutable Eigen::VectorXd tempB2_;
        
        Size findKnotSpan(Real x) const;
        void evaluate(Eigen::Ref<Eigen::VectorXd> B, Real x, Size mu) const;
        void precomputeRkMatrices() const;
        void initializeTempVectors() const;
    };

}

#endif // quantlib_single_segment_bspline_evaluator_hpp