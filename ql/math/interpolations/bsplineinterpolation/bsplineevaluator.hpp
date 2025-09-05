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
    \brief B-spline basis function evaluator for multi-segment curves
*/

#ifndef quantlib_bspline_evaluator_hpp
#define quantlib_bspline_evaluator_hpp

#include "splinesegment.hpp"
#include <ql/shared_ptr.hpp>
#include <ql/types.hpp>
#include <vector>
#include <Eigen/Dense>

namespace QuantLib {

    /*!
     * \brief Evaluates B-spline basis functions for multi-segment curves
     * 
     * This class encapsulates the complex logic for evaluating basis functions
     * across multiple B-spline segments, handling boundary conditions and
     * sidedness correctly. It serves as the single source of truth for
     * basis evaluation, used by both BSplineStructure and StagedProblem.
     */
    class BSplineEvaluator {
    public:
        /*!
         * \brief Constructor with segments
         * \param segments The B-spline segments
         * \param numVariables Total number of variables (sum of segment variables)
         */
        BSplineEvaluator(const std::vector<ext::shared_ptr<BSplineSegment>>& segments,
                         Size numVariables);

        /*!
         * \brief Evaluate all basis functions at a point
         * 
         * This method handles:
         * - Multi-segment curves
         * - Boundary evaluation rules
         * - Sidedness for segment boundaries
         * 
         * \param x The evaluation point
         * \param side The evaluation side (RIGHT or LEFT)
         * \return Vector of basis function values
         */
        Eigen::VectorXd evaluateAll(Real x, 
                                    BSplineSegment::SideEnum side = BSplineSegment::SideRight) const;

        /*!
         * \brief Evaluate a specific derivative at a point
         * 
         * \param x The evaluation point
         * \param nu The derivative order (0 = value, 1 = first derivative, etc.)
         * \param side The evaluation side
         * \return Vector of derivative values
         */
        Eigen::VectorXd evaluateDerivative(Real x, 
                                           Integer nu,
                                           BSplineSegment::SideEnum side = BSplineSegment::SideRight) const;

        /*!
         * \brief Get the range of the multi-segment curve
         * \return Pair of (min, max) x values
         */
        std::pair<Real, Real> range() const;

        /*!
         * \brief Get the number of segments
         */
        Size getNumSegments() const { return segments_.size(); }

        /*!
         * \brief Get the total number of variables
         */
        Size getNumVariables() const { return numVariables_; }

    private:
        std::vector<ext::shared_ptr<BSplineSegment>> segments_;
        Size numVariables_;
        
        /*!
         * \brief Find which segment contains a given x value
         * \param x The point to locate
         * \param side The evaluation side
         * \return Index of the containing segment, or SIZE_MAX if not found
         */
        Size findSegmentIndex(Real x, BSplineSegment::SideEnum side) const;
    };

}

#endif // quantlib_bspline_evaluator_hpp