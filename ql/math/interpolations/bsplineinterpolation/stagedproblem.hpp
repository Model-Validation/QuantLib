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

/*! \file stagedproblem.hpp
    \brief Staged interpolation problem for warm-start optimization
*/

#ifndef quantlib_staged_problem_hpp
#define quantlib_staged_problem_hpp

#include "splineconstraints.hpp"
#include "splinesegment.hpp"
#include "multisegmentbsplineevaluator.hpp"
#include <ql/shared_ptr.hpp>
#include <ql/types.hpp>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace QuantLib {

    /*!
     * \brief Interpolation mode for each point
     */
    enum class InterpolationMode {
        HARD,   // Equality constraint (exact interpolation)
        LS,     // Least squares (soft constraint)
        AUTO    // Automatically determine based on DOF
    };

    /*!
     * \brief Mode specification for a coordinate span
     */
    struct ModeSpan {
        Real start;
        Real end;
        InterpolationMode mode;
        
        // Default constructor for SWIG vector support
        ModeSpan() : start(0.0), end(0.0), mode(InterpolationMode::AUTO) {}
        
        ModeSpan(Real s, Real e, InterpolationMode m) 
            : start(s), end(e), mode(m) {}
    };

    /*!
     * \brief Generate Voronoi-like mode spans based on data density
     * 
     * Creates spans where:
     * - Dense data regions use LS mode (smooth fitting)
     * - Sparse data regions use HARD mode (exact interpolation)
     * - Boundary regions use HARD mode (preserve exact values)
     * 
     * \param dataPoints The x-coordinates of data points
     * \param domainStart Start of the interpolation domain
     * \param domainEnd End of the interpolation domain
     * \param densityThreshold Points per unit length threshold for LS vs HARD
     * \param boundaryWidth Width of boundary regions (from domain edges)
     * \return Vector of ModeSpan objects covering the domain
     */
    std::vector<ModeSpan> generateVoronoiModeSpans(
        const std::vector<Real>& dataPoints,
        Real domainStart,
        Real domainEnd,
        Real densityThreshold = 2.0,
        Real boundaryWidth = 0.1
    );

    /*!
     * \brief Staged problem for B-spline interpolation with warm-start capability
     * 
     * This class separates the base spline constraints from interpolation amendments,
     * allowing the interpolation structure to be built once and reused with different
     * y values. This enables warm-start optimization in the SCS solver.
     */
    class StagedProblem {
    public:
        /*!
         * \brief Constructor with base constraints
         * \param baseConstraints The immutable base spline constraints
         */
        explicit StagedProblem(const ext::shared_ptr<SplineConstraints>& baseConstraints);

        /*!
         * \brief Stage the interpolation structure for given x points and modes
         * 
         * This builds the interpolation layer once. After staging, solve() can be
         * called multiple times with different y values for warm-start benefits.
         * 
         * \param interpolationNodes The x coordinates for interpolation
         * \param modes The mode for each interpolation point
         * \param segments The spline segments for basis evaluation
         */
        void stage(const std::vector<Real>& interpolationNodes,
                   const std::vector<InterpolationMode>& modes,
                   const std::vector<ext::shared_ptr<BSplineSegment>>& segments);

        /*!
         * \brief Stage with mode spans instead of per-point modes
         * 
         * \param interpolationNodes The x coordinates for interpolation
         * \param modeSpans Coordinate ranges with associated modes
         * \param segments The spline segments for basis evaluation
         */
        void stage(const std::vector<Real>& interpolationNodes,
                   const std::vector<ModeSpan>& modeSpans,
                   const std::vector<ext::shared_ptr<BSplineSegment>>& segments);

        /*!
         * \brief Solve the staged problem with given y values
         * 
         * This uses warm-start if called multiple times after staging.
         * The structure built in stage() is reused.
         * 
         * \param values The y values for interpolation
         * \return Solution coefficients
         */
        Eigen::VectorXd solve(const std::vector<Real>& values);

        /*!
         * \brief Check if the problem has been staged
         */
        bool isStaged() const { return staged_; }

        /*!
         * \brief Get the number of hard constraints
         */
        Size getNumHardConstraints() const { return hardIndices_.size(); }

        /*!
         * \brief Get the number of soft constraints
         */
        Size getNumSoftConstraints() const { return softIndices_.size(); }

        /*!
         * \brief Clear the staged structure (forces re-staging)
         */
        void clearStaging();

    private:
        // Layer 1: Base problem (immutable after construction)
        ext::shared_ptr<SplineConstraints> baseConstraints_;
        Size baseNumConstraints_;
        Size baseNumEqualities_;
        Size baseNumInequalities_;
        
        // Layer 2: Interpolation layer (built once in stage())
        Eigen::SparseMatrix<Real> A_hard_;  // Hard constraint matrix
        Eigen::SparseMatrix<Real> Q_soft_;  // Soft points quadratic form
        std::vector<Size> hardIndices_;     // Indices of hard points
        std::vector<Size> softIndices_;     // Indices of soft points
        
        // Combined problem for solving
        ext::shared_ptr<SplineConstraints> combinedConstraints_;
        
        // Staging state
        bool staged_;
        std::vector<Real> lastX_;  // Remember x points for validation
        std::vector<ext::shared_ptr<BSplineSegment>> segments_;  // Remember segments
        ext::shared_ptr<MultiSegmentBSplineEvaluator> evaluator_;  // Evaluator for basis functions
        
        // Helper methods
        InterpolationMode getModeAt(Real x, const std::vector<ModeSpan>& spans) const;
        
        void buildInterpolationMatrices(
            const std::vector<Real>& interpolationNodes,
            const std::vector<InterpolationMode>& modes,
            const std::vector<ext::shared_ptr<BSplineSegment>>& segments);
            
        void combineConstraints();
    };

}

#endif // quantlib_staged_problem_hpp