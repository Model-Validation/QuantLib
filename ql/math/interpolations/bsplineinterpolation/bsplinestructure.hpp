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

/*! \file bsplinestructure.hpp
    \brief Structure to manage composite B-spline segments
*/

#ifndef composite_spline_structure_hpp
#define composite_spline_structure_hpp

#include "bsplineevaluator.hpp"
#include "splineconstraints.hpp"
#include "splinesegment.hpp"
#include "stagedproblem.hpp"
#include <ql/shared_ptr.hpp>
#include <ql/types.hpp>
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace QuantLib {

    /*!
     * \brief Class for managing a composite B-spline structure.
     */
    class BSplineStructure {
      public:
        BSplineStructure(
            const std::vector<ext::shared_ptr<BSplineSegment>>& splineSegments,
            const ext::shared_ptr<SplineConstraints>& splineConstraints,
            bool useSegmentNodes = false,
            bool rejectZeroNode = true);
        BSplineStructure() = default;

        std::vector<Real> transform(const std::vector<Real>& abscissae,
                                    const std::vector<Real>& values,
                                    BSplineSegment::SideEnum side = BSplineSegment::SideRight) const;
        void addInterpolationNodes(const std::vector<Real>& interpolationNodes,
                                   BSplineSegment::SideEnum side = BSplineSegment::SideRight, Size nParameters = 0) const;
        Size getNumVariables() const { return splineConstraints_->getNumVariables(); }
        Integer getNumVariablesSwig() const { return static_cast<Integer>(getNumVariables()); }
        std::vector<Real> evaluateAllSwig(Real x,BSplineSegment::SideEnum side) const {
            Eigen::VectorXd result = evaluateAll(x, side);
            return {result.data(), result.data() + result.size()};
        }


        Real value(const Eigen::VectorXd& coefficients,
                   Real x,
                   Integer nu = 0,
                   BSplineSegment::SideEnum side = BSplineSegment::SideRight) const;
        Real value(const std::vector<Real>& coefficients,
                   Real x,
                   Integer nu = 0,
                   BSplineSegment::SideEnum side = BSplineSegment::SideRight) const;
        const std::vector<ext::shared_ptr<BSplineSegment>>& getSplineSegments() const { return splineSegments_; }
        std::vector<ext::shared_ptr<BSplineSegment>> getSplineSegmentsSwig() const { return splineSegments_; }
        std::pair<Real, Real> range() const;
        Eigen::VectorXd interpolate(const std::vector<Real>& interpolationNodes,
                                    const std::vector<Real>& values);
        std::vector<Real> interpolate_swig(const std::vector<Real>& interpolationNodes,
                                           const std::vector<Real>& values);
        std::vector<std::vector<Real>> get_interpolation_a() const;
        std::vector<Real> get_interpolation_b() const;
        Eigen::VectorXd solve(const std::vector<Real>& parameters) const;
        std::vector<Real> solve_swig(const std::vector<Real>& parameters) const;
        void setConstraints(const ext::shared_ptr<SplineConstraints>& splineConstraints) {
            this->splineConstraints_ = splineConstraints;
        }
        ext::shared_ptr<SplineConstraints> getConstraints() const { return this->splineConstraints_; }
        Eigen::VectorXd evaluateAll(Real x, BSplineSegment::SideEnum side = BSplineSegment::SideRight) const;
        Eigen::VectorXd getSolution() const;
        std::vector<Real> get_solution() const;
        
        // Staged interpolation methods
        void enableStaging(bool enable = true) { useStaging_ = enable; }
        bool isStagingEnabled() const { return useStaging_; }
        void clearStaging() { if (stagedProblem_) stagedProblem_->clearStaging(); }
        
        // Interpolate with mode specification
        Eigen::VectorXd interpolate(const std::vector<Real>& interpolationNodes,
                                    const std::vector<Real>& values,
                                    const std::vector<InterpolationMode>& modes);
        
        // Friend class for unit testing internal state
        friend class BSplineTestAccess;
        
      private:
        std::vector<ext::shared_ptr<BSplineSegment>> splineSegments_;
        ext::shared_ptr<SplineConstraints> splineConstraints_;
        bool rejectZeroNode_ = true;
        bool useSegmentNodes_ = false;
        std::vector<Real> segmentNodes_;
        BSplineEvaluator spline_;
        Real tolerance_ = 0.0e-14;
        Eigen::SparseMatrix<Real> interpolationA_;
        Eigen::VectorXd interpolationBVec_;
        
        // Staging support
        bool useStaging_ = false;
        mutable ext::shared_ptr<StagedProblem> stagedProblem_;
        mutable std::vector<Real> lastStagedX_;
    };
}
#endif // composite_spline_structure_hpp
