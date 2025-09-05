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

#include "bsplineevaluator.hpp"
#include <ql/errors.hpp>

namespace QuantLib {

    BSplineEvaluator::BSplineEvaluator(
        const std::vector<ext::shared_ptr<BSplineSegment>>& segments,
        Size numVariables)
        : segments_(segments), numVariables_(numVariables) {
        
        QL_REQUIRE(!segments_.empty(), "BSplineEvaluator: segments cannot be empty");
        
        // Verify total variables match
        Size totalVars = 0;
        for (const auto& segment : segments_) {
            totalVars += segment->getNumVariables();
        }
        
        QL_REQUIRE(numVariables_ == totalVars,
                   "BSplineEvaluator: numVariables (" << numVariables_ << 
                   ") doesn't match sum of segment variables (" << totalVars << ")");
    }

    Eigen::VectorXd BSplineEvaluator::evaluateAll(Real x, BSplineSegment::SideEnum side) const {
        QL_REQUIRE(side == BSplineSegment::SideRight || side == BSplineSegment::SideLeft,
                   "Sidedness needs to be 'Right' or 'Left'");
        
        const Size nSegments = segments_.size();
        Eigen::VectorXd result = Eigen::VectorXd::Zero(numVariables_);
        
        Size j = 0;  // Current position in result vector
        
        // Handle boundary evaluation correctly
        for (Size i = 0; i < nSegments; ++i) {
            bool inSegment = false;
            const auto& segment = segments_[i];
            const auto segmentRange = segment->range();
            
            if (side == BSplineSegment::SideRight) {
                // For right-sided evaluation:
                // - Include left boundary: x >= range.first
                // - Include right boundary for last segment: x <= range.second
                // - Exclude right boundary for other segments: x < range.second
                if (i == nSegments - 1) {
                    // Last segment: include right boundary
                    inSegment = (segmentRange.first <= x && x <= segmentRange.second);
                } else {
                    // Not last segment: exclude right boundary
                    inSegment = (segmentRange.first <= x && x < segmentRange.second);
                }
            } else {
                // For left-sided evaluation:
                // - Exclude left boundary for non-first segments: x > range.first
                // - Include left boundary for first segment: x >= range.first
                // - Include right boundary: x <= range.second
                if (i == 0) {
                    // First segment: include left boundary
                    inSegment = (segmentRange.first <= x && x <= segmentRange.second);
                } else {
                    // Not first segment: exclude left boundary
                    inSegment = (segmentRange.first < x && x <= segmentRange.second);
                }
            }
            
            if (inSegment) {
                const Eigen::VectorXd segmentResult = segment->evaluateAll(x, -1, side);
                const Size segVars = segment->getNumVariables();
                
                // Safety check before assignment
                QL_REQUIRE(j + segVars <= numVariables_,
                           "Segment assignment would exceed vector bounds: trying to assign " + 
                           std::to_string(segVars) + " values at position " + 
                           std::to_string(j) + " in vector of size " + 
                           std::to_string(numVariables_));
                
                result.segment(j, segVars) = segmentResult;
            }
            j += segment->getNumVariables();
        }
        
        return result;
    }

    Eigen::VectorXd BSplineEvaluator::evaluateDerivative(Real x, Integer nu,
                                                         BSplineSegment::SideEnum side) const {
        QL_REQUIRE(side == BSplineSegment::SideRight || side == BSplineSegment::SideLeft,
                   "Sidedness needs to be 'Right' or 'Left'");
        
        const Size nSegments = segments_.size();
        Eigen::VectorXd result = Eigen::VectorXd::Zero(numVariables_);
        
        Size j = 0;
        
        for (Size i = 0; i < nSegments; ++i) {
            bool inSegment = false;
            const auto& segment = segments_[i];
            const auto segmentRange = segment->range();
            
            // Same boundary logic as evaluateAll
            if (side == BSplineSegment::SideRight) {
                if (i == nSegments - 1) {
                    inSegment = (segmentRange.first <= x && x <= segmentRange.second);
                } else {
                    inSegment = (segmentRange.first <= x && x < segmentRange.second);
                }
            } else {
                if (i == 0) {
                    inSegment = (segmentRange.first <= x && x <= segmentRange.second);
                } else {
                    inSegment = (segmentRange.first < x && x <= segmentRange.second);
                }
            }
            
            if (inSegment) {
                const Eigen::VectorXd segmentResult = segment->evaluateAll(x, nu, side);
                const Size segVars = segment->getNumVariables();
                
                QL_REQUIRE(j + segVars <= numVariables_,
                           "Segment assignment would exceed vector bounds");
                
                result.segment(j, segVars) = segmentResult;
            }
            j += segment->getNumVariables();
        }
        
        return result;
    }

    std::pair<Real, Real> BSplineEvaluator::range() const {
        QL_REQUIRE(!segments_.empty(), "No segments available");
        
        auto firstRange = segments_.front()->range();
        auto lastRange = segments_.back()->range();
        
        return std::make_pair(firstRange.first, lastRange.second);
    }

    Size BSplineEvaluator::findSegmentIndex(Real x, BSplineSegment::SideEnum side) const {
        const Size nSegments = segments_.size();
        
        for (Size i = 0; i < nSegments; ++i) {
            bool inSegment = false;
            const auto segmentRange = segments_[i]->range();
            
            if (side == BSplineSegment::SideRight) {
                if (i == nSegments - 1) {
                    inSegment = (segmentRange.first <= x && x <= segmentRange.second);
                } else {
                    inSegment = (segmentRange.first <= x && x < segmentRange.second);
                }
            } else {
                if (i == 0) {
                    inSegment = (segmentRange.first <= x && x <= segmentRange.second);
                } else {
                    inSegment = (segmentRange.first < x && x <= segmentRange.second);
                }
            }
            
            if (inSegment) {
                return i;
            }
        }
        
        return SIZE_MAX;  // Not found
    }

}