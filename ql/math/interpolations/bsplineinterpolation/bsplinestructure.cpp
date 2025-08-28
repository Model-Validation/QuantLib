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

// TODO: Lets ReSharper relax about the type mismatch between QuantLib::Size and Eigen::Index, should revisit this
// ReSharper disable CppClangTidyBugproneNarrowingConversions
#include "bsplinestructure.hpp"
#include "splineconstraints.hpp"
#include <ql/errors.hpp>
#include <ql/shared_ptr.hpp>
#include <ql/types.hpp>
#include <algorithm>
#include <cmath>
#include <iosfwd>
#include <iostream>
#include <type_traits>
#include <utility>
#include <vector>

namespace QuantLib {
    /*
    Class for a composite spline structure holding the necessary information on knots, degree and constraints
    for a spline.
    */

    // Note that constraints on splineSegments are ignored, you have to pass the composed spline constraints, this
    // can be done on the python side. A rewrite of this would trim down the
    // single BSplineStructure to not include that, and require to wrap it in here with the
    // constraints.
    BSplineStructure::BSplineStructure(
        const std::vector<ext::shared_ptr<BSplineSegment>>& splineSegments,
        const ext::shared_ptr<SplineConstraints>& splineConstraints,
        bool useSegmentNodes,
        bool rejectZeroNode) : splineSegments_(splineSegments),
                               splineConstraints_(splineConstraints),
                               rejectZeroNode_(rejectZeroNode), useSegmentNodes_(useSegmentNodes) {
        // Create a vector to store the second element of the range for all but the last segment
        segmentNodes_.reserve(splineSegments.size() - 1); // Reserve space for efficiency

        for (Size i = 0; i < splineSegments.size() - 1; ++i) {
            auto range = splineSegments[i]->range(); // Assume range() returns std::pair
            segmentNodes_.push_back(range.second); // Add the second element to the vector
        }
    }

    void BSplineStructure::addInterpolationNodes(const std::vector<Real>& interpolationNodes,
                                                 BSplineSegment::SideEnum side,
                                                 Size nParameters) const {
        const Size nInterpolationNodes = interpolationNodes.size();
        const Size nConstraints = splineConstraints_->getNConstraints();
        const Size nVariables = splineConstraints_->getNumVariables();

        for (const Real interpolationNode : interpolationNodes) {
            Eigen::VectorXd row = evaluateAll(interpolationNode, side);
            // Use the new method that maintains SCS ordering
            splineConstraints_->addEqualityConstraintAtBeginning(row, 0.0);
        }

        // We still have the old constraint count
        Eigen::SparseMatrix<Real> B_new(nConstraints + nInterpolationNodes,
                                        nParameters + nInterpolationNodes);
        for (Size i = 0; i < nInterpolationNodes; ++i) {
            B_new.insert(nConstraints + i, nParameters + i) = 1.0;
        }

        Eigen::SparseMatrix<double> C_new(nVariables, nParameters + nInterpolationNodes);

        splineConstraints_->addParameters(nInterpolationNodes, B_new, C_new);
    }

    // Here we interpolate by adding linear constraints that enforce the interpolation
    // If the system would become overdetermined (there already may be conditions for continuity and such,
    // as for Hermite interpolation), then one can instead solve a quadratic problem finding the point
    // minimizing the norm squared of A x - b, rather than solving it. We can indeed take any part of this
    // system and move it to the objective function for the same purpose. This is hacked in here post-hoc, but
    // awaits a better setup.
    Eigen::VectorXd
    BSplineStructure::interpolate(const std::vector<Real>& interpolationNodesOrg,
                                  const std::vector<Real>& valuesOrg) {
        std::vector<Real> interpolationNodes, values;
        // TODO: this is a hack, bootstrapper adds 0.0 as a node, which may conflict with our setup,
        // e.g. if we have an extrapolation to the left
        if (interpolationNodesOrg[0] == 0.0 && rejectZeroNode_) {
            interpolationNodes.reserve(interpolationNodesOrg.size() - 1);
            values.reserve(interpolationNodesOrg.size() - 1);
            interpolationNodes.insert(interpolationNodes.end(), interpolationNodesOrg.begin() + 1,
                                      interpolationNodesOrg.end());
            values.insert(values.end(), valuesOrg.begin() + 1, valuesOrg.end());
        } else {
            interpolationNodes.reserve(interpolationNodesOrg.size());
            values.reserve(interpolationNodesOrg.size());
            interpolationNodes.insert(interpolationNodes.end(), interpolationNodesOrg.begin(),
                                      interpolationNodesOrg.end());
            values.insert(values.end(), valuesOrg.begin(), valuesOrg.end());
        }
        std::vector<Real> transformedValues = transform(interpolationNodes, values);

        // The segment nodes is another semi-hack, the idea to assign special nodes to have a special role actually
        // is a reasonable way to push non-linearity to the bootstrapper and keep the curve problem in each
        // iteration convex. But it is probably better to assign these nodes explicitly in the class. Here
        // we use this to allow for different interpolation spaces in different segments, and tie them
        // by requiring transformed values from left and right to be the same. This was done to incorporate RT
        // interpolation.
        if (useSegmentNodes_) {
            std::vector<Real> segmentNodes, segmentNodeValues;
            //extendedNodes.reserve(interpolationNodes.size());
            //extendedNodes.insert(extendedNodes.end(), interpolationNodes.begin(),
            //                     interpolationNodes.end());
            segmentNodes.reserve(segmentNodes_.size());
            segmentNodeValues.reserve(segmentNodes_.size());

            auto lower = interpolationNodes.begin();

            for (Real segmentNode : segmentNodes_) {
                lower = std::lower_bound(lower, interpolationNodes.end(), segmentNode - tolerance_);
                if (lower != interpolationNodes.end() && *lower <= segmentNode + tolerance_) {
                    const Size index = std::distance(interpolationNodes.begin(), lower);
                    segmentNodes.emplace_back(segmentNode);
                    segmentNodeValues.emplace_back(values[index]);
                }
            }

            //extendedNodes.reserve(extendedNodes.size() + segmentNodes.size());
            //extendedNodes.insert(extendedNodes.end(), segmentNodes.begin(), segmentNodes.end());

            std::vector<Real> transformedValues2 =
                transform(segmentNodes, segmentNodeValues, BSplineSegment::SideLeft);

            transformedValues.reserve(transformedValues.size() + transformedValues2.size());
            transformedValues.insert(transformedValues.end(), transformedValues2.begin(),
                                     transformedValues2.end());
            addInterpolationNodes(segmentNodes, BSplineSegment::SideLeft,
                                  interpolationNodes.size());
        } else {
            transformedValues = transform(interpolationNodes, values);
        }

        const Size nConstraintsBefore = splineConstraints_->getNConstraints();
        splineConstraints_->push();
        addInterpolationNodes(interpolationNodes);

        Eigen::VectorXd solution;
        if (splineConstraints_->fitData_) {
            const Integer firstRow = static_cast<Integer>(nConstraintsBefore);
            const Integer lastRow =
                static_cast<Integer>(nConstraintsBefore + transformedValues.size());
            interpolationA_ = splineConstraints_->getSliceOfA(firstRow, lastRow);

            splineConstraints_->pop();

            interpolationBVec_ = Eigen::Map<const Eigen::VectorXd>(transformedValues.data(),
                transformedValues.size());

            if (firstRow == 0) {
                // No pre-existing constraints - pure unconstrained LS
                Eigen::SparseQR<Eigen::SparseMatrix<Real>, Eigen::COLAMDOrdering<int>> solver;
                solver.compute(interpolationA_);
                solution = solver.solve(interpolationBVec_);
                QL_REQUIRE(solver.info() == Eigen::Success, "Solver failed");
            } else {
                // FIXED: Build LS objective and solve WITH constraints
                // The LS objective is: min 0.5 * x'Px + c'x where P = A'A, c = -A'b
                Eigen::SparseMatrix<Real> P = interpolationA_.transpose() * interpolationA_;
                Eigen::VectorXd c = -(interpolationA_.transpose() * interpolationBVec_);
                
                // Update the objective in the constraint system
                splineConstraints_->setP(P);
                splineConstraints_->setCVector(c);
                
                // CRITICAL FIX: The pre-existing constraints (including inequalities)
                // are still in splineConstraints_ and will be respected by solve()
                // We pass empty parameters since we're using the LS objective
                solution = solve(std::vector<Real>(0));
            }
        } else {
            solution = solve(transformedValues);

            splineConstraints_->pop();
        }
        return solution;
    }

    // ReSharper disable once CppInconsistentNaming
    std::vector<Real> BSplineStructure::interpolate_swig(
        const std::vector<Real>& interpolationNodes,
        const std::vector<Real>& values) {
        Eigen::VectorXd result = interpolate(interpolationNodes, values);
        return {result.data(), result.data() + result.size()};
    }

    std::vector<std::vector<Real>> BSplineStructure::get_interpolation_a() const {
        // Initialize the output vector
        std::vector<std::vector<Real>> result(this->interpolationA_.rows(),
                                              std::vector<Real>(this->interpolationA_.cols(), 0.0));

        // Convert the sparse matrix to a dense std::vector<std::vector<Real>>
        for (int k = 0; k < this->interpolationA_.outerSize(); ++k) {
            for (Eigen::SparseMatrix<Real>::InnerIterator it(this->interpolationA_, k); it; ++it) {
                result[it.row()][it.col()] = it.value();
            }
        }

        return result;
    }

    std::vector<Real> BSplineStructure::get_interpolation_b() const {
        return {interpolationBVec_.data(), interpolationBVec_.data() + interpolationBVec_.size()};
    }

    Eigen::VectorXd BSplineStructure::solve(const std::vector<Real>& parameters) const {
        splineConstraints_->update_b(parameters);
        int status = splineConstraints_->solve();
        QL_REQUIRE(status == 1, "Solution failed, returned " << status << '\n');
        return splineConstraints_->getSolution();
    }

    // ReSharper disable once CppInconsistentNaming
    std::vector<Real> BSplineStructure::solve_swig(const std::vector<Real>& parameters) const {
        Eigen::VectorXd result = solve(parameters);
        return {result.data(), result.data() + result.size()};
    }

    Eigen::VectorXd BSplineStructure::evaluateAll(const Real x,
                                                  const BSplineSegment::SideEnum side) const {
        QL_REQUIRE(side == BSplineSegment::SideRight || side == BSplineSegment::SideLeft,
                   "Sidedness needs to be 'Right' or 'Left");
        
        const Size nSegments = splineSegments_.size();
        const Size nVariables = splineConstraints_->getNumVariables();
        Size j = 0;
        Eigen::VectorXd result = Eigen::VectorXd::Zero(nVariables);
        
        // Handle boundary evaluation correctly
        for (Size i = 0; i < nSegments; ++i) {
            bool inSegment = false;
            const auto& segment = splineSegments_[i];
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
                result.segment(j, segment->getNumVariables()) = segmentResult;
            }
            j += segment->getNumVariables();
        }

        return result;
    }

    Real BSplineStructure::value(const Eigen::VectorXd& coefficients,
                                 Real x,
                                 Integer nu,
                                 BSplineSegment::SideEnum side) const {
        const Size numSegments = this->splineSegments_.size();
        const Size e = nu >= 0 ? 0 : static_cast<Size>(-nu);
        Size mu = 0;
        Size cumVariables = 0;

        Eigen::VectorXd extraCoefficients = Eigen::VectorXd::Zero(e);
        Eigen::VectorXd augmentedCoefficients;
        for (Size i = 0; i < numSegments; ++i) {
            const auto splineSegment = this->splineSegments_[i];
            const auto segmentRange = splineSegment->range();
            bool inSegment = false;
            
            // Apply same boundary logic as in evaluateAll
            if (i == numSegments - 1) {
                // Last segment: include right boundary
                inSegment = (segmentRange.first <= x && x <= segmentRange.second);
            } else {
                // Not last segment: exclude right boundary
                inSegment = (segmentRange.first <= x && x < segmentRange.second);
            }
            
            if (inSegment) {
                mu = i;
                break;
            }
            const Size nVars = splineSegment->getNumVariables();
            augmentedCoefficients = Eigen::VectorXd::Zero(nVars + e);
            augmentedCoefficients.segment(0, nVars) = coefficients.segment(cumVariables, nVars);
            augmentedCoefficients.segment(nVars, e) = extraCoefficients;
            for (Size j = 0; j < e; ++j) {
                extraCoefficients[j] = splineSegment->value(
                    augmentedCoefficients, splineSegment->range().second,
                    nu + static_cast<Integer>(j), BSplineSegment::SideLeft);
            }
            cumVariables += nVars;
        }
        const auto splineSegment = this->splineSegments_[mu];
        const Size nVars = splineSegment->getNumVariables();

        augmentedCoefficients = Eigen::VectorXd::Zero(nVars + e);
        augmentedCoefficients.segment(0, nVars) = coefficients.segment(cumVariables, nVars);
        augmentedCoefficients.segment(nVars, e) = extraCoefficients;

        return splineSegment->value(augmentedCoefficients, x, nu, side);
    }

    Real BSplineStructure::value(const std::vector<Real>& coefficients,
                                 Real x,
                                 Integer nu,
                                 BSplineSegment::SideEnum side) const {
        return value(Eigen::Map<const Eigen::VectorXd>(coefficients.data(), coefficients.size()), x,
                     nu, side);
    }

    std::pair<Real, Real> BSplineStructure::range() const {
        Real left = this->splineSegments_.front()->range().first;
        Real right = this->splineSegments_.back()->range().second;

        return {left, right};
    }

    std::vector<Real> BSplineStructure::transform(const std::vector<Real>& abscissae,
                                                  const std::vector<Real>& values,
                                                  BSplineSegment::SideEnum side) const {

        // Check if abscissae is sorted and vectors are of same length
        // TODO: remove for performance reasons, do this check only where needed
        QL_REQUIRE(std::is_sorted(abscissae.begin(), abscissae.end()),
                   "Abscissae must be sorted in ascending order");
        QL_REQUIRE(abscissae.size() == values.size(),
                   "Abscissae and values must be of the same length");

        std::vector<Real> transformedValues(abscissae.size());

        // Prepare the transformed values vector
        transformedValues.resize(abscissae.size());

        Size segmentIndex = 0;
        const Size numSegments = splineSegments_.size();

        // Iterate over abscissae and values
        for (Size i = 0; i < abscissae.size(); ++i) {
            const Real x = abscissae[i];
            const Real value = values[i];

            // Find the correct spline segment
            while (segmentIndex < numSegments) {
                auto [left, right] = splineSegments_[segmentIndex]->range();
                // TODO: endpoints are not handled correctly, use side variables
                if ((x < left && segmentIndex == 0) || (
                        (x >= left && x < right) && side == BSplineSegment::SideRight) ||
                    ((x > left && x <= right) && side == BSplineSegment::SideLeft) ||
                    (segmentIndex == numSegments - 1 && x >= right)) {
                    // Apply the transformation
                    transformedValues[i] = splineSegments_[segmentIndex]->transform(x, value);
                    break;
                }
                // TODO: More careful here, not sure if this is safe without further checks, e.g. if x < left
                // Move to the next segment
                if (x < left) {
                    QL_FAIL("x = " << x << " is on the left of the range [" << left
                        << ", " << right << "], this should not happen.");
                }
                ++segmentIndex;
            }
        }
        return transformedValues;
    }

    Eigen::VectorXd BSplineStructure::getSolution() const {
        return splineConstraints_->getSolution();
    }

    // ReSharper disable once CppInconsistentNaming
    std::vector<Real> BSplineStructure::get_solution() const {
        Eigen::VectorXd result = splineConstraints_->getSolution();
        return {result.data(), result.data() + result.size()};
    }

}