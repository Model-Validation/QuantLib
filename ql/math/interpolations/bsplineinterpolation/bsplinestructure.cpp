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
                               rejectZeroNode_(rejectZeroNode), 
                               useSegmentNodes_(useSegmentNodes),
                               spline_(splineSegments, splineConstraints->getNumVariables()) {
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
        const Size nEqualitiesBefore = splineConstraints_->getNumEqualities();

        // Mode-specific constraint handling
        if (splineConstraints_->fitData_) {
            // LS MODE: Do NOT add equality constraints
            // Instead, we'll add quadratic penalty terms to the objective
            // The constraints are handled as soft penalties in the objective function
        } else {
            // HARD MODE: Add interpolation equality constraints
            for (const Real interpolationNode : interpolationNodes) {
                Eigen::VectorXd row = evaluateAll(interpolationNode, side);
                
                QL_REQUIRE(row.size() == nVariables,
                          "evaluateAll returned wrong size: " << row.size() <<
                          " expected " << nVariables << " for x=" << interpolationNode);
                
                // Use the new method that maintains SCS ordering
                splineConstraints_->addEqualityConstraintAtBeginning(row, 0.0);
            }
        }

        // Set up parameter mapping for interpolation nodes
        // Hard mode: parameters map to constraint RHS via B matrix  
        // LS mode: parameters map to objective via C matrix
        {
            // Get updated constraint count after potentially adding constraints
            const Size nConstraintsAfter = splineConstraints_->getNConstraints();
            
            if (splineConstraints_->fitData_) {
                // LS MODE: Parameters map to objective via C matrix for warm-start
                // The objective is: min ||Ax - y||² = min x'(A'A)x - 2(A'y)'x
                // We need C matrix such that C*params gives us -A'y term
                
                const Size totalParameters = nParameters + nInterpolationNodes;
                
                // In LS mode, constraints don't have parameters (B is empty)
                // But objective has parameters via C matrix
                Eigen::SparseMatrix<Real> B_new(nConstraintsAfter, totalParameters);  // Empty B
                Eigen::SparseMatrix<Real> C_new(nVariables, totalParameters);
                
                // Build C matrix: C*params should give -A'*y contribution to objective
                // where A is the interpolation constraint matrix and y are the parameter values
                
                // In LS mode, the objective becomes:
                // minimize (1/2) x'Px + c'x - params'A'x
                // This is equivalent to: minimize (1/2) x'Px + (c - A'params)'x
                // So C matrix should be -A' (negative transpose of interpolation matrix)
                
                std::vector<Eigen::Triplet<Real>> C_triplets;
                for (Size i = 0; i < nInterpolationNodes; ++i) {
                    Real x = interpolationNodes[i];
                    Eigen::VectorXd basis = evaluateAll(x, side);
                    
                    for (Size j = 0; j < basis.size(); ++j) {
                        if (std::abs(basis[j]) > 1e-14) {
                            // C matrix: maps parameters to objective linear term
                            // Add -A'[i,j] as coefficient for parameter i affecting variable j
                            C_triplets.emplace_back(j, nParameters + i, -basis[j]);
                        }
                    }
                }
                
                C_new.setFromTriplets(C_triplets.begin(), C_triplets.end());
                
                // TODO: Add A'A term to quadratic objective for least squares
                // This creates the ||Ax - y||² objective: min x'(A'A)x - 2y'Ax
                // Currently not implemented due to SplineConstraints API limitations
                // The A'A term needs to be added to P matrix during SplineConstraints construction
                
                // Add parameters for warm-start capability in LS mode
                splineConstraints_->addParameters(nInterpolationNodes, B_new, C_new);
                
            } else {
                // HARD MODE: Parameters map to constraint RHS via B matrix
                // Constraints: Ax = b + B*params where params are y-values
                
                const Size totalParameters = nParameters + nInterpolationNodes;
                const Size totalConstraints = nConstraintsAfter;
                
                // Validate dimensions
                QL_REQUIRE(nEqualitiesBefore + nInterpolationNodes <= totalConstraints,
                           "Invalid row indices: nEqualitiesBefore=" << nEqualitiesBefore 
                           << " + nInterpolationNodes=" << nInterpolationNodes 
                           << " > totalConstraints=" << totalConstraints);
                
                Eigen::SparseMatrix<Real> B_new(totalConstraints, totalParameters);
                
                // Map parameters to the correct constraint rows (after reordering)
                for (Size i = 0; i < nInterpolationNodes; ++i) {
                    // The i-th interpolation constraint is now at row (nEqualitiesBefore + i)
                    // because addEqualityConstraintAtBeginning inserts after existing equalities
                    Size rowIndex = nEqualitiesBefore + i;
                    Size colIndex = nParameters + i;
                    
                    QL_REQUIRE(rowIndex < totalConstraints && colIndex < totalParameters,
                               "B matrix index out of bounds: (" << rowIndex << "," << colIndex 
                               << ") for matrix " << totalConstraints << "x" << totalParameters);
                    
                    B_new.insert(rowIndex, colIndex) = 1.0;
                }

                Eigen::SparseMatrix<Real> C_new(nVariables, totalParameters);  // Empty C matrix

                
                splineConstraints_->addParameters(nInterpolationNodes, B_new, C_new);
            }
        }
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
        // Check if we should use staging automatically
        if (useStaging_) {
            // Use AUTO mode for all points
            std::vector<InterpolationMode> autoModes(interpolationNodesOrg.size(), InterpolationMode::AUTO);
            return interpolate(interpolationNodesOrg, valuesOrg, autoModes);
        }

        // Legacy path: Use global fitData_ flag to determine mode for all points
        InterpolationMode globalMode = splineConstraints_->fitData_ ? InterpolationMode::LS : InterpolationMode::HARD;
        std::vector<InterpolationMode> globalModes(interpolationNodesOrg.size(), globalMode);
        return interpolate(interpolationNodesOrg, valuesOrg, globalModes);
    }

    Eigen::VectorXd BSplineStructure::interpolate(const std::vector<Real>& interpolationNodesOrg,
                                                  const std::vector<Real>& valuesOrg,
                                                  const std::vector<InterpolationMode>& modes) {
        QL_REQUIRE(interpolationNodesOrg.size() == valuesOrg.size(), 
                   "Number of interpolation nodes must match number of values");
        QL_REQUIRE(interpolationNodesOrg.size() == modes.size(),
                   "Number of interpolation nodes must match number of modes");
        
        // Always use staged approach for proper per-point mode handling
        // The staged path is the only one that correctly implements mixed HARD/LS modes
        if (!stagedProblem_) {
            stagedProblem_ = ext::make_shared<StagedProblem>(splineConstraints_);
        }
        
        // Check if we need to stage or can reuse
        if (lastStagedX_ != interpolationNodesOrg) {
            stagedProblem_->stage(interpolationNodesOrg, modes, splineSegments_, std::vector<Real>());
            lastStagedX_ = interpolationNodesOrg;
        }
        
        return stagedProblem_->solve(valuesOrg);
    }

    Eigen::VectorXd BSplineStructure::interpolate(const std::vector<Real>& interpolationNodesOrg,
                                                  const std::vector<Real>& valuesOrg,
                                                  const std::vector<InterpolationMode>& modes,
                                                  const std::vector<Real>& weights) {
        QL_REQUIRE(interpolationNodesOrg.size() == valuesOrg.size(), 
                   "Number of interpolation nodes must match number of values");
        QL_REQUIRE(interpolationNodesOrg.size() == modes.size(),
                   "Number of interpolation nodes must match number of modes");
        if (!weights.empty()) {
            QL_REQUIRE(interpolationNodesOrg.size() == weights.size(),
                       "Number of interpolation nodes must match number of weights");
        }

        // The staged path is the only one that correctly implements mixed HARD/LS modes with weights
        if (!stagedProblem_) {
            stagedProblem_ = ext::make_shared<StagedProblem>(splineConstraints_);
        }
        
        // Check if we need to stage or can reuse
        if (lastStagedX_ != interpolationNodesOrg) {
            stagedProblem_->stage(interpolationNodesOrg, modes, splineSegments_, weights);
            lastStagedX_ = interpolationNodesOrg;
        }
        
        return stagedProblem_->solve(valuesOrg);
    }

    // ReSharper disable once CppInconsistentNaming
    std::vector<Real> BSplineStructure::interpolate_swig(
        const std::vector<Real>& interpolationNodes,
        const std::vector<Real>& values) {
        Eigen::VectorXd result = interpolate(interpolationNodes, values);
        return {result.data(), result.data() + result.size()};
    }

    std::vector<Real> BSplineStructure::interpolate_swig_modes(
        const std::vector<Real>& interpolationNodes,
        const std::vector<Real>& values,
        const std::vector<InterpolationMode>& modes) {
        Eigen::VectorXd result = interpolate(interpolationNodes, values, modes);
        // Create proper copy - the Eigen vector data might not be contiguous or might be deallocated
        std::vector<Real> output;
        output.reserve(result.size());
        for (int i = 0; i < result.size(); ++i) {
            output.push_back(result(i));
        }
        return output;
    }

    std::vector<Real> BSplineStructure::interpolate_swig_weighted(
        const std::vector<Real>& interpolationNodes,
        const std::vector<Real>& values,
        const std::vector<InterpolationMode>& modes,
        const std::vector<Real>& weights) {
        Eigen::VectorXd result = interpolate(interpolationNodes, values, modes, weights);
        // Create proper copy - the Eigen vector data might not be contiguous or might be deallocated
        std::vector<Real> output;
        output.reserve(result.size());
        for (int i = 0; i < result.size(); ++i) {
            output.push_back(result(i));
        }
        return output;
    }
    
    // Integer overloads for SWIG (enum class conversion workaround)
    std::vector<Real> BSplineStructure::interpolate_swig_modes_int(
        const std::vector<Real>& interpolationNodes,
        const std::vector<Real>& values,
        const std::vector<Integer>& modes) {
        // Convert integer modes to enum
        std::vector<InterpolationMode> enumModes;
        enumModes.reserve(modes.size());
        for (Integer mode : modes) {
            enumModes.push_back(static_cast<InterpolationMode>(mode));
        }
        return interpolate_swig_modes(interpolationNodes, values, enumModes);
    }
    
    std::vector<Real> BSplineStructure::interpolate_swig_weighted_int(
        const std::vector<Real>& interpolationNodes,
        const std::vector<Real>& values,
        const std::vector<Integer>& modes,
        const std::vector<Real>& weights) {
        // Convert integer modes to enum
        std::vector<InterpolationMode> enumModes;
        enumModes.reserve(modes.size());
        for (Integer mode : modes) {
            enumModes.push_back(static_cast<InterpolationMode>(mode));
        }
        return interpolate_swig_weighted(interpolationNodes, values, enumModes, weights);
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
        // Update parameters based on mode
        if (!parameters.empty()) {
            if (splineConstraints_->fitData_) {
                // LS mode: parameters affect objective through C matrix
                splineConstraints_->update_c(parameters);
            } else {
                // Hard mode: parameters affect constraints through B matrix
                splineConstraints_->update_b(parameters);
            }
        }
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
        // Delegate to the MultiSegmentBSplineEvaluator which contains all the complex logic
        return spline_.evaluateAll(x, static_cast<BSplineSide>(side));
    }

    Real BSplineStructure::value(const Eigen::VectorXd& coefficients,
                                 Real x,
                                 Integer nu,
                                 BSplineSegment::SideEnum side) const {
        const Size numSegments = this->splineSegments_.size();
        
        // Early exit for empty structure
        if (numSegments == 0) {
            return 0.0;
        }
        
        const Size e = nu >= 0 ? 0 : static_cast<Size>(-nu);
        
        // Get overall range
        const auto firstSegment = this->splineSegments_[0];
        const auto lastSegment = this->splineSegments_[numSegments - 1];
        const Real xMin = firstSegment->range().first;
        const Real xMax = lastSegment->range().second;
        
        // Check if we need extrapolation
        if (x < xMin) {
            // Left extrapolation
            if (nu > 0) {
                // For derivatives, flat extrapolation means derivative = 0
                return 0.0;
            }
            
            // For nu == 0 (value evaluation)
            // Get the value at the left boundary
            Eigen::VectorXd augmentedCoefficients = Eigen::VectorXd::Zero(firstSegment->getNumVariables() + e);
            augmentedCoefficients.segment(0, firstSegment->getNumVariables()) = 
                coefficients.segment(0, firstSegment->getNumVariables());
            
            Real boundaryValue = firstSegment->value(augmentedCoefficients, xMin, 0, BSplineSegment::SideRight);
            
            // Flat extrapolation: f(x) = f(xMin)
            return boundaryValue;
        }
        else if (x > xMax) {
            // Right extrapolation
            if (nu > 0) {
                // For derivatives, flat extrapolation means derivative = 0
                return 0.0;
            }
            
            // For nu == 0 (value evaluation)
            // Get the value at the right boundary
            Size cumVariables = 0;
            for (Size i = 0; i < numSegments - 1; ++i) {
                cumVariables += this->splineSegments_[i]->getNumVariables();
            }
            
            Eigen::VectorXd augmentedCoefficients = Eigen::VectorXd::Zero(lastSegment->getNumVariables() + e);
            augmentedCoefficients.segment(0, lastSegment->getNumVariables()) = 
                coefficients.segment(cumVariables, lastSegment->getNumVariables());
            
            Real boundaryValue = lastSegment->value(augmentedCoefficients, xMax, 0, BSplineSegment::SideLeft);
            
            // Flat extrapolation: f(x) = f(xMax)
            return boundaryValue;
        }
        
        // Original logic for x within domain
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
        // Delegate to the BSplineEvaluator
        return spline_.range();
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