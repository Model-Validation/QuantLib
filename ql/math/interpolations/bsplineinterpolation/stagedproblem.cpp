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

#include "stagedproblem.hpp"
#include "bsplineevaluator.hpp"
#include <ql/errors.hpp>
#include <algorithm>
#include <numeric>

namespace QuantLib {

    StagedProblem::StagedProblem(const ext::shared_ptr<SplineConstraints>& baseConstraints)
        : baseConstraints_(baseConstraints),
          baseNumConstraints_(baseConstraints->getNumEqualities() + baseConstraints->getNumInequalities()),
          baseNumEqualities_(baseConstraints->getNumEqualities()),
          baseNumInequalities_(baseConstraints->getNumInequalities()),
          staged_(false) {
        
        QL_REQUIRE(baseConstraints_, "Base constraints cannot be null");
    }

    void StagedProblem::stage(const std::vector<Real>& interpolationNodes,
                              const std::vector<InterpolationMode>& modes,
                              const std::vector<ext::shared_ptr<BSplineSegment>>& segments,
                              std::function<Eigen::VectorXd(Real, BSplineSegment::SideEnum)> evaluator) {
        
        QL_REQUIRE(interpolationNodes.size() == modes.size(),
                   "Number of interpolation nodes must match number of modes");
        QL_REQUIRE(!segments.empty(), "Segments cannot be empty");
        
        // Clear previous staging
        clearStaging();
        
        // Store the evaluator if provided
        evaluator_ = evaluator;
        
        // Build interpolation matrices
        buildInterpolationMatrices(interpolationNodes, modes, segments);
        
        // Combine base and interpolation constraints
        combineConstraints();
        
        // Mark as staged and remember x points and segments
        staged_ = true;
        lastX_ = interpolationNodes;
        segments_ = segments;
    }

    void StagedProblem::stage(const std::vector<Real>& interpolationNodes,
                              const std::vector<ModeSpan>& modeSpans,
                              const std::vector<ext::shared_ptr<BSplineSegment>>& segments) {
        
        // Convert spans to per-point modes
        std::vector<InterpolationMode> modes;
        modes.reserve(interpolationNodes.size());
        
        for (Real x : interpolationNodes) {
            modes.push_back(getModeAt(x, modeSpans));
        }
        
        // Call the main stage method
        stage(interpolationNodes, modes, segments);
    }

    Eigen::VectorXd StagedProblem::solve(const std::vector<Real>& values) {
        QL_REQUIRE(staged_, "Problem must be staged before solving");
        QL_REQUIRE(values.size() == lastX_.size(),
                   "Number of y values must match staged x points");
        
        // Build the complete b vector - must match A matrix rows
        // In our case, base has no constraints initially, so b_base is empty
        // We only have hard interpolation constraints
        std::vector<Real> b_current;
        auto b_base = baseConstraints_->get_b_vector();
        
        // Copy base b values (usually empty for our test)
        b_current = b_base;
        
        // Add the y values for hard interpolation constraints
        for (Size i = 0; i < hardIndices_.size(); ++i) {
            b_current.push_back(values[hardIndices_[i]]);
        }
        
        // For now, rebuild constraints with updated b vector
        // (This is not ideal but works until we have dynamic parameter support)
        auto P_matrix = combinedConstraints_->get_p_matrix();
        auto A_matrix = combinedConstraints_->get_a_matrix();
        
        combinedConstraints_ = ext::make_shared<SplineConstraints>(
            baseConstraints_->getNumVariables(),
            P_matrix, A_matrix, b_current
        );
        
        // Solve using the combined constraints
        combinedConstraints_->solve();
        return combinedConstraints_->getSolution();
    }

    void StagedProblem::clearStaging() {
        staged_ = false;
        lastX_.clear();
        segments_.clear();
        hardIndices_.clear();
        softIndices_.clear();
        A_hard_ = Eigen::SparseMatrix<Real>();
        Q_soft_ = Eigen::SparseMatrix<Real>();
        combinedConstraints_ = nullptr;
        evaluator_ = nullptr;
    }

    InterpolationMode StagedProblem::getModeAt(Real x, 
                                               const std::vector<ModeSpan>& spans) const {
        // Find the span containing x
        for (const auto& span : spans) {
            if (x >= span.start && x <= span.end) {
                return span.mode;
            }
        }
        
        // Default to AUTO if not in any span
        return InterpolationMode::AUTO;
    }

    void StagedProblem::buildInterpolationMatrices(
            const std::vector<Real>& interpolationNodes,
            const std::vector<InterpolationMode>& modes,
            const std::vector<ext::shared_ptr<BSplineSegment>>& segments) {
        
        Size n = interpolationNodes.size();
        Size numVariables = baseConstraints_->getNumVariables();
        
        // Separate hard and soft points
        hardIndices_.clear();
        softIndices_.clear();
        
        for (Size i = 0; i < n; ++i) {
            if (modes[i] == InterpolationMode::HARD) {
                hardIndices_.push_back(i);
            } else if (modes[i] == InterpolationMode::LS) {
                softIndices_.push_back(i);
            } else {
                // AUTO mode: decide based on DOF
                Size totalConstraints = baseNumConstraints_ + hardIndices_.size();
                if (totalConstraints < numVariables) {
                    // Have DOF, make it hard
                    hardIndices_.push_back(i);
                } else {
                    // No DOF, make it soft
                    softIndices_.push_back(i);
                }
            }
        }
        
        // Build hard constraint matrix A_hard
        if (!hardIndices_.empty()) {
            std::vector<Eigen::Triplet<Real>> triplets;
            
            for (Size row = 0; row < hardIndices_.size(); ++row) {
                Size pointIdx = hardIndices_[row];
                Real x = interpolationNodes[pointIdx];
                
                // Evaluate basis functions at this point
                Eigen::VectorXd basis;
                if (evaluator_) {
                    // Use the provided evaluator (from BSplineStructure)
                    basis = evaluator_(x, BSplineSegment::SideRight);
                } else {
                    // Fall back to our own evaluation
                    basis = evaluateBasisAt(x, segments);
                }
                
                // Add to triplets
                for (Size col = 0; col < basis.size(); ++col) {
                    if (std::abs(basis[col]) > 1e-14) {
                        triplets.emplace_back(row, col, basis[col]);
                    }
                }
            }
            
            A_hard_.resize(hardIndices_.size(), numVariables);
            A_hard_.setFromTriplets(triplets.begin(), triplets.end());
            
            // Debug: verify size
            QL_REQUIRE(A_hard_.rows() == hardIndices_.size(),
                       "A_hard rows mismatch: " << A_hard_.rows() << " != " << hardIndices_.size());
        }
        
        // Build soft quadratic form Q_soft
        if (!softIndices_.empty()) {
            // Q_soft = sum_i (b_i * b_i'), where b_i is basis at soft point i
            Eigen::SparseMatrix<Real> Q_accumulator(numVariables, numVariables);
            
            for (Size i : softIndices_) {
                Real x = interpolationNodes[i];
                Eigen::VectorXd basis;
                if (evaluator_) {
                    basis = evaluator_(x, BSplineSegment::SideRight);
                } else {
                    basis = evaluateBasisAt(x, segments);
                }
                
                // Add outer product to Q
                for (Size row = 0; row < basis.size(); ++row) {
                    for (Size col = 0; col < basis.size(); ++col) {
                        Q_accumulator.coeffRef(row, col) += basis[row] * basis[col];
                    }
                }
            }
            
            Q_soft_ = Q_accumulator;
        }
    }

    void StagedProblem::combineConstraints() {
        // Create a new SplineConstraints that combines base + interpolation
        
        Size numVariables = baseConstraints_->getNumVariables();
        Size numHard = hardIndices_.size();
        
        // For now, create a simple combined constraint system
        // We'll enhance this when SplineConstraints API is updated
        
        // Get base constraint matrices
        auto P_base = baseConstraints_->get_p_matrix();
        auto A_base = baseConstraints_->get_a_matrix();
        auto b_base = baseConstraints_->get_b_vector();
        
        // Combine A matrices (base + hard interpolation constraints)
        std::vector<std::vector<Real>> A_combined = A_base;
        
        // Debug: Check initial sizes
        Size baseRows = A_base.size();
        
        if (numHard > 0) {
            // Add rows for hard interpolation constraints
            QL_REQUIRE(A_hard_.rows() == numHard, 
                       "A_hard rows " << A_hard_.rows() << " != numHard " << numHard);
            
            for (Size i = 0; i < numHard; ++i) {
                std::vector<Real> row(numVariables, 0.0);
                for (Size j = 0; j < numVariables; ++j) {
                    row[j] = A_hard_.coeff(i, j);
                }
                A_combined.push_back(row);
            }
        }
        
        // Debug: verify final size
        QL_REQUIRE(A_combined.size() == baseRows + numHard,
                   "A_combined size wrong: " << A_combined.size() << 
                   " != " << baseRows << " + " << numHard);
        
        // Extend b vector (will be filled with y values when solving)
        std::vector<Real> b_combined = b_base;
        for (Size i = 0; i < numHard; ++i) {
            b_combined.push_back(0.0);  // Placeholder, will be set in solve()
        }
        
        // Create combined constraints
        // Note: For now, we ignore soft constraints (Q_soft_)
        // This will be added when SplineConstraints supports dynamic objectives
        combinedConstraints_ = ext::make_shared<SplineConstraints>(
            numVariables, P_base, A_combined, b_combined
        );
    }

    Eigen::VectorXd StagedProblem::evaluateBasisAt(
            Real x,
            const std::vector<ext::shared_ptr<BSplineSegment>>& segments,
            BSplineSegment::SideEnum side) const {
        
        // Find the segment containing x
        for (const auto& segment : segments) {
            auto range = segment->range();
            if (x >= range.first && x <= range.second) {
                // Determine sidedness for boundary evaluation
                if (std::abs(x - range.second) < 1e-10) {
                    // Right boundary: use LEFT-sided evaluation
                    return segment->evaluateAll(x, BSplineSegment::SideLeft);
                } else if (std::abs(x - range.first) < 1e-10) {
                    // Left boundary: use RIGHT-sided evaluation
                    return segment->evaluateAll(x, BSplineSegment::SideRight);
                } else {
                    // Interior: use specified side (default RIGHT)
                    return segment->evaluateAll(x, side);
                }
            }
        }
        
        QL_FAIL("Point " << x << " not found in any segment");
    }

}