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
                              const std::vector<ext::shared_ptr<BSplineSegment>>& segments) {
        
        QL_REQUIRE(interpolationNodes.size() == modes.size(),
                   "Number of interpolation nodes must match number of modes");
        QL_REQUIRE(!segments.empty(), "Segments cannot be empty");
        
        // Clear previous staging
        clearStaging();
        
        // Create our own MultiSegmentBSplineEvaluator
        evaluator_ = ext::make_shared<MultiSegmentBSplineEvaluator>(segments, baseConstraints_->getNumVariables());
        
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
        
        // Use proper parameter mapping for warm-start optimization
        // This preserves the SCS solver state and enables warm-start
        
        // Build parameter vector from y values
        // Parameters correspond to interpolation points in the order they were added
        std::vector<Real> parameters;
        
        // Add y values for hard interpolation points
        for (Size i = 0; i < hardIndices_.size(); ++i) {
            parameters.push_back(values[hardIndices_[i]]);
        }
        
        // Add y values for soft interpolation points  
        for (Size i = 0; i < softIndices_.size(); ++i) {
            parameters.push_back(values[softIndices_[i]]);
        }
        
        // Update parameters in the constraint system
        // This updates both B matrix (for hard constraints) and C matrix (for LS objectives)
        combinedConstraints_->update_b(parameters);  // Updates RHS: b := b_base + B*parameters
        combinedConstraints_->update_c(parameters);  // Updates objective: c := c_base + C*parameters
        
        // Solve using updated parameters - this enables warm-start
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
                
                // Evaluate basis functions at this point using MultiSegmentBSplineEvaluator
                Eigen::VectorXd basis = evaluator_->evaluateAll(x, SideRight);
                
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
                Eigen::VectorXd basis = evaluator_->evaluateAll(x, SideRight);
                
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
        
        // Create combined quadratic form that includes soft constraints
        std::vector<std::vector<Real>> P_combined = P_base;
        
        // Add soft constraint quadratic form to P matrix
        if (!softIndices_.empty() && Q_soft_.rows() > 0) {
            // Q_soft_ contains the sum of outer products: sum_i (b_i * b_i')
            // where b_i is the basis vector at soft point i
            // Add this to P matrix for least squares fitting
            for (Size i = 0; i < numVariables; ++i) {
                if (i >= P_combined.size()) {
                    P_combined.resize(i + 1);
                }
                if (P_combined[i].size() < numVariables) {
                    P_combined[i].resize(numVariables, 0.0);
                }
                
                for (Size j = 0; j < numVariables; ++j) {
                    P_combined[i][j] += Q_soft_.coeff(i, j);
                }
            }
        }
        
        // Create combined constraints with parameter support
        combinedConstraints_ = ext::make_shared<SplineConstraints>(
            numVariables, P_combined, A_combined, b_combined
        );
        
        // Add parameter mapping for warm-start capability
        Size totalParameters = hardIndices_.size() + softIndices_.size();
        if (totalParameters > 0) {
            // Build B matrix for hard constraints (maps parameters to RHS)
            Eigen::SparseMatrix<Real> B_matrix(A_combined.size(), totalParameters);
            std::vector<Eigen::Triplet<Real>> B_triplets;
            
            // Hard constraints: parameter i affects constraint (baseRows + i) with coefficient 1.0
            for (Size i = 0; i < hardIndices_.size(); ++i) {
                B_triplets.emplace_back(baseRows + i, i, 1.0);
            }
            B_matrix.setFromTriplets(B_triplets.begin(), B_triplets.end());
            
            // Build C matrix for LS objective (maps parameters to objective linear term)
            Eigen::SparseMatrix<Real> C_matrix(numVariables, totalParameters);
            std::vector<Eigen::Triplet<Real>> C_triplets;
            
            // For soft points, add -A'*y contribution to objective
            for (Size i = 0; i < softIndices_.size(); ++i) {
                Size softIdx = softIndices_[i];
                Real x = lastX_[softIdx];
                Eigen::VectorXd basis = evaluator_->evaluateAll(x, SideRight);
                
                for (Size j = 0; j < basis.size(); ++j) {
                    if (std::abs(basis[j]) > 1e-14) {
                        // Map soft parameter (hardIndices_.size() + i) to variable j
                        C_triplets.emplace_back(j, hardIndices_.size() + i, -basis[j]);
                    }
                }
            }
            C_matrix.setFromTriplets(C_triplets.begin(), C_triplets.end());
            
            // Add parameters to the constraint system
            combinedConstraints_->addParameters(totalParameters, B_matrix, C_matrix);
        }
    }

    std::vector<ModeSpan> generateVoronoiModeSpans(
        const std::vector<Real>& dataPoints,
        Real domainStart,
        Real domainEnd,
        Real densityThreshold,
        Real boundaryWidth) {
        
        QL_REQUIRE(domainStart < domainEnd, "Domain start must be less than domain end");
        QL_REQUIRE(densityThreshold > 0, "Density threshold must be positive");
        QL_REQUIRE(boundaryWidth >= 0, "Boundary width must be non-negative");
        
        std::vector<ModeSpan> spans;
        
        // Sort data points for processing
        std::vector<Real> sortedPoints = dataPoints;
        std::sort(sortedPoints.begin(), sortedPoints.end());
        
        // Remove duplicates and filter to domain
        auto last = std::unique(sortedPoints.begin(), sortedPoints.end());
        sortedPoints.erase(last, sortedPoints.end());
        
        // Filter to domain range
        auto domainBegin = std::lower_bound(sortedPoints.begin(), sortedPoints.end(), domainStart);
        auto domainEnd_it = std::upper_bound(sortedPoints.begin(), sortedPoints.end(), domainEnd);
        std::vector<Real> domainPoints(domainBegin, domainEnd_it);
        
        if (domainPoints.empty()) {
            // No data points in domain - use HARD mode for entire domain
            spans.emplace_back(domainStart, domainEnd, InterpolationMode::HARD);
            return spans;
        }
        
        // Define boundary regions (always HARD mode for boundary stability)
        Real leftBoundaryEnd = domainStart + boundaryWidth;
        Real rightBoundaryStart = domainEnd - boundaryWidth;
        
        // Process the domain in segments
        Real currentPos = domainStart;
        const Real minSpanWidth = 0.01; // Minimum span width to avoid tiny spans
        
        while (currentPos < domainEnd) {
            Real spanStart = currentPos;
            Real spanEnd;
            InterpolationMode spanMode;
            
            // Determine if we're in a boundary region
            bool inLeftBoundary = (currentPos < leftBoundaryEnd);
            bool inRightBoundary = (currentPos >= rightBoundaryStart);
            
            if (inLeftBoundary) {
                // Left boundary region - use HARD mode
                spanEnd = std::min(leftBoundaryEnd, domainEnd);
                spanMode = InterpolationMode::HARD;
            } else if (inRightBoundary) {
                // Right boundary region - use HARD mode
                spanEnd = domainEnd;
                spanMode = InterpolationMode::HARD;
            } else {
                // Interior region - analyze local density
                
                // Find data points in a local window around current position
                Real windowSize = std::min(0.2, (domainEnd - domainStart) / 5.0); // Adaptive window
                Real windowStart = currentPos;
                Real windowEnd = std::min(currentPos + windowSize, rightBoundaryStart);
                
                // Count points in window
                auto windowBegin = std::lower_bound(domainPoints.begin(), domainPoints.end(), windowStart);
                auto windowEnd_it = std::upper_bound(domainPoints.begin(), domainPoints.end(), windowEnd);
                Size pointCount = std::distance(windowBegin, windowEnd_it);
                
                // Calculate local density
                Real windowWidth = windowEnd - windowStart;
                Real localDensity = (windowWidth > 1e-10) ? pointCount / windowWidth : 0.0;
                
                // Determine mode based on density
                spanMode = (localDensity >= densityThreshold) ? InterpolationMode::LS : InterpolationMode::HARD;
                
                // Extend span to include similar-density regions
                spanEnd = windowEnd;
                Real extendStep = windowSize * 0.5;
                
                while (spanEnd < rightBoundaryStart && (spanEnd - spanStart) < (domainEnd - domainStart) * 0.5) {
                    Real nextWindowEnd = std::min(spanEnd + extendStep, rightBoundaryStart);
                    
                    // Check density in extended region
                    auto extWindowBegin = std::lower_bound(domainPoints.begin(), domainPoints.end(), spanEnd);
                    auto extWindowEnd_it = std::upper_bound(domainPoints.begin(), domainPoints.end(), nextWindowEnd);
                    Size extPointCount = std::distance(extWindowBegin, extWindowEnd_it);
                    
                    Real extWindowWidth = nextWindowEnd - spanEnd;
                    Real extDensity = (extWindowWidth > 1e-10) ? extPointCount / extWindowWidth : 0.0;
                    
                    // Check if extended region has similar density characteristics
                    bool extShouldLS = (extDensity >= densityThreshold);
                    bool currentShouldLS = (spanMode == InterpolationMode::LS);
                    
                    if (extShouldLS == currentShouldLS) {
                        // Similar density - extend the span
                        spanEnd = nextWindowEnd;
                    } else {
                        // Different density - end current span
                        break;
                    }
                }
            }
            
            // Ensure minimum span width
            if (spanEnd - spanStart < minSpanWidth && spanEnd < domainEnd) {
                spanEnd = std::min(spanStart + minSpanWidth, domainEnd);
            }
            
            // Add span if it has meaningful width
            if (spanEnd > spanStart + 1e-12) {
                spans.emplace_back(spanStart, spanEnd, spanMode);
            }
            
            currentPos = spanEnd;
        }
        
        // Merge adjacent spans with same mode to avoid fragmentation
        if (spans.size() > 1) {
            std::vector<ModeSpan> mergedSpans;
            mergedSpans.push_back(spans[0]);
            
            for (Size i = 1; i < spans.size(); ++i) {
                if (spans[i].mode == mergedSpans.back().mode && 
                    std::abs(spans[i].start - mergedSpans.back().end) < 1e-10) {
                    // Same mode and adjacent - merge
                    mergedSpans.back().end = spans[i].end;
                } else {
                    // Different mode or gap - add as new span
                    mergedSpans.push_back(spans[i]);
                }
            }
            spans = std::move(mergedSpans);
        }
        
        return spans;
    }


}