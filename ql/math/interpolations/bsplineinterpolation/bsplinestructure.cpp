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

// TODO Lets ReSharper relax about the type mismatch between QuantLib::Size and Eigen::Index, should revisit this
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
        const ext::shared_ptr<SplineConstraints>& splineConstraints) :
        splineSegments_(splineSegments),
        splineConstraints_(splineConstraints) {
        std::cout << "Structure with " << splineSegments.size() << " segments" << std::endl;
    }

    void BSplineStructure::addInterpolationNodes(const std::vector<Real>& interpolationNodes) const {
        const Size nInterpolationNodes = interpolationNodes.size();
        const Size nConstraints = splineConstraints_->getNConstraints();
        const Size nVariables = splineConstraints_->getNumVariables();

        for (const Real interpolationNode : interpolationNodes) {
            Eigen::VectorXd row = evaluateAll(interpolationNode);
            splineConstraints_->addLinearConstraint(row, 0.0, SplineConstraints::ConstraintType::Equal);
        }

        // We still have the old constraint count
        Eigen::SparseMatrix<Real> B_new(nConstraints + nInterpolationNodes, nInterpolationNodes);
        for (Size i = 0; i < nInterpolationNodes; ++i) {
            B_new.insert(nConstraints + i, i) = 1.0;
        }

        Eigen::SparseMatrix<double> C_new(nVariables, nInterpolationNodes);

        splineConstraints_->addParameters(nInterpolationNodes, B_new, C_new);
    }

    Eigen::VectorXd
    BSplineStructure::interpolate(const std::vector<Real>& interpolationNodes,
                                  const std::vector<Real>& values) const {
        splineConstraints_->push();
        addInterpolationNodes(interpolationNodes);
        Eigen::VectorXd solution = solve(transform(interpolationNodes, values));
        splineConstraints_->pop();
        return solution;
    }

    Eigen::VectorXd BSplineStructure::solve(const std::vector<Real>& parameters) const {
        splineConstraints_->update_b(parameters);
        int status = splineConstraints_->solve();
        QL_REQUIRE(status == 1, "Solution failed, returned " << status << '\n');
        return splineConstraints_->getSolution();
    }

    Eigen::VectorXd BSplineStructure::evaluateAll(Real x) const {
        Size nSegments = splineSegments_.size();
        Size nVariables = splineConstraints_->getNumVariables();
        Size j = 0;
        Eigen::VectorXd result = Eigen::VectorXd::Zero(nVariables);

        // TODO endpoints are not handled correctly, use side variables
        for (Size i = 0; i < nSegments; ++i) {
            if (splineSegments_[i]->range().first <= x && x < splineSegments_[i]->range().second) {
                Eigen::VectorXd segmentResult = splineSegments_[i]->evaluateAll(x);
                result.segment(j, splineSegments_[i]->getNumVariables()) = segmentResult;
            }
            j += splineSegments_[i]->getNumVariables();
        }

        return result;
    }

    Real BSplineStructure::value(const Eigen::VectorXd& coefficients,
                                            Real x,
                                            BSplineSegment::Side side) const {
        Size numSegments = splineSegments_.size();
        Size mu = 0;
        Size cumVariables = 0;
        for (Size i = 0; i < numSegments; ++i) {
            if (splineSegments_[i]->range().first <= x && x < splineSegments_[i]->range().second) {
                mu = i;
                break;
            }
            cumVariables += splineSegments_[i]->getNumVariables();
        }
        return splineSegments_[mu]->value(
            coefficients.segment(cumVariables, splineSegments_[mu]->getNumVariables()), x, side);
    }

    std::pair<Real, Real> BSplineStructure::range() const {
        Real left = splineSegments_.front()->range().first;
        Real right = splineSegments_.back()->range().second;

        return {left, right};
    }

    void BSplineStructure::setConstraints(const ext::shared_ptr<SplineConstraints>& splineConstraints) {
        splineConstraints_ = splineConstraints;
    }

    std::vector<Real> BSplineStructure::transform(const std::vector<Real>& abscissae,
                                                  const std::vector<Real>& values) const {

        // Check if abscissae is sorted and vectors are of same length
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
                // TODO endpoints are not handled correctly, use side variables
                if ((x < left && segmentIndex == 0) || (x >= left && x < right) ||
                    (segmentIndex == numSegments - 1 && x >= right)) {
                    // Apply the transformation
                    transformedValues[i] = splineSegments_[segmentIndex]->transform(x, value);
                    break;
                }
                // TODO More careful here, not sure this is safe without further checks, e.g. if x < left
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

}
