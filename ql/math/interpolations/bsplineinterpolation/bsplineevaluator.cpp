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

#include "bsplineevaluator.hpp"
#include <stdexcept>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <ql/errors.hpp>
#include <ql/types.hpp>

//  TODO Lets ReSharper relax about the type mismatch between QuantLib::Size and Eigen::Index, should revisit this
//  ReSharper disable CppClangTidyBugproneNarrowingConversions

namespace QuantLib {
    BSplineEvaluator::BSplineEvaluator() : knots_({0.0, 1.0}), degree_(0), numBasisFunctions_(1) {}

    BSplineEvaluator::BSplineEvaluator(const std::vector<double>& knots, Size degree)
    : knots_(knots), degree_(degree), numBasisFunctions_(knots.size() - degree - 1) {
        precomputeRkMatrices();
        initializeTempVectors();
    }

    void BSplineEvaluator::precomputeRkMatrices() const {
        Rk_matrices_.resize(degree_ + 1);
        for (Size k = 1; k <= degree_; ++k) {
            Rk_matrices_[k] = Eigen::SparseMatrix<double>(k + 1, k);
            Rk_matrices_[k].reserve(Eigen::VectorXi::Constant(k, 2)); // Each column will have at most 2 non-zero entries
        }
    }

    void BSplineEvaluator::initializeTempVectors() const {
        tempB1_.resize(degree_ + 1);
        tempB2_.resize(degree_ + 1);
    }

    // Pre: knots have at least d+1 repeats of the first and last entries
    // Pre: knots_.begin() <= x <= knots_.end()
    // knots_ are t_0, t_1, ..., t_{p+d+1} where p is the number of basis functions
    // lower_bound returns the first element that is not less than x, which is then t_{\mu+1}
    // in particular \mu+1>=d+1

    Size BSplineEvaluator::findKnotSpan(double x) const {
        // TODO need to rethink this a bit, should maybe return the last one in the second case
        // This returns the index of the first knot that is greater than x
        auto it = std::upper_bound(knots_.begin(), knots_.end(), x);  // Should be binary search
        if (it != knots_.end()) {
            return std::distance(knots_.begin(), it);
        } else if (x >= knots_.back()) {
            return knots_.size() - degree_ - 1; // TODO Not -2, right?
        } else {
            //TODO This never happens, right?
            QL_FAIL("x = " << x << " is outside the range of the knots [" << knots_.front() << ", "
                           << knots_.back() << "]");
        }
    }

    Eigen::VectorXd BSplineEvaluator::evaluateAll(double x) const {
        Size mu = findKnotSpan(x); // So mu >= degree + 1
        Eigen::VectorXd B = Eigen::VectorXd::Zero(numBasisFunctions_);

        evaluate(B.segment(mu - degree_ - 1, degree_ + 1), x, mu);

        return B;
    }

    /*
     * This implementation is inspired from Lyche-Morken "Spline Methods" book, chapter 2.4
     */
    void BSplineEvaluator::evaluate(Eigen::Ref<Eigen::VectorXd> B, double x, Size mu) const {
        Size d = degree_;
        //Eigen::VectorXd* tempB1 = &tempB1_;
        //Eigen::VectorXd* tempB2 = &tempB2_;

        // Initialize B_0
        //tempB1->head(1).coeffRef(0) = 1.0;
        tempB1_.head(1).coeffRef(0) = 1.0;

        // Compute B_k for k = 1, ..., d
        for (Size k = 1; k <= d; ++k) {
            auto& Rk = Rk_matrices_[k];
            for (Size i = 0; i < k; ++i) {
                if (knots_[i + mu - k] != knots_[i + mu]) {
                    Rk.coeffRef(i, i) =
                        (knots_[i + mu] - x) / (knots_[i + mu] - knots_[i + mu - k]);
                    Rk.coeffRef(i + 1, i) =
                        (x - knots_[i + mu - k]) / (knots_[i + mu] - knots_[i + mu - k]);
                }
            }
            //tempB2->head(k + 1) = Rk * tempB1->head(k);
            // The noalias is needed to avoid aliasing issues, meaning that the result is stored
            // in a temporary result and then copied to the destination
            tempB2_.head(k + 1).noalias() = Rk * tempB1_.head(k);
            tempB1_.swap(tempB2_);
            //std::swap(tempB1, tempB2);

        }
        B.head(d + 1) = tempB1_.head(d + 1);
    }

    // Pre: x is within domain of the spline
    double BSplineEvaluator::value(const Eigen::VectorXd& coefficients, double x) const {
        QL_REQUIRE(coefficients.size() == numBasisFunctions_,
                   "The size of coefficients vector must match the number of basis functions.");

        Size mu = findKnotSpan(x);
        evaluate(tempB1_, x, mu);

        // Compute the inner product with only the relevant coefficients
        return tempB1_.dot(coefficients.segment(mu - degree_ - 1, degree_ + 1));
    }
} // namespace QuantLib
