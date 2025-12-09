/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2024-2025 SEB AB Sverrir Thorvaldsson

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.
*/

#include "bsplineevaluator.hpp"
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <ql/errors.hpp>
#include <ql/types.hpp>

namespace QuantLib {
    
    BSplineEvaluator::BSplineEvaluator() : knots_({0.0, 1.0}), degree_(0), numBasisFunctions_(1) {}

    BSplineEvaluator::BSplineEvaluator(const std::vector<Real>& knots, Integer degree)
    : knots_(knots), degree_(degree), numBasisFunctions_(knots.size() - degree - 1) {
        QL_REQUIRE(!knots_.empty(), "Knots vector cannot be empty");
        QL_REQUIRE(degree_ >= 0, "Degree must be non-negative");
        precomputeRkMatrices();
        initializeTempVectors();
    }

    void BSplineEvaluator::precomputeRkMatrices() const {
        Rk_matrices_.resize(degree_ + 1);
        for (Size k = 1; k <= degree_; ++k) {
            Rk_matrices_[k] = Eigen::SparseMatrix<Real>(k + 1, k);
            Rk_matrices_[k].reserve(Eigen::VectorXi::Constant(k, 2)); // Each column will have at most 2 non-zero entries
        }
    }

    void BSplineEvaluator::initializeTempVectors() const {
        tempB1_.resize(degree_ + 1);
        tempB2_.resize(degree_ + 1);
    }

    Size BSplineEvaluator::findKnotSpan(Real x) const {
        // This returns the index of the first knot that is greater than x
        auto it = std::upper_bound(knots_.begin(), knots_.end(), x);
        if (it != knots_.end()) {
            return std::distance(knots_.begin(), it);
        } else if (x >= knots_.back()) {
            return knots_.size() - degree_ - 1;
        } else {
            QL_FAIL("x = " << x << " is outside the range of the knots [" << knots_.front() << ", "
                           << knots_.back() << "]");
        }
    }

    Eigen::VectorXd BSplineEvaluator::evaluateAll(Real x) const {
        Size mu = findKnotSpan(x);
        Eigen::VectorXd B = Eigen::VectorXd::Zero(numBasisFunctions_);
        evaluate(B.segment(mu - degree_ - 1, degree_ + 1), x, mu);
        return B;
    }

    /*
     * This implementation is inspired from Lyche-Morken "Spline Methods" book, chapter 2.4
     * This is the efficient Cox-de Boor recursion using pre-allocated matrices
     */
    void BSplineEvaluator::evaluate(Eigen::Ref<Eigen::VectorXd> B, Real x, Size mu) const {
        Size d = degree_;

        // Initialize B_0
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
            // The noalias is needed to avoid aliasing issues
            tempB2_.head(k + 1).noalias() = Rk * tempB1_.head(k);
            tempB1_.swap(tempB2_);
        }
        B.head(d + 1) = tempB1_.head(d + 1);
    }

    Real BSplineEvaluator::value(const Eigen::VectorXd& coefficients, Real x) const {
        QL_REQUIRE(static_cast<Size>(coefficients.size()) == this->numBasisFunctions_,
                   "The size of coefficients vector must match the number of basis functions.");

        const Size mu = findKnotSpan(x);
        evaluate(tempB1_, x, mu);

        // Compute the inner product with only the relevant coefficients
        return tempB1_.dot(coefficients.segment(mu - this->degree_ - 1, this->degree_ + 1));
    }

}