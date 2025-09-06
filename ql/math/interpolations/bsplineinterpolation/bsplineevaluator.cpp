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

    BSplineEvaluator::BSplineEvaluator(const std::vector<Real>& knots, Integer degree)
        : knots_(knots), degree_(degree) {
        QL_REQUIRE(!knots_.empty(), "Knots vector cannot be empty");
        QL_REQUIRE(degree_ >= 0, "Degree must be non-negative");
    }

    Eigen::VectorXd BSplineEvaluator::evaluateAll(Real x) const {
        // This is a placeholder implementation
        // In practice, this should implement the Cox-de Boor recursion formula
        // or use another method to evaluate B-spline basis functions
        
        Size n = knots_.size() - degree_ - 1;
        if (n <= 0) {
            return Eigen::VectorXd::Zero(1);
        }
        
        Eigen::VectorXd result = Eigen::VectorXd::Zero(n);
        
        // Simple placeholder: find the interval and set one basis function to 1
        // This is NOT a correct B-spline evaluation, just enough to compile
        for (Size i = 0; i < knots_.size() - 1; ++i) {
            if (x >= knots_[i] && x < knots_[i + 1] && i < n) {
                result[i] = 1.0;
                break;
            }
        }
        
        return result;
    }

    Real BSplineEvaluator::value(const Eigen::VectorXd& coefficients, Real x) const {
        Eigen::VectorXd basis = evaluateAll(x);
        if (basis.size() != coefficients.size()) {
            return 0.0;
        }
        return basis.dot(coefficients);
    }

}