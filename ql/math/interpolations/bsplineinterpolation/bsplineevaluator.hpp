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

#ifndef bsplineevaluator_hpp
#define bsplineevaluator_hpp


#include <Eigen/Sparse>
#include <vector>
#include <ql/types.hpp>

namespace QuantLib {
    class BSplineEvaluator {
      public:
        BSplineEvaluator();

        BSplineEvaluator(const std::vector<double>& knots, Size degree);

        // Public method to evaluate all basis functions at x
        Eigen::VectorXd evaluateAll(double x) const;

        // Public method to compute the value of the spline given coefficients
        double value(const Eigen::VectorXd& coeffs, double x) const;


      private:
        std::vector<double> knots_;
        Size degree_;
        Size numBasisFunctions_;
        mutable std::vector<Eigen::SparseMatrix<double>> Rk_matrices_;
        mutable Eigen::VectorXd tempB1_;
        mutable Eigen::VectorXd tempB2_;

        Size findKnotSpan(double x) const;
        void precomputeRkMatrices();
        void initializeTempVectors();

        // Private method to evaluate non-zero basis functions at x
        void evaluate(Eigen::Ref<Eigen::VectorXd> B, double x, Size mu) const;
    };
}

#endif bsplineevaluator_hpp
