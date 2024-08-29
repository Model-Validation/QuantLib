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
#include "splineconstraints.hpp"
#include "splinestructure.hpp"
#include <algorithm>
#include <cmath>
#include <iosfwd>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>
#include <ql/errors.hpp>
#include <ql/types.hpp>
#include <ql/shared_ptr.hpp>

namespace QuantLib {
    /*
    Class for a spline structure holding the necessary information on knots, degree and constraints
    for a spline.
    */
    BSplineStructure::BSplineStructure(const std::vector<double>& simpleKnots,
                                     Size degree,
                                     const std::vector<Integer>& knotIndices,
                                     ext::shared_ptr<SplineConstraints>& splineConstraints,
                                     InterpolationSmoothness smoothness,
                                     Side side,
                                     //const std::vector<double>& interpolationNodes,
                                     //double knotTolerance,
                                     //double rankTolerance,
                                     int requiredPoints,
                                     bool isGlobal)
    : simpleKnots_(simpleKnots), degree_(degree), nSimpleKnots_(simpleKnots.size()),
      knotIndices_(knotIndices),
      splineConstraints_(std::move(splineConstraints)), side_(side), 
        interpolationSmoothness_(smoothness),
      //knotTolerance_(knotTolerance), rankTolerance_(rankTolerance),
      requiredPoints_(requiredPoints), isGlobal_(isGlobal),
      startPoint_(simpleKnots.front()),
      endPoint_(simpleKnots.back()) {

        QL_REQUIRE(std::is_sorted(simpleKnots_.begin(), simpleKnots_.end(), std::less<>()),
                   "The x (simpleKnots) must be strictly increasing.");

        if (knotIndices_.empty()) {
            defaultKnotIndices();
        } else {
            QL_REQUIRE(
                checkKnotIndices(knotIndices_, degree_, nSimpleKnots_ - 1),
                "The provided knots are invalid for the spline of the given degree.\n"
                "The sequence must be non-decreasing and cover all indexes, no index repeated "
                "more than degree+1 times, and the\n"
                "end knots must have exactly that count.");
        }

        nKnots_ = knotIndices_.size();

        knots_.resize(nKnots_);
        for (Size i = 0; i < nKnots_; ++i) {
            knots_[i] = simpleKnots_[knotIndices_[i]];
        }

        spline_ = BSplineEvaluator(knots_, degree_);
    }

    // Copy constructor
    BSplineStructure::BSplineStructure(const BSplineStructure& other) {
        *this = other;
        splineConstraints_ = ext::make_shared<SplineConstraints>(*other.splineConstraints_);
        spline_ = BSplineEvaluator(knots_, degree_);
    }


    // Accessor for knot range
    const std::pair<double, double> BSplineStructure::range() const {
        return {startPoint_, endPoint_};
    }


    // Accessor for knots_
    const std::vector<double>& BSplineStructure::knots() const {
        return knots_;
    }

    // Accessor for degree_
    Size BSplineStructure::degree() const {
        return degree_;
    }

    void BSplineStructure::setConstraints(ext::shared_ptr<SplineConstraints>& splineConstraints) {
        splineConstraints_ = splineConstraints;
    }

    void BSplineStructure::addInterpolationNodes(const std::vector<double>& interpolationNodes) {
        Size nInterpolationNodes = interpolationNodes.size();
        Size nConstraints = splineConstraints_->getNConstraints();
        Size nParameters = splineConstraints_->getNParameters();
        Size nVariables = splineConstraints_->getNVariables();
        interpolationNodes_ = interpolationNodes;

        for (Size i = 0; i < interpolationNodes.size(); ++i) {
            Eigen::VectorXd row = evaluateAll(interpolationNodes[i], degree_, side_);
            splineConstraints_->addLinearConstraint(row, 0.0);
        }

        // We still have the old constraint count
        Eigen::SparseMatrix<double> B_new(nConstraints + nInterpolationNodes, nInterpolationNodes);
        for (int i = 0; i < nInterpolationNodes; ++i) {
            B_new.insert(nConstraints + i, i) = 1.0;
        }

        Eigen::SparseMatrix<double> C_new(nVariables, nInterpolationNodes);

        splineConstraints_->addParameters(nInterpolationNodes, B_new, C_new);
    }

    Eigen::VectorXd BSplineStructure::interpolate(const std::vector<double>& interpolationNodes, const std::vector<double> values) {
        splineConstraints_->push();
        addInterpolationNodes(interpolationNodes);
        Eigen::VectorXd solution = solve(values);
        splineConstraints_->pop();
        return solution;
    }

    Eigen::VectorXd BSplineStructure::solve(const std::vector<double> parameters) {
        splineConstraints_->update_b(parameters);
        int status = splineConstraints_->solve();
        QL_REQUIRE(status == 1, "Solution failed, returned " << status << std::endl);
        return splineConstraints_->getSolution();
    }

    Eigen::VectorXd BSplineStructure::getSolution() {
        return splineConstraints_->getSolution();
    }


    /*
    This function calculates value at x of all spline basis functions given the knot sequence, and also
    allows designating a degree exceeding the degree of the spline, this is needed for primitive
    value.
    */
    Eigen::VectorXd
    BSplineStructure::evaluateAll(double x, Size degree, Side side) const {
        Size p = (degree != -1) ? degree : degree_;
        BSplineEvaluator spline;

        // Excess degree
        Size e = (p > degree_) ? p - degree_ : 0;

        std::vector<double> knotsVector(knots_.begin(), knots_.end());

        // Pad the knot sequence at the ends if the degree is higher than degree_, this may happen when
        // integrating
        if (p > degree_) {
            knotsVector.insert(knotsVector.begin(), e, startPoint_);
            knotsVector.insert(knotsVector.end(), e, endPoint_);
        }

        // The "basis" size (some are then 0 functions)
        // TODO: we could also do the padding of 0 is then answer and not burden the evaluation
        Size n = knotsVector.size() - p - 1 + e * 2;

        // TODO the sidedness could be taken care of in the evaluation function, it just affects the
        // mu
        if (side == Side::Right) {
            spline = BSplineEvaluator(knotsVector, p);
        } else if (side == Side::Left) {
            std::reverse(knotsVector.begin(), knotsVector.end());
            for (auto& knot : knotsVector) {
                knot = -knot;
            }
            x = -x;
            spline = BSplineEvaluator(knotsVector, p);
        } else {
            throw std::invalid_argument("Side must be either 'left' or 'right'");
        }

        Eigen::VectorXd basisValues = spline.evaluateAll(x);

        if (side == Side::Left) {
            std::reverse(basisValues.begin(), basisValues.end());
        }

        // return std::vector<double>(basisValues.begin(), basisValues.end());
        return basisValues;
    }

    // TODO this should only be sensitive to the smoothness, not RT or not, nor degree
    void BSplineStructure::defaultKnotIndices() {
        if (interpolationSmoothness_ == InterpolationSmoothness::Default) {
            nKnots_ = nSimpleKnots_ + 2ULL * degree_;
            knotIndices_.assign(degree_, 0);
            for (Size i = 0; i < nSimpleKnots_; ++i) {
                knotIndices_.push_back(i);
            }
            knotIndices_.insert(knotIndices_.end(), degree_, nSimpleKnots_ - 1);
            //structure_.assign(1, degree_ + 1);
            //structure_.insert(structure_.end(), nSimpleKnots_ - 2, 1);
            //structure_.push_back(degree_ + 1);
        } else if (interpolationSmoothness_ ==
                   InterpolationSmoothness::Hermite) {
            nKnots_ = 2 * nSimpleKnots_ + 2ULL * (degree_ - 1);
            knotIndices_.assign(degree_ - 1, 0);
            for (Size i = 0; i < 2 * nSimpleKnots_; ++i) {
                knotIndices_.push_back(i / 2);
            }
            knotIndices_.insert(knotIndices_.end(), degree_ - 1, nSimpleKnots_ - 1);
            //structure_.assign(1, degree_ + 1);
            //structure_.insert(structure_.end(), nSimpleKnots_ - 2, 2);
            //structure_.push_back(degree_ + 1);
        } else {
            std::ostringstream oss;
            oss << "Interpolation smoothness: " << int(interpolationSmoothness_) << " not supported yet";
            throw std::runtime_error(oss.str());
        }
    }

    Eigen::VectorXd
    BSplineStructure::derivative_(double x, int nu, Size degree, double x0, Side side) const {
        return Eigen::VectorXd();
    }

    Eigen::SparseMatrix<double>
    BSplineStructure::singleDerivativeMatrix(Size degree, bool differenceOperator) const {
        const auto& knots = knots_;

        Size p = (degree != -1) ? degree : degree_;
        Size n = knots.size() - p - 1;

        if (p == 0) {
            return Eigen::SparseMatrix<double>(n + 1, n);
        }

        // Compute differences at p points apart
        Eigen::VectorXd differences(n + 1);
        for (Size i = 0; i <= n; ++i) {
            differences[i] = knots[p + i] - knots[i];
        }

        // Compute reciprocals for non-zero differences
        Eigen::VectorXd reciprocals = Eigen::VectorXd::Zero(n + 1);
        for (Size i = 0; i <= n; ++i) {
            if (differences[i] != 0) {
                reciprocals[i] = static_cast<double>(p) / differences[i];
            }
        }

        // Construct sparse pseudoinverse
        Eigen::SparseMatrix<double> t_matrix(n + 1, n + 1);
        for (Size i = 0; i <= n; ++i) {
            if (reciprocals[i] != 0) {
                t_matrix.insert(i, i) = reciprocals[i];
            }
        }

        // Construct the delta matrix
        Eigen::SparseMatrix<double> delta_matrix(n + 1, n);
        for (Size i = 0; i < n; ++i) {
            delta_matrix.insert(i, i) = 1.0;
            delta_matrix.insert(i + 1, i) = -1.0;
        }

        if (differenceOperator) {
            return delta_matrix;
        } else {
            return t_matrix * delta_matrix;
        }
    }

    Eigen::MatrixXd BSplineStructure::singleAntiDerivativeMatrix(Size degree,
                                                                double t0,
                                                                bool differenceOperator) const {

        if (std::isnan(t0)) {
            t0 = startPoint_;
        }

        Size p = (degree != -1) ? degree : degree_;
        Size e = (p > degree_) ? p - degree_ + 1 : 0;
        std::vector<double> knotsVector = knots_;

        if (p > degree_) {
            knotsVector.insert(knotsVector.begin(), e, startPoint_);
            knotsVector.insert(knotsVector.end(), e, endPoint_);
        }

        Size n = knotsVector.size() - p - 1 + 2 * e;

        // Linear constraint sets the anti-derivative to be 0 at t0
        Eigen::VectorXd shiftConstraint = evaluateAll(t0, p + 1).segment(0, n - 1);

        // Compute differences at p+1 points apart
        Eigen::VectorXd differences = Eigen::Map<Eigen::VectorXd>(knotsVector.data() + (p + 1), n) -
                                      Eigen::Map<Eigen::VectorXd>(knotsVector.data(), n);

        Eigen::VectorXd reciprocals = Eigen::VectorXd::Zero(differences.size());
        reciprocals = differences.cwiseInverse() * (p + 1);

        // Construct sparse matrices
        Eigen::SparseMatrix<double> tMatrix(n, n);
        for (Size i = 0; i < n; ++i) {
            if (reciprocals[i] != 0) {
                tMatrix.insert(i, i) = reciprocals[i];
            }
        }

        // Mask for non-zero splines
        Eigen::VectorXd mask = differences.unaryExpr([](double d) { return d != 0 ? 1.0 : 0.0; });

        Eigen::SparseMatrix<double> deltaMatrix(n, n - 1);
        for (Size i = 0; i < n - 1; ++i) {
            deltaMatrix.insert(i, i) = 1.0;
            deltaMatrix.insert(i + 1, i) = -1.0;
        }

        Eigen::MatrixXd resultMatrix;
        if (differenceOperator) {
            resultMatrix = (Eigen::MatrixXd(deltaMatrix) * mask.asDiagonal()).transpose();
        } else {
            resultMatrix = (Eigen::MatrixXd(tMatrix * deltaMatrix) * mask.asDiagonal()).transpose();
        }

        resultMatrix.row(0) = shiftConstraint.transpose();

        return resultMatrix;
    }

    //EigenMatrix BSplineStructure::derivativeMatrix(int nu,
    //                                              Size degree,
    //                                              bool differenceOperator,
    //                                              double t0) const {
    //    if (std::isnan(t0)) {
    //        t0 = startPoint_;
    //    }

    //    Size p = (degree != -1) ? degree : degree_;

    //    if (nu == 0) {
    //        Size n = (p <= degree_) ? knots_.size() - p - 1 : knots_.size() - p - 1 + (p - degree_) * 2;
    //        Eigen::SparseMatrix<double> identityMatrix(n, n);
    //        identityMatrix.setIdentity();
    //        return identityMatrix;
    //    } else if (nu >= 1) {
    //        Eigen::SparseMatrix<double> der =
    //            std::get<Eigen::SparseMatrix<double>>(derivativeMatrix(nu - 1, p - 1));
    //        Eigen::SparseMatrix<double> singleDer = singleDerivativeMatrix(p, differenceOperator);
    //        return (der * singleDer).eval();
    //    } else { // nu <= -1, i.e. antiderivative
    //        Eigen::MatrixXd singleAntiDer = singleAntiDerivativeMatrix(p, t0, differenceOperator);
    //        // if (std::isnan(t0)) {
    //        //     return derivativeMatrix(nu + 1, p + 1, differenceOperator) * singleAntiDer;
    //        // } else {
    //        if (nu == -1) {
    //            return (std::get<Eigen::SparseMatrix<double>>(
    //                        derivativeMatrix(nu + 1, p + 1, differenceOperator)) *
    //                    singleAntiDer)
    //                .eval();
    //        } else {
    //            return (std::get<Eigen::MatrixXd>(
    //                        derivativeMatrix(nu + 1, p + 1, differenceOperator, t0)) *
    //                    singleAntiDer)
    //                .eval();
    //        }
    //        //}
    //    }
    //}

    //Eigen::VectorXd BSplineStructure::derivative_(
    //    double x, int nu, Size degree, double x0, Splines::Side side) const {
    //    Size p = (degree != -1) ? degree : degree_;
    //    Splines::Side side_ = side;

    //    Eigen::VectorXd evaluate_all;

    //    if (std::find(knots_.begin(), knots_.end(), x) == knots_.end()) {
    //        evaluate_all = evaluateAll(x, p - nu, Splines::Side::Right);
    //    } else if (side_ == Splines::Side::Left || side_ == Splines::Side::Right) {
    //        evaluate_all = evaluateAll(x, p - nu, side_);
    //    } else if (side_ == Splines::Side::Average) {
    //        evaluate_all = (evaluateAll(x, p - nu, Splines::Side::Left) +
    //                        evaluateAll(x, p - nu, Splines::Side::Right)) /
    //                       2.0;
    //    } else { // Actual
    //        evaluate_all = evaluateAll(x, p - nu, Splines::Side::Right);
    //        Eigen::VectorXd evaluate_all_left = evaluateAll(x, p - nu, Splines::Side::Left);
    //        // Find indices where values disagree
    //        for (int i = 0; i < evaluate_all.size(); ++i) {
    //            if (evaluate_all[i] != evaluate_all_left[i]) {
    //                evaluate_all[i] = std::numeric_limits<double>::quiet_NaN();
    //            }
    //        }
    //    }

    //    if (nu >= 0) {
    //        return evaluate_all * std::get<Eigen::SparseMatrix<double>>(derivativeMatrix(nu, p));
    //    } else {
    //        return evaluate_all * std::get<Eigen::MatrixXd>(derivativeMatrix(nu, p));
    //    }
    //}

    static bool checkKnotIndices(const std::vector<Integer>& knotIndices, Integer degree, Size n) {

        if (knotIndices.front() != 0 || knotIndices.back() != n) {
            return false;
        }

        for (Size i = 0; i < knotIndices.size() - 1; ++i) {
            if (knotIndices[i] > knotIndices[i + 1] || knotIndices[i + 1] - knotIndices[i] > 1) {
                return false;
            }
        }

        std::unordered_map<Integer, Integer> knotCountMap;
        for (Integer index : knotIndices) {
            knotCountMap[index]++;
        }

        std::vector<Integer> counts;
        for (const auto& pair : knotCountMap) {
            counts.push_back(pair.second);
        }

        if (std::any_of(counts.begin(), counts.end(),
                        [degree](Integer value) { return value > degree + 1; })) {
            return false;
        }

        if (counts.front() != degree + 1 || counts.back() != degree + 1) {
            return false;
        }

        return true;
    }


}
