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
#include "bsplineevaluator.hpp"
#include "splineconstraints.hpp"
#include "splinesegment.hpp"
#include "ql/math/factorial.hpp"
#include "ql/termstructures/iterativebootstrap.hpp"
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
#include <iostream>

inline std::string vectorToString(const std::vector<Integer>& matrix) {
    std::ostringstream oss;

    oss << "{"; // Start each row with an opening brace

    bool first = true;
    for (auto value : matrix) {
        if (first == false) {
            oss << ", ";
        }
        oss << value;
        first = false;
    }

    oss << "}";

    return oss.str(); // Convert the string stream to a string
}


namespace QuantLib {
    /*
    Class for a spline structure holding the necessary information on knots, degree and constraints
    for a spline.
    */
    BSplineSegment::BSplineSegment(const std::vector<Real>& simpleKnots,
                                       Integer degree,
                                       const std::vector<Integer>& knotIndices,
                                       InterpolationSmoothness smoothness,
                                       InterpolationTransform interpolationTransform,
                                       Side side,
                                       const Size requiredPoints,
                                       const bool isGlobal)
    : simpleKnots_(simpleKnots), degree_(degree), knotIndices_(knotIndices),
      interpolationSmoothness_(smoothness),
    interpolationTransform_(interpolationTransform),
        side_(side), requiredPoints_(requiredPoints),
      //knotTolerance_(knotTolerance), rankTolerance_(rankTolerance),
      isGlobal_(isGlobal), startPoint_(simpleKnots.front()),
      endPoint_(simpleKnots.back()),
      nSimpleKnots_(static_cast<Integer>(simpleKnots.size())) {
        QL_REQUIRE(std::is_sorted(simpleKnots_.begin(), simpleKnots_.end(), std::less<>()),
                   "The x (simpleKnots) must be strictly increasing.");

//        std::cout << "Transform is " << static_cast<Integer>(interpolationTransform) << std::endl;
 
        if (knotIndices_.empty()) {
            defaultKnotIndices();
        } else {
            QL_REQUIRE(
                checkKnotIndices(knotIndices_, degree_, nSimpleKnots_ - 1),
                "The provided knots are invalid for the spline of the given degree.\n"
                "The sequence must be non-decreasing and cover all indexes, no index repeated "
                "more than degree+1 times, and the\n"
                "end knots must have exactly that count.\n"
                << vectorToString(knotIndices_) << "\nnSimpleKnots: " << nSimpleKnots_ 
                << "\t Degree: " << degree_ << "\n");
        }

        nKnots_ = static_cast<Integer>(knotIndices_.size());

        knots_.resize(nKnots_);
        for (Integer i = 0; i < nKnots_; ++i) {
            knots_[i] = simpleKnots_[knotIndices_[i]];
        }

        spline_ = BSplineEvaluator(knots_, degree_);
        //InterpolationTransform interpolationTransform = InterpolationTransform::Default;
        // Initialize transform and inverse based on the transform parameter
        switch (interpolationTransform) {
            case InterpolationTransform::Log:
                transform = [](Real t, Real x) { return std::log(x); };
                transformDerivative = [](Real t, Real x, Natural n) {
                    if (n == 0)
                        return std::log(x);
                    return ((n % 2 == 0) ? -1.0 : 1.0) * Factorial::get(n-1) /
                           std::pow(x, static_cast<Real>(n));
                };
                inverse = [](Real t, Real x) { return std::exp(x); };
                break;
            case InterpolationTransform::Exp:
                transform = [](Real t, Real x) { return std::exp(x); };
                transformDerivative = [](Real t, Real x, Integer n) {
                    return std::exp(x);
                };
                inverse = [](Real t, Real x) {
                    return std::log(x);
                };
                break;
            case InterpolationTransform::RateTime:
                transform = [](Real t, Real x) { return t * x; };
                transformDerivative = [](Real t, Real x, Integer n) {
                    return (n == 0) ? t * x : ((n == 1) ? t : 0.0);
                };
                inverse = [this](Real t, Real x) {
                    if (t > 0.0)
                        return x / t;
                    std::cout << "Using x0 " << this->x0_ << "\n";
                    return this->x0_;
                };
                break;
            case InterpolationTransform::Default:
                transform = [](Real t, Real x) { return x; };
                transformDerivative = [](Real t, Real x, Integer n) {
                    return (n == 0) ? x : ((n == 1) ? 1.0 : 0.0);
                };
                inverse = [](Real t, Real x) { return x; };
                break;
        }
    }

    // Copy constructor
    BSplineSegment::BSplineSegment(const BSplineSegment& other)
        : degree_(other.degree_),
          knotIndices_(other.knotIndices_),
          interpolationSmoothness_(other.interpolationSmoothness_),
          interpolationTransform_(other.interpolationTransform_),
          side_(other.side_),
          requiredPoints_(other.requiredPoints_),
          isGlobal_(other.isGlobal_),
          startPoint_(other.startPoint_), endPoint_(other.endPoint_), knots_(other.knots_),
          nSimpleKnots_(other.nSimpleKnots_), nKnots_(other.nKnots_) {
        spline_ = BSplineEvaluator(knots_, degree_); // Recreate spline evaluator
    }

    // Assignment operator
    BSplineSegment& BSplineSegment::operator=(const BSplineSegment& other) {
        // Check for self-assignment
        if (this != &other) {
            // Assign each member variable from the other object
            degree_ = other.degree_;
            knotIndices_ = other.knotIndices_;
            interpolationSmoothness_ = other.interpolationSmoothness_;
            interpolationTransform_ = other.interpolationTransform_;
            side_ = other.side_;
            requiredPoints_ = other.requiredPoints_;
            isGlobal_ = other.isGlobal_;
            startPoint_ = other.startPoint_;
            endPoint_ = other.endPoint_;
            knots_ = other.knots_;
            nSimpleKnots_ = other.nSimpleKnots_;
            nKnots_ = other.nKnots_;

            // Recreate the spline evaluator with the new state
            spline_ = BSplineEvaluator(knots_, degree_);
        }
        return *this; // Return *this to allow chained assignments
    }

    // Accessor for knot range
    std::pair<double, double> BSplineSegment::range() const {
        return {startPoint_, endPoint_};
    }

    // Accessor for knots_
    const std::vector<double>& BSplineSegment::knots() const {
        return knots_;
    }

    // Accessor for degree_
    Size BSplineSegment::degree() const {
        return degree_;
    }

    /*
    This function calculates value at x of all spline basis functions given the knot sequence, and also
    allows designating a degree exceeding the degree of the spline, this is needed for primitive
    value.
    */
    Eigen::VectorXd
    BSplineSegment::evaluateAll(Real x, Size degree, Side side) const {
        Size p = (degree != static_cast<Size>(-1)) ? degree : static_cast<Size>(degree_);
        BSplineEvaluator spline;

        // Excess degree
        Size e = (p > static_cast<Size>(degree_)) ? p - static_cast<Size>(degree_) : 0;

        std::vector<double> knotsVector(knots_.begin(), knots_.end());

        // Pad the knot sequence at the ends if the degree is higher than degree_, this may happen when
        // calculating anti-derivatives
        if (p > static_cast<Size>(degree_)) {
            knotsVector.insert(knotsVector.begin(), e, startPoint_);
            knotsVector.insert(knotsVector.end(), e, endPoint_);
        }

        // The "basis" size (some are then 0 functions)
        // TODO: we could also do the padding of 0 is then answer and not burden the evaluation
        // Size n = knotsVector.size() - p - 1 + e * 2;

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
    void BSplineSegment::defaultKnotIndices() {
        if (interpolationSmoothness_ == InterpolationSmoothness::Default) {
            nKnots_ = nSimpleKnots_ + 2 * degree_;
            knotIndices_.assign(degree_, 0);
            for (Integer i = 0; i < nSimpleKnots_; ++i) {
                knotIndices_.push_back(i);
            }
            knotIndices_.insert(knotIndices_.end(), degree_, nSimpleKnots_ - 1);
            //structure_.assign(1, degree_ + 1);
            //structure_.insert(structure_.end(), nSimpleKnots_ - 2, 1);
            //structure_.push_back(degree_ + 1);
        } else if (interpolationSmoothness_ ==
                   InterpolationSmoothness::Hermite) {
            nKnots_ = 2 * nSimpleKnots_ + 2 * (degree_ - 1);
            knotIndices_.assign(degree_ - 1, 0);
            for (Integer i = 0; i < 2 * nSimpleKnots_; ++i) {
                knotIndices_.push_back(i / 2);
            }
            knotIndices_.insert(knotIndices_.end(), degree_ - 1, nSimpleKnots_ - 1);
            //structure_.assign(1, degree_ + 1);
            //structure_.insert(structure_.end(), nSimpleKnots_ - 2, 2);
            //structure_.push_back(degree_ + 1);
        } else {
            std::ostringstream oss;
            oss << "Interpolation smoothness: " << static_cast<int>(interpolationSmoothness_) << " not supported yet";
            throw std::runtime_error(oss.str());
        }
    }

    Eigen::VectorXd
    BSplineSegment::derivative(Real x, Integer nu, Size degree, Real x0, Side side) const {
        return {};
    }

    // TODO implement sidedness
    Real BSplineSegment::value(const Eigen::VectorXd& coefficients, const Real t, Side side)  {
        if (t == 0.0 && interpolationTransform_ == InterpolationTransform::RateTime) {
            // TODO this works for Linear RT only, general solution requires derivative
            return value(coefficients, simpleKnots_[1]) / simpleKnots_[1];
        }
        //std::cout << t << ", " << spline_.value(coefficients, t) << ", "
        //          << inverse(t, spline_.value(coefficients, t)) << "\n"
        return inverse(t, spline_.value(coefficients, t));
    }

    Eigen::SparseMatrix<Real>
    BSplineSegment::singleDerivativeMatrix(Size degree, bool differenceOperator) const {
        const auto& knots = knots_;

        Size p = (degree != static_cast<Size>(-1)) ? degree : degree_;
        Size n = knots.size() - p - 1;

        if (p == 0) {
            return {static_cast<Eigen::Index>(n+1), static_cast<Eigen::Index>(n)};
        }

        // Compute differences at p points apart
        Eigen::VectorXd differences(n + 1);
        for (Size i = 0; i <= n; ++i) {
            differences[i] = knots[p + i] - knots[i];
        }

        // Compute reciprocals for non-zero differences
        Eigen::VectorXd reciprocals = Eigen::VectorXd::Zero(n + 1);
        for (Size i = 0; i <= n; ++i) {
            if (differences[i] != 0.0) {
                reciprocals[i] = static_cast<Real>(p) / differences[i];
            }
        }

        // Construct sparse pseudo-inverse
        Eigen::SparseMatrix<double> t_matrix(n + 1, n + 1);
        for (Size i = 0; i <= n; ++i) {
            if (reciprocals[i] != 0.0) {
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

    Eigen::MatrixXd BSplineSegment::singleAntiDerivativeMatrix(Size degree,
                                                                double t0,
                                                                bool differenceOperator) const {

        if (std::isnan(t0)) {
            t0 = startPoint_;
        }

        Size p = (degree != static_cast<Size>(-1)) ? degree : static_cast<Size>(degree_);
        Size e = (p > static_cast<Size>(degree_)) ? p - static_cast<Size>(degree_) + 1 : 0;
        std::vector<double> knotsVector = knots_;

        if (p > static_cast<Size>(degree_)) {
            knotsVector.insert(knotsVector.begin(), e, startPoint_);
            knotsVector.insert(knotsVector.end(), e, endPoint_);
        }

        Size n = knotsVector.size() - p - 1 + 2 * e;

        // Linear constraint sets the antiderivative to be 0 at t0
        Eigen::VectorXd shiftConstraint = evaluateAll(t0, p + 1).segment(0, n - 1);

        // Compute differences at p+1 points apart
        Eigen::VectorXd differences = Eigen::Map<Eigen::VectorXd>(knotsVector.data() + (p + 1), n) -
                                      Eigen::Map<Eigen::VectorXd>(knotsVector.data(), n);

        Eigen::VectorXd reciprocals = Eigen::VectorXd::Zero(differences.size());
        reciprocals = differences.cwiseInverse() * (p + 1);

        // Construct sparse matrices
        Eigen::SparseMatrix<double> tMatrix(n, n);
        for (Size i = 0; i < n; ++i) {
            if (reciprocals[i] != 0.0) {
                tMatrix.insert(i, i) = reciprocals[i];
            }
        }

        // Mask for non-zero splines
        Eigen::VectorXd mask = differences.unaryExpr([](Real d) { return d != 0.0 ? 1.0 : 0.0; });

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

    //EigenMatrix BSplineSegment::derivativeMatrix(int nu,
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

    //Eigen::VectorXd BSplineSegment::derivative(
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
}

bool checkKnotIndices(const std::vector<Integer>& knotIndices, Integer degree, Integer n) {

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
    counts.reserve(knotCountMap.size());
    for (const auto& pair : knotCountMap) {
        counts.emplace_back(pair.second);
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
