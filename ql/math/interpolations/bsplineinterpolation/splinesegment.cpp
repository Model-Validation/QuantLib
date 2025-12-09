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
#include "bsplineevaluator.hpp"
#include "splineconstraints.hpp"
#include "splinesegment.hpp"
#include "rw.h"
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

// Force recompilation - timestamp issue
namespace QuantLib {
    /*
    Class for a spline structure holding the necessary information on knots, degree and constraints
    for a spline.
    */
    BSplineSegment::BSplineSegment(const std::vector<Real>& simpleKnots,
                                       Integer degree,
                                       const std::vector<Integer>& knotIndices,
                                       InterpolationSmoothnessEnum smoothness,
                                       InterpolationTransformEnum interpolationTransform,
                                       SideEnum side,
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


        if (this->knotIndices_.empty()) {
            defaultKnotIndices();
        } else {
            QL_REQUIRE(
                checkKnotIndices(this->knotIndices_, this->degree_, this->nSimpleKnots_ - 1),
                "The provided knots are invalid for the spline of the given degree.\n"
                "The sequence must be non-decreasing and cover all indexes, no index repeated "
                "more than degree+1 times, and the\n"
                "end knots must have exactly that count.\n"
                << vectorToString(knotIndices_) << "\nnSimpleKnots: " << nSimpleKnots_
                << "\t Degree: " << degree_ << "\n");
        }

        this->nKnots_ = static_cast<Integer>(this->knotIndices_.size());

        knots_.resize(nKnots_);
        for (Integer i = 0; i < nKnots_; ++i) {
            this->knots_[i] = this->simpleKnots_[this->knotIndices_[i]];
        }

        this->spline_ = BSplineEvaluator(this->knots_, this->degree_);

        // InterpolationTransformEnum interpolationTransform = TransformDefault;
        // Initialize transform and inverse based on the transform parameter
        switch (interpolationTransform) {
            case TransformLog:
                this->transform = [](Real t, Real x) { return std::log(x); };
                this->transformDerivative = [](Real t, Real x, Natural n) {
                    if (n == 0)
                        return std::log(x);
                    return ((n % 2 == 0) ? -1.0 : 1.0) * Factorial::get(n-1) /
                           std::pow(x, static_cast<Real>(n));
                };
                this->inverse = [](Real t, Real x) { return std::exp(x); };
                break;
            case TransformExp:
                this->transform = [](Real t, Real x) { return std::exp(x); };
                this->transformDerivative = [](Real t, Real x, Integer n) {
                    return std::exp(x);
                };
                this->inverse = [](Real t, Real x) {
                    return std::log(x);
                };
                break;
            case TransformRateTime:
                this->transform = [](Real t, Real x) { return t * x; };
                this->transformDerivative = [](Real t, Real x, Integer n) {
                    return (n == 0) ? t * x : ((n == 1) ? t : 0.0);
                };
                this->inverse = [this](Real t, Real x) {
                    if (t > 0.0)
                        return x / t;
                    // Debug output removed to prevent console popups
                    // std::cout << "Using rateTimeX0 " << this->rateTimeX0 << "\n";
                    return this->rateTimeX0;
                };
                break;
            case TransformRateTimeAnnualToContinuous:
                this->transform = [](Real t, Real x) { return t * std::log1p(x); };
                this->transformDerivative = [](Real t, Real x, Integer n) {
                    if (n == 0)
                        return t * std::log1p(x);
                    return t * ((n % 2 == 0) ? -1.0 : 1.0) * Factorial::get(n - 1) /
                           std::exp(std::log1p(x) * static_cast<Real>(n));
                           // std::pow(x + 1, static_cast<Real>(n));
                };
                this->inverse = [this](Real t, Real x) {
                    if (t > 0.0)
                        return std::expm1(x / t);
                    // Debug output removed to prevent console popups
                    // std::cout << "Using rateTimeX0 " << this->rateTimeX0 << "\n";
                    return this->rateTimeX0; // TODO fix this, need to take a derivative or generically return nan to be picked up by caller
                };
                break;
            case TransformContinuousToAnnual:
                this->transform = [](Real t, Real x) { return std::expm1(x); };
                this->transformDerivative = [](Real t, Real x, Integer n) {
                    if (n == 0)
                        return std::expm1(x);
                    return std::exp(x);
                };
                this->inverse = [this](Real t, Real x) {
                    return std::log1p(x);
                };
                break;
            case TransformContinuousToSimple:
                this->transform = [](Real t, Real x) {
                    if (t == 0.0) {
                        return x;
                    }
                    return std::expm1(t * x) / t;
                };
                this->transformDerivative = [](Real t, Real x, Integer n) {
                    if (n == 0) {
                        if (t == 0.0) {
                            return x;
                        }
                        return std::expm1(t * x) / t;
                    }
                    return std::exp(t * x) * std::pow(t, n - 1);
                };

                this->inverse = [this](Real t, Real x) {
                    if (t == 0.0) {
                        return x;
                    }
                    return std::log1p(t * x) / t;
                };
                break;
            case TransformDefault:
                this->transform = [](Real t, Real x) { return x; };
                this->transformDerivative = [](Real t, Real x, Integer n) {
                    return (n == 0) ? x : ((n == 1) ? 1.0 : 0.0);
                };
                this->inverse = [](Real t, Real x) { return x; };
                break;
        }
    }

    // Copy constructor
    BSplineSegment::BSplineSegment(const BSplineSegment& other)
        : simpleKnots_(other.simpleKnots_),
          degree_(other.degree_),
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
    std::pair<Real, Real> BSplineSegment::range() const {
        return {startPoint_, endPoint_};
    }

    // Accessor for knots_
    const std::vector<Real>& BSplineSegment::knots() const {
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
    BSplineSegment::evaluateAll(Real t, Size degree, SideEnum side) const {
        const Size p = (degree != static_cast<Size>(-1)) ? degree : static_cast<Size>(degree_);
        BSplineEvaluator spline;

        // Excess degree
        const Size e = (p > static_cast<Size>(degree_)) ? p - static_cast<Size>(degree_) : 0;

        std::vector<Real> knotsVector(knots_.begin(), knots_.end());
        

        // Adjust knot sequence for different degree
        if (p > static_cast<Size>(degree_)) {
            // Anti-derivative: pad the knot sequence at the ends
            knotsVector.insert(knotsVector.begin(), e, startPoint_);
            knotsVector.insert(knotsVector.end(), e, endPoint_);
        } else if (p < static_cast<Size>(degree_)) {
            // Derivative: need to handle knot vector adjustment properly
            Size deficit = static_cast<Size>(degree_) - p;
            
            // Use the standard reduction strategy for ALL cases (including degree 0)
            QL_REQUIRE(knotsVector.size() > 2 * deficit,
                       "Insufficient knots to reduce degree from " << degree_ 
                       << " to " << p << ". Knot vector has size " << knotsVector.size()
                       << " but need to remove " << 2 * deficit << " knots.");
            
            // Remove 'deficit' knots from each end - this works for all degrees including 0
            knotsVector.erase(knotsVector.begin(), knotsVector.begin() + deficit);
            knotsVector.erase(knotsVector.end() - deficit, knotsVector.end());
        }

        // The "basis" size (some are then 0 functions)
        // TODO: we could also do the padding of 0 is then answer and not burden the evaluation
        // Size n = knotsVector.size() - p - 1 + e * 2;

        // Implement the reverse/negate transformation for Left side
        QL_ASSERT(side == SideRight || side == SideLeft,
                  "Side must be either 'Left' or 'Right'");
        
        Real x = t;
        if (side == SideLeft) {
            // TEMPORARY DEBUG: Print what we're transforming
            std::ostringstream debugMsg;
            debugMsg << "DEBUG Left transformation for t=" << t << ":\n";
            debugMsg << "  Original knots: ";
            for (auto k : knotsVector) debugMsg << k << " ";
            debugMsg << "\n";
            
            // Transform knot vector for left-sided evaluation
            std::reverse(knotsVector.begin(), knotsVector.end());
            
            debugMsg << "  After reverse: ";
            for (auto k : knotsVector) debugMsg << k << " ";
            debugMsg << "\n";
            
            for (auto& knot : knotsVector) {
                knot = -knot;
            }
            x = -t;
            
            debugMsg << "  After negate: ";
            for (auto k : knotsVector) debugMsg << k << " ";
            debugMsg << "\n  x = " << x;
            
            // Output debug info
            std::cout << debugMsg.str() << std::endl;
        }
        
        // Create evaluator with (possibly transformed) knot vector
        spline = BSplineEvaluator(knotsVector, p);
        
        // DEBUG: Verify the evaluator was created with correct dimensions
        Size expectedNumBasis = knotsVector.size() - p - 1;
        // The evaluator should have this many basis functions
        
        Eigen::VectorXd basisValues = spline.evaluateAll(x);
        
        // DEBUG: Report actual sizes
        if (side == SideLeft) {
            std::ostringstream oss;
            oss << "Left-side evaluation debug:\n"
                << "  knotsVector.size() = " << knotsVector.size() << "\n"
                << "  degree (p) = " << p << "\n"  
                << "  expectedNumBasis = " << expectedNumBasis << "\n"
                << "  basisValues.size() = " << basisValues.size() << "\n"
                << "  First element = " << (basisValues.size() > 0 ? basisValues[0] : -999) << "\n"
                << "  Last element = " << (basisValues.size() > 0 ? basisValues[basisValues.size()-1] : -999);
            throw std::runtime_error(oss.str());
        }
        
        if (side == SideLeft) {
            // Reverse the basis values to account for the transformation
            // Wait - what if reverseInPlace() is somehow truncating?
            // Let's try manual reverse
            Eigen::VectorXd reversed(basisValues.size());
            for (Eigen::Index i = 0; i < basisValues.size(); ++i) {
                reversed[i] = basisValues[basisValues.size() - 1 - i];
            }
            basisValues = reversed;
        }

        // return std::vector<double>(basisValues.begin(), basisValues.end());
        return basisValues;
    }

    // TODO: this should only be sensitive to the smoothness, not RT or not, nor degree
    void BSplineSegment::defaultKnotIndices() {
        if (interpolationSmoothness_ == SmoothnessDefault) {
            nKnots_ = nSimpleKnots_ + 2 * static_cast<Integer>(degree_);
            knotIndices_.assign(degree_, 0);
            for (Integer i = 0; i < nSimpleKnots_; ++i) {
                knotIndices_.push_back(i);
            }
            knotIndices_.insert(knotIndices_.end(), degree_, nSimpleKnots_ - 1);
            //structure_.assign(1, degree_ + 1);
            //structure_.insert(structure_.end(), nSimpleKnots_ - 2, 1);
            //structure_.push_back(degree_ + 1);
        } else if (interpolationSmoothness_ ==
                   SmoothnessHermite) {
            nKnots_ = 2 * nSimpleKnots_ + 2 * static_cast<Integer>((degree_ - 1));
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

    Eigen::VectorXd BSplineSegment::derivativeFunctional(Real x,
                                               Integer nu,
                                               Size degree,
                                               Real x0,
                                               SideEnum side)
        const // Assuming Side is an enum with values Left, Right, Average, Actual, etc.
    {
        QL_REQUIRE(this->interpolationTransform_ == TransformDefault,
                   "Derivative for transformed curves not implemented yet.");
        // Determine the degree of the spline
        const Size p = (degree != static_cast<Size>(-1)) ? degree : this->degree_;
        const SideEnum actualSide = (side != SideDefault) ? side : this->side_;

        Eigen::VectorXd evaluateAll;

        // Handle the logic for evaluating the spline derivative
        if (std::find(this->knots_.begin(), this->knots_.end(), x) == this->knots_.end()) {
            // x is not in the knots
            evaluateAll = this->evaluateAll(x, p - nu, SideRight);
        } else if (actualSide == SideLeft || actualSide == SideRight) {
            // Evaluate on the specified side (left or right)
            evaluateAll = this->evaluateAll(x, p - nu, actualSide);
        } else if (actualSide == SideAverage) {
            // Average the evaluations on the left and right
            const Eigen::VectorXd evaluateLeft = this->evaluateAll(x, p - nu, SideLeft);
            const Eigen::VectorXd evaluateRight = this->evaluateAll(x, p - nu, SideRight);
            evaluateAll = (evaluateLeft + evaluateRight) / 2.0;
        } else { // SideActual
            // Evaluate on the right and compare with the left side
            Eigen::VectorXd evaluateRight = this->evaluateAll(x, p - nu, SideRight);
            Eigen::VectorXd evaluateLeft = this->evaluateAll(x, p - nu, SideLeft);

            // Find indices where the values differ
            Eigen::VectorXd evaluateAllResult = evaluateRight;
            for (int i = 0; i < evaluateAllResult.size(); ++i) {
                if (evaluateRight[i] != evaluateLeft[i]) {  // NOLINT(clang-diagnostic-float-equal)
                    evaluateAllResult[i] = std::numeric_limits<Real>::quiet_NaN();
                }
            }
            evaluateAll = evaluateAllResult;
        }

        // Multiply by the derivative matrix
        if (nu >= 0) {
            const Eigen::SparseMatrix<Real> derivative = this->derivativeMatrix(nu, p);
            Eigen::VectorXd product = derivative.transpose() * evaluateAll;
            return product;
        }
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        const Eigen::SparseMatrix<Real> antiDerivative = this->antiDerivativeMatrix(nu, p, x0).transpose();

        solver.analyzePattern(antiDerivative);
        solver.factorize(antiDerivative);


        Eigen::VectorXd product = solver.solve(evaluateAll);
        return product;
    }


    // TODO: implement sidedness
    Real BSplineSegment::value(const Eigen::VectorXd& coefficients,
                               Real t,
                               Integer nu,
                               SideEnum side) {
        QL_ASSERT(nu == 0 || this->interpolationTransform_ == TransformDefault,
                  "Derivative value not implemented for transformed curves");
        if (t == 0.0 && this->interpolationTransform_ == TransformRateTime) {
            // TODO: this works for Linear RT only, general solution requires derivative
            return value(coefficients, this->simpleKnots_[1], nu, SideLeft) /
                   this->simpleKnots_[1];
        }

        if (nu == 0) {
            return inverse(t, this->spline_.value(coefficients, t));
            //return coefficients.dot(derivativeFunctional(t, nu, -1, this->startPoint_, side));
        }
        // Multiply by the derivative matrix
        if (nu > 0) {
            if (nu > static_cast<Integer>(this->degree_)) {
                return 0.0;
            }
            if (this->derivativeCache_.find(nu) == this->derivativeCache_.end()) {
                const auto derivative = ext::make_shared<Eigen::SparseMatrix<Real>>(this->derivativeMatrix(nu, this->degree_));
                const std::vector<Real> knotVector(this->knots_.begin() + nu, this->knots_.end() - nu);
                const auto spline = ext::make_shared<BSplineEvaluator>(knotVector, this->degree_ - nu);
                this->splineCache_.emplace(nu, spline);
                this->derivativeCache_.emplace(nu, derivative);
                return spline->value(derivative->operator*(coefficients), t);
            }
            return this->splineCache_.at(nu)->value(
                this->derivativeCache_.at(nu)->operator*(coefficients), t);
        } else {
            if (this->antiDerivativeCache_.find(nu) == this->antiDerivativeCache_.end()) {
                const auto antiDerivative =
                    ext::make_shared<Eigen::SparseMatrix<Real>>(this->antiDerivativeMatrix(nu, this->degree_, this->startPoint_));
                //std::vector<Real> knotVector(this->knots_.size() - static_cast<Size>(2) * nu);
                std::vector<Real> knotVector;
                knotVector.reserve(this->knots_.size() - static_cast<Size>(2) * nu); // Pre-allocate space

                // Add the front value
                knotVector.insert(knotVector.end(), -nu, this->knots_.front());
                // Add the original vector
                knotVector.insert(knotVector.end(), this->knots_.begin(), this->knots_.end());
                // Add the back value
                knotVector.insert(knotVector.end(), -nu, this->knots_.back());

                const auto spline(
                    ext::make_shared<BSplineEvaluator>(knotVector, this->degree_ - nu));
                this->splineCache_.emplace(nu, spline);
                auto solver = ext::make_shared<Eigen::SparseLU<Eigen::SparseMatrix<Real>>>();
                solver->analyzePattern(*antiDerivative);
                solver->factorize(*antiDerivative);
                this->antiDerivativeCache_.emplace(nu, solver);
                return spline->value(solver->solve(coefficients), t);
            }
            return this->splineCache_.at(nu)->value(
                this->antiDerivativeCache_.at(nu)->solve(coefficients), t);
        }
    }

    Eigen::VectorXd BSplineSegment::valueFunctional(Real t, SideEnum side) const {
        QL_ASSERT(this->interpolationTransform_ == TransformDefault,
                  "Value operator not implemented for transformed curves");
        QL_ASSERT(side == SideRight || side == SideLeft,
                  "valueFunctional only accepts sides 'Right' and 'Left', consider using "
                  "derivativeFunctional with nu=0");

        return evaluateAll(t, static_cast<Size>(-1), side);
    }

    Eigen::SparseMatrix<Real>
    BSplineSegment::singleDerivativeMatrix(Size degree, bool differenceOperator) const {
        const Integer p = static_cast<Integer>((degree != static_cast<Size>(-1)) ? degree : this->degree_);
        const Integer e = static_cast<Integer>(p) - static_cast<Integer>(this->degree_);
        // Imply the size of proper knot vector size by removing/adding 2*e knots at ends of the sequence corresponding to degree_
        const Integer n = static_cast<Integer>(this->knots_.size()) -
                          2 * static_cast<Integer>(this->degree_) + p - 1;
        const Integer m = static_cast<Integer>(this->knots_.size());

        if (p == 0) {
            return {n-1, n};
        }

        // Construct the delta matrix
        Eigen::SparseMatrix<Real> deltaMatrix(n - 1, n);
        for (Integer i = 0; i < n - 1; ++i) {
            deltaMatrix.insert(i, i) = -1.0;
            deltaMatrix.insert(i, i + 1) = 1.0;
        }

        if (differenceOperator) {
            return deltaMatrix;
        }

        // Compute differences at p points apart
        //Is there a quick way to create a diagonal sparse matrix in Eigen from an Eigen Vector?
        Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(n - 1);
        for (Integer i = 1; i < n; ++i) {
            if (const Real diff =
                    this->knots_[std::min(p + i - e, m - 1)] - this->knots_[std::max(i - e, 0)];
                diff != 0.0) {
                diagonal[i - 1] = static_cast<Real>(p) / diff;
            }
        }

        // This returns an (n-1)xn sparse matrix
        Eigen::SparseMatrix<Real> result =
            Eigen::DiagonalMatrix<Real, Eigen::Dynamic>(diagonal) * deltaMatrix;
        result.makeCompressed();
        return result;
    }

    Eigen::SparseMatrix<Real> BSplineSegment::singleAntiDerivativeMatrix(
        Size degree, double t0, bool differenceOperator) const {

        // Handle t0 as the start point if it is NaN
        if (std::isnan(t0)) {
            t0 = this->startPoint_;
        }

        // Determine the degree and any required endpoint padding
        const Integer p = static_cast<Integer>((degree != static_cast<Size>(-1)) ? degree : this->degree_);
        const Integer e = static_cast<Integer>(p+1) - static_cast<Integer>(this->degree_);
        //const Size p = (degree != static_cast<Size>(-1)) ? degree : this->degree_;
        //const Size e = (p > this->degree_) ? p - this->degree_: 0;
        //std::vector<Real> knotsVector = this->knots_;

        //// Pad the knot vector for larger degree
        // if (e > 0) {
        //     knotsVector.insert(knotsVector.begin(), e, this->startPoint_);
        //     knotsVector.insert(knotsVector.end(), e, this->endPoint_);
        // }
        const Integer n = static_cast<Integer>(this->knots_.size()) - p - 2 + 2 * e;

        Eigen::SparseMatrix<Real> augmentedDerivativeMatrix =
            singleDerivativeMatrix(p + 1, differenceOperator);
        // This will be (n-)xn matrix
        QL_ASSERT(augmentedDerivativeMatrix.rows() == n - 1 &&
                      augmentedDerivativeMatrix.cols() == n,
                  "Derivative matrix has dimensions ("
                      << augmentedDerivativeMatrix.rows() << "," << augmentedDerivativeMatrix.cols()
                      << "). Expected (" << n - 1 << "," << n << ").");

        augmentedDerivativeMatrix.reserve(augmentedDerivativeMatrix.nonZeros() + p + 2);
        augmentedDerivativeMatrix.conservativeResize(n, n);
        Eigen::VectorXd shiftConstraint = evaluateAll(t0, p + 1);

        for (Integer j = 0; j < n; ++j) {
            if (shiftConstraint[j] != 0.0) {
                augmentedDerivativeMatrix.insert(n - 1, j) = shiftConstraint[j];
            }
        }

        augmentedDerivativeMatrix.makeCompressed();

        return augmentedDerivativeMatrix;
        // Linear
        //// Compute differences at p+1 points apart
        //Eigen::VectorXd differences(n);
        //for (Size i = 0; i < n; ++i) {
        //    differences[i] = knotsVector[p + 1 + i] - knotsVector[i];
        //}

        //// Construct deltaMatrix (a difference operator matrix)
        //Eigen::SparseMatrix<double> deltaMatrix(n, n - 1);
        //for (Size i = 0; i < n - 1; ++i) {
        //    deltaMatrix.insert(i, i) = 1.0;
        //    deltaMatrix.insert(i + 1, i) = -1.0;
        //}

        //// Masking matrix
        //Eigen::MatrixXd zMatrix = Eigen::MatrixXd::Zero(n, n);
        //for (Size i = 0; i < n; ++i) {
        //    if (differences[i] != 0.0) {
        //        zMatrix(i, i) = 1.0;
        //    }
        //}
        //// Handle trimming based on `f`
        //const Size f = (e > 0) ? 1 : 0;
        //zMatrix = zMatrix.block(0, f, n, n - static_cast<Size>(2) * f);

        //// Construct the final result matrix
        //// Use pseudo-inverse of (shiftConstraint stacked with deltaMatrix)
        //Eigen::MatrixXd combinedMatrix(n, n);
        //combinedMatrix.row(0) = shiftConstraint;
        //if (differenceOperator) {
        //    combinedMatrix.block(1, 0, n - 1, n - 1) = deltaMatrix.toDense();
        //} else {
        //    Eigen::SparseMatrix<double> tMatrix(n, n);
        //    for (Size i = 0; i < n; ++i) {
        //        if (differences[i] != 0.0) {
        //            tMatrix.insert(i, i) = static_cast<Real>(p + 1) / differences[i];
        //        }
        //    }
        //    // Use pseudo-inverse of (shiftConstraint stacked with tMatrix * deltaMatrix)
        //    combinedMatrix.block(1, 0, n - 1, n - 1) = (tMatrix * deltaMatrix).toDense();
        //}
        //Eigen::MatrixXd pseudoInverse =
        //    combinedMatrix.completeOrthogonalDecomposition().pseudoInverse();
        //return pseudoInverse.block(0, 1, n - 1, n) * zMatrix;
    }


    Eigen::SparseMatrix<Real> BSplineSegment::derivativeMatrix(
        const Integer nu, const Size degree, const bool differenceOperator) const {
        QL_ASSERT(nu >= 0, "Only non-negative nu supported, got nu=" << nu);

        const Size p = (degree != static_cast<Size>(-1)) ? degree : this->degree_;

        // Handle the nu == 0 case
        if (nu == 0) {
            const Size n = this->knots_.size() - 2 * this->degree_ + p - 1;
            //}
            Eigen::SparseMatrix<Real> identityMatrix(n, n);
            identityMatrix.setIdentity();
            return identityMatrix;
        }

        const Eigen::SparseMatrix<Real> singleDerivative =
            singleDerivativeMatrix(p, differenceOperator);

        if (p == 0) {
            return singleDerivative;
        }

        // Handle positive derivative and positive degree
        const Eigen::SparseMatrix<Real> recursionDerivative = derivativeMatrix(nu - 1, p - 1);
        Eigen::SparseMatrix<Real> result = recursionDerivative * singleDerivative;
        result.makeCompressed();
        return result;
    }

    Eigen::SparseMatrix<Real>
    BSplineSegment::antiDerivativeMatrix(Integer nu,
                                         Size degree,
                                         const Eigen::VectorXd& t0,
                                         bool differenceOperator) const {
        QL_ASSERT(nu <= 0, "Only non-positive nu supported, got nu=" << nu);

        const Size p = (degree != static_cast<Size>(-1)) ? degree : this->degree_;
        const Size n = this->knots_.size() - 2 * this->degree_ + p - 1;

        // Handle the nu == 0 case
        if (nu == 0) {
            Eigen::SparseMatrix<Real> identityMatrix(n, n);
            identityMatrix.setIdentity();
            return identityMatrix;
        }

        // Handle anti-derivative (nu < 0)
        // Handle default t0 value
        Real t0Value;
        Eigen::VectorXd t1Value;

        // If t0 is an Eigen::VectorXd, extract the first element and prepare the tail
        if (t0.size() >= 1) {
            t0Value = t0[0]; // Use the first element for this call
            t1Value = t0.tail(t0.size() - 1); // t1Value is the tail (t0[1:]), can be length 0
        } else {
            t0Value = std::numeric_limits<Real>::quiet_NaN();
            t1Value = Eigen::VectorXd(1);
            t1Value << std::numeric_limits<Real>::quiet_NaN();
        }

        // Recursive call with t1

        const Eigen::SparseMatrix<Real> antiDerivative =
            antiDerivativeMatrix(nu + 1, p + 1, t1Value, differenceOperator);

        const Size m = antiDerivative.rows();

        // Call singleAntiDerivativeMatrix using t0Value (Real) or quiet_NaN() if t1Vector is empty
        Eigen::SparseMatrix<Real> singleAntiDerivative =
            singleAntiDerivativeMatrix(p, t0Value, differenceOperator);
        singleAntiDerivative.conservativeResize(m, m);
        for (Size i = n + 1; i < m; ++i) {
            singleAntiDerivative.insert(i, i) = 1.0;
        }
        Eigen::SparseMatrix<Real> result = antiDerivative * singleAntiDerivative;
        result.makeCompressed();
        return result;
    }

    std::string vectorToString(const std::vector<QuantLib::Integer>& matrix) {
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


    bool checkKnotIndices(const std::vector<Integer>& knotIndices, Size degree, Integer n) {

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
                        [degree](const Size value) { return value > degree + 1; })) {
            return false;
        }

        if (static_cast<Size>(counts.front()) != degree + 1 || static_cast<Size>(counts.back()) != degree + 1) {
            return false;
        }

        return true;
    }
}
