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

// ReSharper disable CppClangTidyBugproneNarrowingConversions
#include "splineconstraints.hpp"
#include "eigenutilities.hpp"
#include <ql/errors.hpp>
#include <ql/types.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <scs_types.h>
#include <tuple>
#include <utility>
#include <vector>

namespace QuantLib {
    SplineConstraints::SplineConstraints(Size numVariables,
                                         const std::vector<std::vector<double>>& P_quadForm,
                                         const std::vector<std::vector<double>>& A_constraints,
                                         const std::vector<double>& b_rhs,
                                         const std::vector<double>& c_linearForm,
                                         Size numEqualities,
                                         Size numInequalities,
                                         bool fitData,
                                         double epsAbsolute,
                                         double epsRelative,
                                         double epsInfeasible)
    : numVariables_(numVariables), fitData_(fitData), epsAbsolute_(epsAbsolute),
      epsRelative_(epsRelative), epsInfeasible_(epsInfeasible) {
        
        // SCS-ordered constructor: constraints are already in the correct order
        // Equalities first, then inequalities
        
        numConstraints_ = numEqualities + numInequalities;
        QL_REQUIRE(A_constraints.size() == numConstraints_,
                   "A_constraints size must equal numEqualities + numInequalities. "
                   "Provided: " << A_constraints.size() << ", expected: " << numConstraints_);
        
        numEqualities_ = numEqualities;
        numInequalities_ = numInequalities;
        
        // Set up P matrix
        if (!P_quadForm.empty()) {
            P_ = EigenUtilities::convertToEigenSparseMatrix(P_quadForm);
        } else {
            P_ = Eigen::SparseMatrix<double>(numVariables, numVariables); // zero matrix
        }
        
        QL_REQUIRE(static_cast<Size>(P_.rows()) == numVariables &&
                       static_cast<Size>(P_.cols()) == numVariables,
                   "Matrix P_quadForm must be " << numVariables << "x" << numVariables
                                                << ". Provided: " << P_.rows() << "x" << P_.cols());
        
        // Set up A matrix (already in correct order)
        A_triplets_ = EigenUtilities::convertToTriplets(A_constraints);
        
        // Set up b vector
        QL_REQUIRE(b_rhs.size() == numConstraints_,
                   "b_rhs size must equal numConstraints");
        b_list_ = std::vector<double>(b_rhs);
        b_ = Eigen::Map<Eigen::VectorXd>(b_list_.data(), b_list_.size());
        
        // Set up c vector
        c_list_ = c_linearForm.empty() ? std::vector<Real>(numVariables, 0.0) :
                                         std::vector<Real>(c_linearForm);
        c_ = Eigen::Map<Eigen::VectorXd>(c_list_.data(), c_list_.size());
        
        // Build constraint types based on counts
        constraintTypes_.clear();
        constraintTypes_.reserve(numConstraints_);
        for (Size i = 0; i < numEqualities_; ++i) {
            constraintTypes_.push_back(ConstraintType::Equal);
        }
        for (Size i = 0; i < numInequalities_; ++i) {
            constraintTypes_.push_back(ConstraintType::LessEqual);
        }
        
        // Mark as already ordered (SCS ordering is the contract)
        isOrdered_ = true;
        
        // No permutation needed - we maintain strict SCS ordering
        
        numParameters_ = 0;
        hasParameters_ = false;
    }
    
    SplineConstraints::SplineConstraints(Size numVariables,
                                         const std::vector<std::vector<double>>& P_quadForm,
                                         const std::vector<std::vector<double>>& A_constraints,
                                         const std::vector<double>& b_rhs,
                                         const std::vector<double>& c_linearForm,
                                         const std::vector<ConstraintType>& constraintTypes,
                                         bool fitData,
                                         double epsAbsolute,
                                         double epsRelative,
                                         double epsInfeasible)
    : fitData_(fitData), numVariables_(numVariables), numEqualities_(0), numInequalities_(0),
      numParameters_(0), scsDataIsUpToDate_(false), isSolved_(false), epsAbsolute_(epsAbsolute),
      epsRelative_(epsRelative), epsInfeasible_(epsInfeasible) {
        if (!P_quadForm.empty()) {
            P_ = EigenUtilities::convertToEigenSparseMatrix(P_quadForm);
        } else {
            P_ = Eigen::SparseMatrix<double>(numVariables, numVariables); // zero matrix
        }

        QL_REQUIRE(static_cast<Size>(P_.rows()) == numVariables &&
                       static_cast<Size>(P_.cols()) == numVariables,
                   "Matrix P_quadForm must be " << numVariables << "x" << numVariables
                                                << ". Provided: " << P_.rows() << "x" << P_.cols());

        QL_REQUIRE(A_constraints.empty() ||
                       std::all_of(A_constraints.begin(), A_constraints.end(),
                                   [numVariables](const std::vector<double>& row) {
                                       return row.size() == numVariables;
                                   }),
                   "Matrix A_constraints must have " << numVariables
                                                     << " columns, but at least one row does not.");
        numConstraints_ = A_constraints.size();
        A_triplets_ = EigenUtilities::convertToTriplets(A_constraints);

        QL_REQUIRE(b_rhs.size() == numConstraints_,
                   "Vector b_rhs must have the same length as the number of rows in "
                   "A_constraints. Provided: "
                       << b_rhs.size() << ", expected: " << numConstraints_);
        b_list_ = std::vector<double>(b_rhs);
        b_ = Eigen::Map<Eigen::VectorXd>(b_list_.data(), b_list_.size());

        QL_REQUIRE(c_linearForm.empty() || c_linearForm.size() == numVariables,
                   "Vector c_linearForm must be empty or have length "
                       << numVariables << ". Provided: " << c_linearForm.size());

        c_list_ = c_linearForm.empty() ? std::vector<Real>(numVariables, 0.0) :
                                         std::vector<Real>(c_linearForm);
        c_ = Eigen::Map<Eigen::VectorXd>(c_list_.data(), c_list_.size());

        QL_REQUIRE(constraintTypes.size() == numConstraints_ || constraintTypes.empty(),
                   "Vector constraintTypes must be empty or have the same length as the number of "
                   "rows in A_constraints.\nProvided: "
                       << constraintTypes.size() << ", expected: " << numConstraints_);

        if (constraintTypes.empty()) {
            constraintTypes_ = std::vector<ConstraintType>(numConstraints_, ConstraintType::Equal);
            numEqualities_ = numConstraints_;
            numInequalities_ = 0;
            isOrdered_ = true;
        } else {
            constraintTypes_ = std::vector<ConstraintType>(constraintTypes);
            numEqualities_ =
                std::count(constraintTypes_.begin(), constraintTypes_.end(), ConstraintType::Equal);
            numInequalities_ = numConstraints_ - numEqualities_;
            // Check if constraints are already in SCS order
            bool inOrder = true;
            bool seenInequality = false;
            for (const auto& type : constraintTypes_) {
                if (type == ConstraintType::LessEqual) {
                    seenInequality = true;
                } else if (seenInequality) {
                    // Found equality after inequality - not in SCS order
                    inOrder = false;
                    break;
                }
            }
            isOrdered_ = inOrder;
            QL_REQUIRE(isOrdered_, "Constraints must be in SCS order (equalities first, then inequalities)");
        }

        parameters_ = Eigen::VectorXd::Zero(numParameters_);
        hasParameters_ = false;
    }

    SplineConstraints::SplineConstraints(Size numVariables,
                                         const Eigen::SparseMatrix<double>& P_quadForm,
                                         const std::vector<Eigen::Triplet<double>>& A_triplets,
                                         const std::vector<double>& b_rhs,
                                         const std::vector<double>& c_linearForm,
                                         const std::vector<ConstraintType>& constraintTypes,
                                         bool fitData,
                                         double epsAbsolute,
                                         double epsRelative,
                                         double epsInfeasible)
    : fitData_(fitData), numVariables_(numVariables), numEqualities_(0), numInequalities_(0),
      numParameters_(0), P_(P_quadForm), A_triplets_(A_triplets), b_list_(b_rhs),
      scsDataIsUpToDate_(false), isSolved_(false), epsAbsolute_(epsAbsolute), epsRelative_(epsRelative),
      epsInfeasible_(epsInfeasible) {

        QL_REQUIRE(static_cast<Size>(P_.rows()) == numVariables &&
                       static_cast<Size>(P_.cols()) == numVariables,
                   "Matrix P_quadForm must be " << numVariables << "x" << numVariables
                                                << ". Provided: " << P_.rows() << "x" << P_.cols());


        numConstraints_ = 0;
        // Find the highest row number, this is the number of constraints. The [&] capture in the
        // lambda makes all variables in the scope available by reference.
        std::for_each(
            A_triplets.begin(), A_triplets.end(), [&](const Eigen::Triplet<double>& triplet) {
                numConstraints_ = std::max(numConstraints_, static_cast<Size>(triplet.row()));
            });
        numConstraints_++;

        A_ = Eigen::SparseMatrix<double>(numConstraints_, numVariables);
        A_.setFromTriplets(A_triplets_.begin(), A_triplets_.end());

        QL_REQUIRE(b_rhs.size() == numConstraints_,
                   "Vector b_rhs must have the same length as the number of rows in "
                   "A_constraints. Provided: "
                       << b_rhs.size() << ", expected: " << numConstraints_);

        b_ = Eigen::Map<Eigen::VectorXd>(b_list_.data(), b_list_.size());

        QL_REQUIRE(c_.size() == 0 || static_cast<Size>(c_.size()) == numVariables,
                   "Vector c_linearForm must be empty or have length "
                       << numVariables << ". Provided: " << c_.size());

        c_ = Eigen::Map<const Eigen::VectorXd>(c_linearForm.data(), c_linearForm.size());
        if (c_.size() == 0) {
            c_ = Eigen::VectorXd::Zero(numVariables);
        }

        c_list_ = c_linearForm;

        QL_REQUIRE(constraintTypes.size() == numConstraints_ || constraintTypes.empty(),
                   "Vector constraintTypes must be empty or have the same length as the number of "
                   "rows in A_constraints.\nProvided: "
                       << constraintTypes.size() << ", expected: " << numConstraints_);

        if (constraintTypes.empty()) {
            constraintTypes_ = std::vector<ConstraintType>(numConstraints_, ConstraintType::Equal);
            numEqualities_ = numConstraints_;
            numInequalities_ = 0;
            isOrdered_ = true;
        } else {
            constraintTypes_ = std::vector<ConstraintType>(constraintTypes);
            numEqualities_ =
                std::count(constraintTypes_.begin(), constraintTypes_.end(), ConstraintType::Equal);
            numInequalities_ = std::count(constraintTypes_.begin(), constraintTypes_.end(),
                                          ConstraintType::LessEqual);
            // Check if constraints are already in SCS order
            bool inOrder = true;
            bool seenInequality = false;
            for (const auto& type : constraintTypes_) {
                if (type == ConstraintType::LessEqual) {
                    seenInequality = true;
                } else if (seenInequality) {
                    // Found equality after inequality - not in SCS order
                    inOrder = false;
                    break;
                }
            }
            isOrdered_ = inOrder;
            QL_REQUIRE(isOrdered_, "Constraints must be in SCS order (equalities first, then inequalities)");
        }

        parameters_ = Eigen::VectorXd::Zero(numParameters_);
        hasParameters_ = false;
    }

    // SplineConstraints::SplineConstraints(const SplineConstraints& other) {
    //     *this = other;
    // }

    // Removed updateOrdering - we maintain strict SCS ordering

    // Removed reorderByConstraints - we maintain strict SCS ordering

    // void SplineConstraints::setParameterMatrixB(
    //     Size nParameters, const std::vector<std::vector<double>>& B_parameterMatrix) {
    //     numParameters_ = nParameters;
    //     hasParameters_ = true;

    //    QL_REQUIRE(B_parameterMatrix.empty() ||
    //                   (B_parameterMatrix.size() == numConstraints_ &&
    //                    std::all_of(B_parameterMatrix.begin(), B_parameterMatrix.end(),
    //                                [nParameters](const std::vector<double>& row) {
    //                                    return row.size() == nParameters;
    //                                })),
    //               "Matrix B_parameterMatrix must be empty or be "
    //                   << numConstraints_ << "x" << nParameters
    //                   << ". Provided: " << B_parameterMatrix.size() << " rows.");

    //    B_triplets_ = EigenUtilities::convertToTriplets(B_parameterMatrix);
    //}

    void SplineConstraints::addParameters(Size nNewParameters,
                                          const Eigen::SparseMatrix<double>& B_new,
                                          const Eigen::SparseMatrix<double>& C_new) {
        numParameters_ += nNewParameters;
        parameters_list_.resize(numParameters_);
        parameters_ = Eigen::Map<Eigen::VectorXd>(
            parameters_list_.data(),
            numParameters_); // Zero initialized (Eigen::VectorXd::Zero(numParameters_)
        hasParameters_ = true;

        // QL_REQUIRE(static_cast<Size>(B_new.rows()) == numConstraints_ &&
        // static_cast<Size>(B_new.cols()) == nNewParameters,
        //            "Matrix B_new must be " << numConstraints_ << "x" << nNewParameters
        //                                    << ". Provided: " << B_new.rows() << "x"
        //                                    << B_new.cols());

        // QL_REQUIRE(static_cast<Size>(C_new.rows()) == numVariables_ &&
        // static_cast<Size>(C_new.cols()) == nNewParameters,
        //            "Matrix C_new must be " << numVariables_ << "x" << nNewParameters << ".
        //            Provided: "
        //                                    << C_new.rows() << "x" << C_new.cols());


        std::vector<Eigen::Triplet<double>> B_new_triplets =
            EigenUtilities::convertToTriplets(B_new);
        std::vector<Eigen::Triplet<double>> C_new_triplets =
            EigenUtilities::convertToTriplets(C_new);

        B_triplets_.insert(B_triplets_.end(), B_new_triplets.begin(), B_new_triplets.end());
        C_triplets_.insert(C_triplets_.end(), C_new_triplets.begin(), C_new_triplets.end());

        // Parameters don't affect ordering
        scsDataIsUpToDate_ = false;
    }

    void SplineConstraints::push() {
        constraintStack_.emplace(numConstraints_, numParameters_, A_triplets_.size(),
                                 B_triplets_.size(), C_triplets_.size(), numEqualities_,
                                 numInequalities_);
    }

    void SplineConstraints::pop() {
        if (!constraintStack_.empty()) {
            Size nConstraints = std::get<0>(constraintStack_.top());
            Size nParameters = std::get<1>(constraintStack_.top());
            Size nATriplets = std::get<2>(constraintStack_.top());
            Size nBTriplets = std::get<3>(constraintStack_.top());
            Size nCTriplets = std::get<4>(constraintStack_.top());
            numEqualities_ = std::get<5>(constraintStack_.top());
            numInequalities_ = std::get<6>(constraintStack_.top());
            constraintStack_.pop();

            b_list_.resize(nConstraints);
            b_.resize(nConstraints);
            parameters_list_.resize(nParameters);
            parameters_.resize(nParameters);
            A_triplets_.resize(nATriplets);
            B_triplets_.resize(nBTriplets);
            C_triplets_.resize(nCTriplets);
            numConstraints_ = nConstraints;
            numParameters_ = nParameters;
            // Pop doesn't affect SCS ordering
            scsDataIsUpToDate_ = false;
        }
    }

    Eigen::SparseMatrix<Real> SplineConstraints::getSliceOfA(Integer firstRow,
                                                             Integer lastRow) const {
        std::vector<Eigen::Triplet<double>> subTriplets;

        // Filter triplets for rows >= firstRow
        for (const auto& triplet : A_triplets_) {
            if (triplet.row() >= firstRow && triplet.row() < lastRow) {
                subTriplets.emplace_back(triplet.row() - firstRow, // Adjust row index
                                         triplet.col(), triplet.value());
            }
        }

        // Create the sub-matrix
        Eigen::SparseMatrix<double> subMatrix(lastRow - firstRow, numVariables_);
        subMatrix.setFromTriplets(subTriplets.begin(), subTriplets.end());
        return subMatrix;
    }

    // void SplineConstraints::setParameterMatrixC(
    //     Size nParameters, const std::vector<std::vector<double>>& C_parameterMatrix) {
    //     numParameters_ = nParameters;
    //     hasParameters_ = true;

    //    QL_REQUIRE(C_parameterMatrix.empty() ||
    //                   (C_parameterMatrix.size() == numVariables_ &&
    //                    std::all_of(C_parameterMatrix.begin(), C_parameterMatrix.end(),
    //                                [nParameters](const std::vector<double>& row) {
    //                                    return row.size() == nParameters;
    //                                })),
    //               "Matrix C_parameterMatrix must be empty or be "
    //                   << numVariables_ << "x" << nParameters
    //                   << ". Provided: " << C_parameterMatrix.size() << "rows.");

    //    C_triplets_ = EigenUtilities::convertToTriplets(C_parameterMatrix);
    //}

    void SplineConstraints::updateParameterValues(const Eigen::VectorXd& parameters) {
        parameters_ = parameters;
        scsDataIsUpToDate_ = false;
    }

    void SplineConstraints::updateParameterValues(const std::vector<Real>& parameters) {
        parameters_list_ = std::vector<Real>(parameters);
        parameters_ = Eigen::Map<Eigen::VectorXd>(parameters_list_.data(), parameters_list_.size());
        scsDataIsUpToDate_ = false;
    }

    void SplineConstraints::addLinearConstraint(const Eigen::VectorXd& constraint,
                                                double rhs,
                                                ConstraintType constraintType) {
        Size rowIndex = numConstraints_; // Use numConstraints_ to determine the row index

        // TODO: lazy programming, should do this without copying vector etc.
        Eigen::SparseVector<double> sparseRow = constraint.sparseView();
        for (Eigen::SparseVector<double>::InnerIterator it(sparseRow); it; ++it) {
            A_triplets_.emplace_back(rowIndex, it.index(), it.value());
        }
        b_list_.push_back(rhs);
        b_ = Eigen::Map<Eigen::VectorXd>(b_list_.data(), b_list_.size());
        constraintTypes_.push_back(constraintType);
        if (constraintType == ConstraintType::Equal) {
            numEqualities_++;
        } else {
            numInequalities_++;
        }
        numConstraints_++; // Increase the number of constraints
        QL_REQUIRE(numConstraints_ == numEqualities_ + numInequalities_,
                   "Wrong parity for equalities and inequalities");

        scsDataIsUpToDate_ = false;
        // Check if we're still in SCS order after adding constraint
        if (constraintType == ConstraintType::Equal && numInequalities_ > 0) {
            // Added equality after inequalities - order broken
            isOrdered_ = false;
            QL_FAIL("Cannot add equality constraint after inequality constraints. Use addEqualityConstraintAtBeginning() instead.");
        }
    }
    
    void SplineConstraints::addEqualityConstraintAtBeginning(const Eigen::VectorXd& constraint,
                                                              double rhs) {
        // To maintain SCS ordering, insert equality constraints before inequalities
        // New row index for this equality is at position numEqualities_ (after existing equalities)
        Size insertPosition = numEqualities_;
        
        // Shift all inequality constraint indices up by 1
        std::vector<Eigen::Triplet<Real>> newTriplets;
        newTriplets.reserve(A_triplets_.size() + constraint.nonZeros());
        
        // Copy equality constraints (rows 0 to numEqualities_-1)
        for (const auto& triplet : A_triplets_) {
            if (static_cast<Size>(triplet.row()) < numEqualities_) {
                newTriplets.push_back(triplet);
            } else {
                // Shift inequality constraints up by 1
                newTriplets.emplace_back(triplet.row() + 1, triplet.col(), triplet.value());
            }
        }
        
        // Add new equality constraint at position numEqualities_
        Eigen::SparseVector<double> sparseRow = constraint.sparseView();
        for (Eigen::SparseVector<double>::InnerIterator it(sparseRow); it; ++it) {
            newTriplets.emplace_back(insertPosition, it.index(), it.value());
        }
        
        // Update triplets
        A_triplets_ = std::move(newTriplets);
        
        // Insert RHS value at correct position
        b_list_.insert(b_list_.begin() + insertPosition, rhs);
        b_ = Eigen::Map<Eigen::VectorXd>(b_list_.data(), b_list_.size());
        
        // Insert constraint type at correct position
        constraintTypes_.insert(constraintTypes_.begin() + insertPosition, ConstraintType::Equal);
        
        // Update counts
        numEqualities_++;
        numConstraints_++;
        
        QL_REQUIRE(numConstraints_ == numEqualities_ + numInequalities_,
                   "Wrong parity for equalities and inequalities");
        
        // System remains ordered (SCS order maintained)
        isOrdered_ = true;
        scsDataIsUpToDate_ = false;
    }


    // void SplineConstraints::addLinearConstraint(const Eigen::VectorXd& constraint,
    //                                             bool isEquality,
    //                                             double rhs) {
    //     int rowIndex = numConstraints_; // Use numConstraints_ to determine the row index

    //    Eigen::SparseVector<double> sparseRow = constraint.sparseView();
    //    for (Eigen::SparseVector<double>::InnerIterator it(sparseRow); it; ++it) {
    //        A_triplets_.emplace_back(rowIndex, it.index(), it.value());
    //    }
    //    b_.push_back(rhs);
    //    numConstraints_++; // Increase the number of constraints

    //    scsDataIsUpToDate_ = false;
    //}

    // void SplineConstraints::addObjectiveFunction(const std::vector<std::vector<double>>& P,
    //                                              const std::vector<double>& c) {
    //     Size n = numVariables_;
    //     Eigen::SparseMatrix<double> P_eigen;
    //     if (!P.empty()) {
    //         P_eigen = EigenUtilities::convertToEigenSparseMatrix(P);
    //     } else {
    //         P_eigen = Eigen::SparseMatrix<double>(n, n); // zero matrix of correct dimensions
    //     }

    //    QL_ASSERT(P_.rows() == n && P_.cols() == n, "Matrix P_quadForm must " << n << "x" << n
    //                                                                 << ". Provided: " <<
    //                                                                 P_.rows()
    //                                                                 << "x" << P_.cols());

    //    QL_ASSERT(c.empty() || c.size() == n,
    //              "Vector c_linearForm must be empty or have length " << n << ". Provided: " <<
    //              c.size());

    //    Eigen::VectorXd c_eigen = !c.empty() ? EigenUtilities::convertToEigenVector(c) :
    //    Eigen::VectorXd::Zero(n);

    //    addObjectiveFunction(P_eigen, c_eigen);
    //}

    //// Empty and zero structures must be explicit, all dimensions need to check
    // void SplineConstraints::addObjectiveFunction(const Eigen::SparseMatrix<double>& P,
    //                                              const Eigen::VectorXd& c) {
    //     P_ += P;
    //     c_ += c;

    //    scsDataIsUpToDate_ = false;
    //}

    // I use these in the debugger sometimes TODO find a better solution
    [[maybe_unused]] static std::string matrixToString(
        const Eigen::SparseMatrix<double>& matrix) { // NOLINT(misc-use-anonymous-namespace)
        std::ostringstream oss;

        // Add the initial opening brace
        oss << "{\n";

        // Loop over the matrix and write its entries to the string
        for (int i = 0; i < matrix.rows(); ++i) {
            oss << "  {"; // Start each row with an opening brace

            for (int j = 0; j < matrix.cols(); ++j) {
                double value = matrix.coeff(i, j); // Get the value at (i, j)
                // oss << std::fixed << std::setprecision(1) << value;
                oss << value;

                // Add a comma and space unless it's the last column
                if (j != matrix.cols() - 1) {
                    oss << ", ";
                }
            }

            // End the row with a closing brace and a comma
            if (i < matrix.rows() - 1) {
                oss << "},\n";
            }
        }

        // Add the final closing brace
        oss << "}\n}";

        return oss.str(); // Convert the string stream to a string
    }

    void SplineConstraints::updateScsData() {
        // No reordering needed - always maintain SCS order
        A_.resize(numConstraints_, numVariables_);
        A_.setFromTriplets(A_triplets_.begin(), A_triplets_.end());

        if (hasParameters_) {
            B_.resize(numConstraints_, numParameters_);
            B_.setFromTriplets(B_triplets_.begin(), B_triplets_.end());
            C_.resize(numVariables_, numParameters_);
            C_.setFromTriplets(C_triplets_.begin(), C_triplets_.end());
            Eigen::VectorXd bp = b_ + B_ * parameters_;
            Eigen::VectorXd cp = c_ + C_ * parameters_;
            scsData_ = new SCS::SCSSolver(P_, A_, bp, cp, numEqualities_, numInequalities_,
                                          epsAbsolute_, epsRelative_, epsInfeasible_);
        } else {
            scsData_ = new SCS::SCSSolver(P_, A_, b_, c_, numEqualities_, numInequalities_,
                                          epsAbsolute_, epsRelative_, epsInfeasible_);
        }

        // TODO: Remove this, just for debug
        // matrixToString(A_);
        scsDataIsUpToDate_ = true;
        warmStart_ = 0;
    }

    int SplineConstraints::solve() {
        if (!scsDataIsUpToDate_) {
            updateScsData();
        }

        // Solve the SCS problem
        int status = scsData_->solve(warmStart_);
        warmStart_ = 1;

        if (status != 1) {
            // Don't print to stderr - it causes popups in some environments
            // The error will be handled by the exception thrown in the calling code
            // std::cerr << "Solver returned error: " << status << '\n';
        } else {
            isSolved_ = true;
        }

        return status;
    }

    Eigen::VectorXd SplineConstraints::getSolution() const {
        QL_REQUIRE(isSolved_, "No solution present, call solve() first.");
        return scsData_->solution_x();
    }

    Size SplineConstraints::getNumVariables() const {
        return numVariables_;
    }

    Size SplineConstraints::getNConstraints() const {
        return numConstraints_;
    }

    Size SplineConstraints::getNParameters() const {
        return numParameters_;
    }

    void SplineConstraints::update_b(const std::vector<Real>& parameters) {
        if (!scsDataIsUpToDate_) {
            updateScsData();
        }
        parameters_list_ = std::vector<double>(parameters);
        parameters_ = Eigen::Map<Eigen::VectorXd>(parameters_list_.data(), parameters_list_.size());

        Eigen::VectorXd bp = b_ + B_ * parameters_;
        Eigen::VectorXd cp = Eigen::VectorXd(0);
        scs_int status = scsData_->update(bp, cp);

        if (status != 0) {
            // Don't print to stderr - it causes popups in some environments
            // std::cerr << "Update returned error: " << status << "\n";
        }
    }

    void SplineConstraints::update_c(const std::vector<double>& parameters) {
        if (!scsDataIsUpToDate_) {
            updateScsData();
        }
        parameters_list_ = std::vector<double>(parameters);
        parameters_ = Eigen::Map<Eigen::VectorXd>(parameters_list_.data(), parameters_list_.size());

        Eigen::VectorXd bp = Eigen::VectorXd(0);
        Eigen::VectorXd cp = c_ + C_ * parameters_;
        scs_int status = scsData_->update(bp, cp);

        if (status != 0) {
            // Don't print to stderr - it causes popups in some environments
            // std::cerr << "Update returned error: " << status << "\n";
        }
    }
}
