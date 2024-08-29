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

#include "splineconstraints.hpp"
#include "eigenutilities.hpp"
#include <algorithm>
#include <iostream>
#include <malloc.h>
#include <numeric>
#include <vector>
#include <tuple>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <ql/errors.hpp>
#include <ql/types.hpp>
#include <scs.h>
#include <scs_types.h>

namespace QuantLib {
    SplineConstraints::SplineConstraints(
        Size nVariables,
        const std::vector<std::vector<double>>& P_quadForm,
        const std::vector<std::vector<double>>& A_constraints,
        const std::vector<double>& b_rhs,
        const std::vector<double>& c_linearForm,
        const std::vector<ConstraintType>& constraintTypes
        // Size nParameters,
        // const std::vector<double>& parameters,
        // const std::vector<std::vector<double>>& B_parameterMatrix,
        // const std::vector<std::vector<double>>& C_parameterMatrix)
        )
    : nVariables_(nVariables), nParameters_(0), scsDataIsUpToDate_(false), isSolved_(false), nEqualities_(0), nInequalities_(0) {
        if (!P_quadForm.empty()) {
            P_ = EigenUtilities::convertToEigenSparseMatrix(P_quadForm);
        } else {
            P_ = Eigen::SparseMatrix<double>(nVariables, nVariables); // zero matrix
        }

        QL_REQUIRE(P_.rows() == nVariables && P_.cols() == nVariables,
                   "Matrix P_quadForm must be " << nVariables << "x" << nVariables
                                                << ". Provided: " << P_.rows() << "x" << P_.cols());

        QL_REQUIRE(A_constraints.empty() ||
                       std::all_of(A_constraints.begin(), A_constraints.end(),
                                   [nVariables](const std::vector<double>& row) {
                                       return row.size() == nVariables;
                                   }),
                   "Matrix A_constraints must have " << nVariables
                                                     << " columns, but at least one row does not.");
        nConstraints_ = A_constraints.size();
        A_triplets_ = EigenUtilities::convertToTriplets(A_constraints);

        QL_REQUIRE(b_rhs.size() == nConstraints_,
                   "Vector b_rhs must have the same length as the number of rows in "
                   "A_constraints. Provided: "
                       << b_rhs.size() << ", expected: " << nConstraints_);
        b_list_ = std::vector<double>(b_rhs);
        b_ = Eigen::Map<Eigen::VectorXd>(b_list_.data(), b_list_.size());

        QL_REQUIRE(c_linearForm.empty() || c_linearForm.size() == nVariables,
                   "Vector c_linearForm must be empty or have length "
                       << nVariables << ". Provided: " << c_linearForm.size());

        c_list_ = c_linearForm.empty() ? std::vector<double>(nVariables, 0.0) :
                                         std::vector<double>(c_linearForm);
        c_ = Eigen::Map<Eigen::VectorXd>(c_list_.data(), c_list_.size());

        QL_REQUIRE(constraintTypes.size() == nConstraints_ || constraintTypes.empty(),
                   "Vector constraintTypes must be empty or have the same length as the number of "
                   "rows in A_constraints.\nProvided: "
                       << constraintTypes.size() << ", expected: " << nConstraints_);

        if (constraintTypes.empty()) {
            constraintTypes_ = std::vector<ConstraintType>(nConstraints_, ConstraintType::Equal);
            nEqualities_ = nConstraints_;
            nInequalities_ = 0;
            isOrdered_ = true;
        } else {
            constraintTypes_ = std::vector<ConstraintType>(constraintTypes);
            nEqualities_ = std::count(constraintTypes_.begin(), constraintTypes_.end(),
									  ConstraintType::Equal);
            nInequalities_ = nConstraints_ - nEqualities_;
            isOrdered_ = false;
        }

        parameters_ = Eigen::VectorXd::Zero(nParameters_);
        hasParameters_ = false;
    }

    SplineConstraints::SplineConstraints(const SplineConstraints& other) {
        *this = other;
    }

    void SplineConstraints::updateOrdering() {
        // Create trivial permutation
        permutation_ = std::vector<int>(nConstraints_);
        std::iota(permutation_.begin(), permutation_.end(), 0);

        // Sort permutation stably using the constraint type enum
        std::stable_sort(permutation_.begin(), permutation_.end(), [this](int i, int j) {
            return constraintTypes_[i] < constraintTypes_[j];
        });
    }

    void SplineConstraints::reorderByConstraints() {
        // transform the triplets renumbering the rows according to permutation, both for A and B
        // matrices
        updateOrdering();
        std::vector<Eigen::Triplet<double>> reorderedTriplets;
        reorderedTriplets.reserve(A_triplets_.size());
        for (const auto& triplet : A_triplets_) {
            reorderedTriplets.emplace_back(permutation_[triplet.row()], triplet.col(),
                                           triplet.value());
        }
        A_.resize(nConstraints_, nVariables_);
        A_.setFromTriplets(reorderedTriplets.begin(), reorderedTriplets.end());

        if (hasParameters_) {
            reorderedTriplets.clear();
            reorderedTriplets.reserve(B_triplets_.size());
            for (const auto& triplet : B_triplets_) {
                reorderedTriplets.emplace_back(permutation_[triplet.row()], triplet.col(),
                                               triplet.value());
            }
            B_.resize(nConstraints_, nParameters_);
            B_.setFromTriplets(reorderedTriplets.begin(), reorderedTriplets.end());
        }

        //auto begin = PermutationIterator::make_permutation_iterator(A_triplets_.begin(), permutation_);
        //auto end = PermutationIterator::make_permutation_iterator(A_triplets_.end(), permutation_);
        //A_.setFromTriplets(begin, end);

        //if (hasParameters_) {
        //    auto B_begin = make_permutation_iterator(B_triplets_.begin(), permutation_);
        //    auto B_end = make_permutation_iterator(B_triplets_.end(), permutation_);
        //    B_.setFromTriplets(B_begin, B_end);
        //}
        // Reorder rhs values and constraints
        std::vector<double> reorderedB(nConstraints_);
        //std::vector<ConstraintType> reorderedConstraintTypes(nConstraints_);
        nEqualities_ = nInequalities_ = 0;
        for (int i = 0; i < nConstraints_; ++i) {
            reorderedB[i] = b_list_[permutation_[i]];
            //reorderedConstraintTypes[i] = constraintTypes_[permutation_[i]];
            if (constraintTypes_[i] == ConstraintType::Equal) {
                nEqualities_++;
            } else {
                nInequalities_++;
            }
        }
        b_list_ = reorderedB;
        b_ = Eigen::Map<Eigen::VectorXd>(b_list_.data(), nConstraints_);

        //constraintTypes_ = reorderedConstraintTypes;

        //QL_ASSERT(
        //    (nEqualities_ == 0 || constraintTypes_[nEqualities_ - 1] == ConstraintType::Equal) &&
        //        (nInequalities_ == 0 ||
        //         constraintTypes_[nEqualities_] == ConstraintType::LessEqual),
        //    "Something went wrong with reordering constraints.");
        isOrdered_ = true;
    }

    //void SplineConstraints::setParameterMatrixB(
    //    Size nParameters, const std::vector<std::vector<double>>& B_parameterMatrix) {
    //    nParameters_ = nParameters;
    //    hasParameters_ = true;

    //    QL_REQUIRE(B_parameterMatrix.empty() ||
    //                   (B_parameterMatrix.size() == nConstraints_ &&
    //                    std::all_of(B_parameterMatrix.begin(), B_parameterMatrix.end(),
    //                                [nParameters](const std::vector<double>& row) {
    //                                    return row.size() == nParameters;
    //                                })),
    //               "Matrix B_parameterMatrix must be empty or be "
    //                   << nConstraints_ << "x" << nParameters
    //                   << ". Provided: " << B_parameterMatrix.size() << " rows.");

    //    B_triplets_ = EigenUtilities::convertToTriplets(B_parameterMatrix);
    //}

    void SplineConstraints::addParameters(Size nNewParameters,
                                          Eigen::SparseMatrix<double>& B_new,
                                          Eigen::SparseMatrix<double>& C_new) {
        nParameters_ += nNewParameters;
        parameters_list_.resize(nParameters_);
        parameters_ = Eigen::Map<Eigen::VectorXd>(parameters_list_.data(), nParameters_); // Zero initialized (Eigen::VectorXd::Zero(nParameters_)
        hasParameters_ = true;

        QL_REQUIRE(B_new.rows() == nConstraints_ && B_new.cols() == nNewParameters,
                   "Matrix B_new must be " << nConstraints_ << "x" << nNewParameters
                                           << ". Provided: " << B_new.rows() << "x"
                                           << B_new.cols());

        QL_REQUIRE(C_new.rows() == nVariables_ && C_new.cols() == nNewParameters,
                   "Matrix C_new must be " << nVariables_ << "x" << nNewParameters << ". Provided: "
                                           << C_new.rows() << "x" << C_new.cols());


        std::vector<Eigen::Triplet<double>> B_new_triplets =
            EigenUtilities::convertToTriplets(B_new);
        std::vector<Eigen::Triplet<double>> C_new_triplets =
            EigenUtilities::convertToTriplets(C_new);

        B_triplets_.insert(B_triplets_.end(), B_new_triplets.begin(), B_new_triplets.end());
        C_triplets_.insert(C_triplets_.end(), C_new_triplets.begin(), C_new_triplets.end());

        isOrdered_ = false;
        scsDataIsUpToDate_ = false;
    }

    void SplineConstraints::push() {
        constraintStack_.push(std::make_tuple(nConstraints_, nParameters_, A_triplets_.size(),
                                              B_triplets_.size(), C_triplets_.size()));
    }

    void SplineConstraints::pop() {
        if (!constraintStack_.empty()) {
            Size nConstraints = std::get<0>(constraintStack_.top());
            Size nParameters = std::get<1>(constraintStack_.top());
            Size nATriplets = std::get<2>(constraintStack_.top());
            Size nBTriplets = std::get<3>(constraintStack_.top());
            Size nCTriplets = std::get<4>(constraintStack_.top());
            constraintStack_.pop();

            if (nConstraints <= nConstraints_) {
                b_list_.resize(nConstraints);
                b_.resize(nConstraints);
            }
            if (nParameters <= nParameters_) {
                parameters_list_.resize(nParameters);
                parameters_.resize(nParameters);
            }
            if (nATriplets <= A_triplets_.size()) {
                A_triplets_.resize(nATriplets);
            }
            if (nBTriplets <= B_triplets_.size()) {
                B_triplets_.resize(nBTriplets);
            }
            if (nCTriplets <= C_triplets_.size()) {
                C_triplets_.resize(nCTriplets);
            }
            nConstraints_ = nConstraints;
            nParameters_ = nParameters;
            isOrdered_ = false;
            scsDataIsUpToDate_ = false;
        }
    }

    //void SplineConstraints::setParameterMatrixC(
    //    Size nParameters, const std::vector<std::vector<double>>& C_parameterMatrix) {
    //    nParameters_ = nParameters;
    //    hasParameters_ = true;

    //    QL_REQUIRE(C_parameterMatrix.empty() ||
    //                   (C_parameterMatrix.size() == nVariables_ &&
    //                    std::all_of(C_parameterMatrix.begin(), C_parameterMatrix.end(),
    //                                [nParameters](const std::vector<double>& row) {
    //                                    return row.size() == nParameters;
    //                                })),
    //               "Matrix C_parameterMatrix must be empty or be "
    //                   << nVariables_ << "x" << nParameters
    //                   << ". Provided: " << C_parameterMatrix.size() << "rows.");

    //    C_triplets_ = EigenUtilities::convertToTriplets(C_parameterMatrix);
    //}

    void SplineConstraints::updateParameterValues(const Eigen::VectorXd& parameters) {
        parameters_ = parameters;
        scsDataIsUpToDate_ = false;
    }

    void SplineConstraints::updateParameterValues(const std::vector<double> parameters) {
        parameters_list_ = std::vector<double>(parameters);
        parameters_ = Eigen::Map<Eigen::VectorXd>(parameters_list_.data(), parameters_list_.size());
        scsDataIsUpToDate_ = false;
    }

    void SplineConstraints::addLinearConstraint(const Eigen::VectorXd& constraint,
                                                double rhs,
                                                ConstraintType constraintType) {
        int rowIndex = nConstraints_; // Use nConstraints_ to determine the row index

        // TODO lazy programming, should do this without copying vector etc.
        Eigen::SparseVector<double> sparseRow = constraint.sparseView();
        for (Eigen::SparseVector<double>::InnerIterator it(sparseRow); it; ++it) {
            A_triplets_.emplace_back(rowIndex, it.index(), it.value());
        }
        b_list_.push_back(rhs);
        b_ = Eigen::Map<Eigen::VectorXd>(b_list_.data(), b_list_.size());
        constraintTypes_.push_back(constraintType);
        nConstraints_++; // Increase the number of constraints

        scsDataIsUpToDate_ = false;
        isOrdered_ = false;
    }


    // void SplineConstraints::addLinearConstraint(const Eigen::VectorXd& constraint,
    //                                             bool isEquality,
    //                                             double rhs) {
    //     int rowIndex = nConstraints_; // Use nConstraints_ to determine the row index

    //    Eigen::SparseVector<double> sparseRow = constraint.sparseView();
    //    for (Eigen::SparseVector<double>::InnerIterator it(sparseRow); it; ++it) {
    //        A_triplets_.emplace_back(rowIndex, it.index(), it.value());
    //    }
    //    b_.push_back(rhs);
    //    nConstraints_++; // Increase the number of constraints

    //    scsDataIsUpToDate_ = false;
    //}

    // void SplineConstraints::addObjectiveFunction(const std::vector<std::vector<double>>& P,
    //                                              const std::vector<double>& c) {
    //     Size n = nVariables_;
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

    //// Empty and zero structures must be explict, all dimensions need to check
    // void SplineConstraints::addObjectiveFunction(const Eigen::SparseMatrix<double>& P,
    //                                              const Eigen::VectorXd& c) {
    //     P_ += P;
    //     c_ += c;

    //    scsDataIsUpToDate_ = false;
    //}

    void SplineConstraints::updateScsData() {
        if (!isOrdered_) {
            reorderByConstraints();
        } else {
            A_.resize(nConstraints_, nVariables_);
            A_.setFromTriplets(A_triplets_.begin(), A_triplets_.end());
        }

        if (hasParameters_) {
            B_.resize(nConstraints_, nParameters_);
            B_.setFromTriplets(B_triplets_.begin(), B_triplets_.end());
            C_.resize(nVariables_, nParameters_);
            C_.setFromTriplets(C_triplets_.begin(), C_triplets_.end());
            Eigen::VectorXd bp = b_ + B_ * parameters_;
            Eigen::VectorXd cp = c_ + C_ * parameters_;
            //scsData_ = SCS::SCSSolver(P_, A_, bp, cp, nEqualities_, nInequalities_); 
            scsData_ = new SCS::SCSSolver(P_, A_, bp, cp, nEqualities_, nInequalities_);
        } else {
            //scsData_ = SCS::SCSSolver(P_, A_, b_, c_, nEqualities_, nInequalities_);
            scsData_ = new SCS::SCSSolver(P_, A_, b_, c_, nEqualities_, nInequalities_);
        }

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
            std::cerr << "Solver returned error: " << status << std::endl;
        } else {
            isSolved_ = true;
        }

        return status;
    }

    Eigen::VectorXd SplineConstraints::getSolution() {
        QL_REQUIRE(isSolved_, "No solution present, call solve() first.");
        return scsData_->solution_x();
    }

    Size SplineConstraints::getNVariables() const {
        return Size(nVariables_);
    }

    Size SplineConstraints::getNConstraints() const {
        return Size(nConstraints_);
    }

    Size SplineConstraints::getNParameters() const {
        return Size(nParameters_);
    }

    void SplineConstraints::update_b(const std::vector<double> params) {
        if (!scsDataIsUpToDate_) {
            updateScsData();
        }
        parameters_list_ = std::vector<double>(params);
        parameters_ = Eigen::Map<Eigen::VectorXd>(parameters_list_.data(), parameters_list_.size());

        Eigen::VectorXd bp = b_ + B_ * parameters_;
        Eigen::VectorXd cp = Eigen::VectorXd(0);
        scs_int status = scsData_->update(bp, cp);

        if (status != 0) {
            std::cerr << "Update returned error: " << status << std::endl;
        }
    }

    void SplineConstraints::update_c(std::vector<double> params) {
        if (!scsDataIsUpToDate_) {
            updateScsData();
        }
        parameters_list_ = std::vector<double>(params);
        parameters_ = Eigen::Map<Eigen::VectorXd>(parameters_list_.data(), parameters_list_.size());

        Eigen::VectorXd bp = Eigen::VectorXd(0);
        Eigen::VectorXd cp = c_ + C_ * parameters_;
        scs_int status = scsData_->update(bp, cp);

        if (status != 0) {
            std::cerr << "Update returned error: " << status << std::endl;
        }
    }
}


// Note that in the SCS structure we use pointers into std::vector-s, careful not to destruct these!
// TODO fix this, safer to just copy.
namespace SCS {
    // Implementing the constructors and destructor
    SCSSolver::SCSSolver()
    : nVariables_(0), nConstraints_(0), nEqualities_(0), nInequalities_(0), scs_data_(),
      scs_cone_(), scs_sol_(), scs_settings_(), scs_info_(), scs_work_(nullptr),
      is_updated_b_but_data_not_reset_(false), is_updated_c_but_data_not_reset_(false) {
        scs_data_.P = nullptr;
        scs_data_.A = nullptr;
    }

    SCSSolver::SCSSolver(const Eigen::SparseMatrix<double>& P_quadForm,
                         const Eigen::SparseMatrix<double>& A_constraints,
                         const Eigen::VectorXd& bp_rhs,
                         const Eigen::VectorXd& cp_linearForm,
                         QuantLib::Integer nEqualities,
                         QuantLib::Integer nInequalities) {
        convertEigenToCSC(P_quadForm.triangularView<Eigen::Upper>(), scs_P_x_, scs_P_i_, scs_P_p_);
        convertEigenToCSC(A_constraints, scs_A_x_, scs_A_i_, scs_A_p_);
        scs_b_ = bp_rhs;
        scs_c_ = cp_linearForm;
        nVariables_ = static_cast<scs_int>(P_quadForm.rows());
        nConstraints_ = static_cast<scs_int>(bp_rhs.size());
        nEqualities_ = static_cast<scs_int>(nEqualities);
        nInequalities_ = static_cast<scs_int>(nInequalities);
        initializeSCS();
    }

    void SCSSolver::convertEigenToCSC(const Eigen::SparseMatrix<double>& mat,
                                      std::vector<scs_float>& scs_x,
                                      std::vector<scs_int>& scs_i,
                                      std::vector<scs_int>& scs_p) {
        scs_x.clear();
        scs_i.clear();
        scs_p.clear();

        for (int k = 0; k < mat.outerSize(); ++k) {
            scs_p.push_back(static_cast<scs_int>(scs_x.size()));
            for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it) {
                scs_x.push_back(it.value());
                scs_i.push_back(static_cast<scs_int>(it.row()));
            }
        }
        scs_p.push_back(static_cast<scs_int>(scs_x.size()));
    }

    void SCSSolver::initializeSCS() {
        scs_data_.m = nConstraints_;
        scs_data_.n = nVariables_;
        scs_data_.P = createSCSMatrix(scs_data_.n, scs_data_.n, scs_P_x_.size(), scs_P_x_.data(),
                                      scs_P_i_.data(), scs_P_p_.data());
        scs_data_.A = createSCSMatrix(scs_data_.m, scs_data_.n, scs_A_x_.size(), scs_A_x_.data(),
                                      scs_A_i_.data(), scs_A_p_.data());
        scs_data_.b = createCVector(scs_b_);
        scs_data_.c = createCVector(scs_c_);

        scs_cone_.z = nEqualities_;   /* number of linear equalities */
        scs_cone_.l = nInequalities_; /* number of linear inequalities */
        // No other cones for the time being
        scs_cone_.bu = NULL;
        scs_cone_.bl = NULL;
        scs_cone_.bsize = 0;
        scs_cone_.q = NULL;
        scs_cone_.qsize = 0;
        scs_cone_.s = NULL;
        scs_cone_.ssize = 0;
        scs_cone_.ep = 0;
        scs_cone_.ed = 0;
        scs_cone_.p = NULL;
        scs_cone_.psize = 0;

        // initialize these with zero these for safety
        scs_sol_.x = (scs_float*)calloc(nVariables_, sizeof(scs_float));
        scs_sol_.y = (scs_float*)calloc(nConstraints_, sizeof(scs_float));
        scs_sol_.s = (scs_float*)calloc(
            nConstraints_,
            sizeof(scs_float)); // Seems to be for all constraints, not just inequality constraints

        scs_set_default_settings(&scs_settings_);
        // Disable verbose output
        scs_settings_.verbose = 0;

        scs_work_ = scs_init(&scs_data_, &scs_cone_, &scs_settings_);

        QL_REQUIRE(scs_work_ != nullptr, "Initialization of conic solver failed.");
    }

    QuantLib::Integer SCS::SCSSolver::solve() {
        return static_cast<QuantLib::Integer>(scs_solve(scs_work_, &scs_sol_, &scs_info_,
                                                        scs_settings_.warm_start));
    }

    QuantLib::Integer SCS::SCSSolver::solve(int warm_start) {
        return static_cast<QuantLib::Integer>(
            scs_solve(scs_work_, &scs_sol_, &scs_info_, warm_start));
    }

    QuantLib::Integer SCS::SCSSolver::update(Eigen::VectorXd& b, Eigen::VectorXd& c) {
        return static_cast<QuantLib::Integer>(scs_update(
            scs_work_, b.size() == 0 ? SCS_NULL : b.data(), c.size() == 0 ? SCS_NULL : c.data()));
    }

    Eigen::VectorXd SCS::SCSSolver::solution_x() {
        std::vector<double> sol(scs_sol_.x, scs_sol_.x + nVariables_);
        return Eigen::Map<Eigen::VectorXd>(sol.data(), nVariables_);
    }

    Eigen::VectorXd SCS::SCSSolver::solution_y() {
        std::vector<double> sol(scs_sol_.y, scs_sol_.y + nConstraints_);
        return Eigen::Map<Eigen::VectorXd>(sol.data(), nConstraints_);
    }

    Eigen::VectorXd SCS::SCSSolver::solution_s() {
        std::vector<double> sol(scs_sol_.s, scs_sol_.s + nConstraints_);
        return Eigen::Map<Eigen::VectorXd>(sol.data(), nConstraints_);
    }

    SCSSolver::~SCSSolver() {
        // scs structs other than scs_work_ are on the stack
        scs_finish(scs_work_);
        free(scs_data_.P);
        free(scs_data_.A);
        free(scs_data_.b);
        free(scs_data_.c);
        free(scs_sol_.s);
        free(scs_sol_.x);
        free(scs_sol_.y);
    }

    ScsMatrix* SCSSolver::createSCSMatrix(
        scs_int rows, scs_int cols, scs_int nnz, scs_float* x, scs_int* i, scs_int* p) {
        ScsMatrix* mat = (ScsMatrix*)malloc(sizeof(ScsMatrix));
        QL_ASSERT(mat, "Memory allocation failed for SCS matrix.");

        mat->m = rows;
        mat->n = cols;
        mat->x = x;
        mat->i = i;
        mat->p = p;
        return mat;
    }

    scs_float* SCSSolver::createCVector(const Eigen::VectorXd& src) {
        scs_float* dest = (scs_float*)malloc(src.size() * sizeof(scs_float));
        QL_ASSERT(dest, "Memory allocation failed for SCS vector.");

        std::copy(src.begin(), src.end(), dest);
        return dest;
    }
}