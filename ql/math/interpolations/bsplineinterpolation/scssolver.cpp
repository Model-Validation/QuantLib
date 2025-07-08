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

using namespace QuantLib;

// Note that in the SCS structure we use pointers into std::vector-s, careful not to destruct these!
// TODO: fix this, safer to just copy.
namespace SCS {
    // Implementing the constructors and destructor
    SCSSolver::SCSSolver()
    : nVariables_(0), nConstraints_(0), nEqualities_(0), nInequalities_(0), scs_data_(),
      scs_cone_(), scs_settings_(), scs_info_(), scs_work_(nullptr), scs_sol_() {
        scs_data_.P = nullptr;
        scs_data_.A = nullptr;
    }

    SCSSolver::SCSSolver(const Eigen::SparseMatrix<double>& P_quadForm,
                         const Eigen::SparseMatrix<double>& A_constraints,
                         const Eigen::VectorXd& bp_rhs,
                         const Eigen::VectorXd& cp_linearForm,
                         Size nEqualities,
                         Size nInequalities) {
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

    // ReSharper disable once CppMemberFunctionMayBeStatic
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
        scs_data_.P =
            createSCSMatrix(scs_data_.n, scs_data_.n, static_cast<scs_int>(scs_P_x_.size()),
                            scs_P_x_.data(), scs_P_i_.data(), scs_P_p_.data());
        scs_data_.A =
            createSCSMatrix(scs_data_.m, scs_data_.n, static_cast<scs_int>(scs_A_x_.size()),
                            scs_A_x_.data(), scs_A_i_.data(), scs_A_p_.data());
        scs_data_.b = createCVector(scs_b_);
        scs_data_.c = createCVector(scs_c_);

        scs_cone_.z = nEqualities_;   /* number of linear equalities */
        scs_cone_.l = nInequalities_; /* number of linear inequalities */
        // No other cones for the time being
        scs_cone_.bu = nullptr;
        scs_cone_.bl = nullptr;
        scs_cone_.bsize = 0;
        scs_cone_.q = nullptr;
        scs_cone_.qsize = 0;
        scs_cone_.s = nullptr;
        scs_cone_.ssize = 0;
        scs_cone_.ep = 0;
        scs_cone_.ed = 0;
        scs_cone_.p = nullptr;
        scs_cone_.psize = 0;

        // initialize these with zero these for safety
        scs_sol_.x = static_cast<scs_float*>(calloc(nVariables_, sizeof(scs_float)));
        scs_sol_.y = static_cast<scs_float*>(calloc(nConstraints_, sizeof(scs_float)));
        scs_sol_.s = static_cast<scs_float*>(calloc(
            nConstraints_,
            sizeof(scs_float))); // Seems to be for all constraints, not just inequality constraints

        scs_set_default_settings(&scs_settings_);
        // Disable verbose output
        scs_settings_.verbose = 0;
        // TODO: One should be able to set these
        scs_settings_.eps_abs = 1e-12;
        scs_settings_.eps_rel = 1e-12;
        scs_settings_.eps_infeas = 1e-13;

        scs_work_ = scs_init(&scs_data_, &scs_cone_, &scs_settings_);

        QL_REQUIRE(scs_work_ != nullptr, "Initialization of conic solver failed.");
    }

    Integer SCS::SCSSolver::solve() {
        return static_cast<Integer>(
            scs_solve(scs_work_, &scs_sol_, &scs_info_, scs_settings_.warm_start));
    }

    Integer SCS::SCSSolver::solve(int warm_start) {
        return static_cast<Integer>(
            scs_solve(scs_work_, &scs_sol_, &scs_info_, warm_start));
    }

    Integer SCS::SCSSolver::update(Eigen::VectorXd& b, Eigen::VectorXd& c) const {
        return static_cast<Integer>(scs_update(
            scs_work_, b.size() == 0 ? nullptr : b.data(), c.size() == 0 ? nullptr : c.data()));
    }

    Eigen::VectorXd SCS::SCSSolver::solution_x() const {
        std::vector<double> sol(scs_sol_.x, scs_sol_.x + nVariables_);
        return Eigen::Map<Eigen::VectorXd>(sol.data(), nVariables_);
    }

    Eigen::VectorXd SCS::SCSSolver::solution_y() const {
        std::vector<double> sol(scs_sol_.y, scs_sol_.y + nConstraints_);
        return Eigen::Map<Eigen::VectorXd>(sol.data(), nConstraints_);
    }

    Eigen::VectorXd SCS::SCSSolver::solution_s() const {
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

    // ReSharper disable once CppMemberFunctionMayBeStatic
    ScsMatrix* SCSSolver::createSCSMatrix(
        scs_int rows, scs_int cols, scs_int nnz, scs_float* x, scs_int* i, scs_int* p) {
        ScsMatrix* mat = static_cast<ScsMatrix*>(malloc(sizeof(ScsMatrix)));
        QL_ASSERT(mat, "Memory allocation failed for SCS matrix.");

        mat->m = rows;
        mat->n = cols;
        mat->x = x;
        mat->i = i;
        mat->p = p;
        return mat;
    }

    // ReSharper disable once CppMemberFunctionMayBeStatic
    scs_float* SCSSolver::createCVector(const Eigen::VectorXd& src) {
        scs_float* dest = static_cast<scs_float*>(malloc(src.size() * sizeof(scs_float)));
        QL_ASSERT(dest, "Memory allocation failed for SCS vector.");

        std::copy(src.begin(), src.end(), dest);
        return dest;
    }
}