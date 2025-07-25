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

/*! \file scssolver.hpp
    \brief SCS solver integration with Eigen.
*/
// ReSharper disable CppInconsistentNaming
#ifndef scs_solver_hpp
#define scs_solver_hpp

#include <ql/types.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <scs.h>
#include <scs_types.h>

namespace SCS {

    /*!
     * \brief Class for integrating the SCS solver with Eigen.
     */
    class SCSSolver {  // NOLINT(cppcoreguidelines-special-member-functions)
      public:
        /*!
         * \brief Default constructor for SCSSolver.
         */
        SCSSolver();

        /*!
         * \brief Constructor for SCSSolver with specified parameters.
         * \param P_quadForm Quadratic form matrix.
         * \param A_constraints Constraint matrix.
         * \param bp_rhs Right-hand side vector.
         * \param cp_linearForm Linear form vector.
         * \param nEqualities Number of equality constraints.
         * \param nInequalities Number of inequality constraints.
         * \param eps_abs Absolute tolerance for convergence (default: 1e-12).
         * \param eps_rel Relative tolerance for convergence (default: 1e-12).
         * \param eps_infeas Infeasibility tolerance (default: 1e-13).
         */
        SCSSolver(const Eigen::SparseMatrix<double>& P_quadForm,
                  const Eigen::SparseMatrix<double>& A_constraints,
                  const Eigen::VectorXd& bp_rhs,
                  const Eigen::VectorXd& cp_linearForm,
                  QuantLib::Size nEqualities,
                  QuantLib::Size nInequalities,
                  double eps_abs = 1e-12,
                  double eps_rel = 1e-12,
                  double eps_infeas = 1e-13);

        /*!
         * \brief Destructor for SCSSolver.
         */
        ~SCSSolver();

        /*!
         * \brief Convert Eigen structures to SCS format.
         * \param mat The Eigen sparse matrix to convert.
         * \param scs_x Vector to store the SCS values.
         * \param scs_i Vector to store the SCS row indices.
         * \param scs_p Vector to store the SCS column pointers.
         */
        void convertEigenToCSC(const Eigen::SparseMatrix<double>& mat,
                               std::vector<scs_float>& scs_x,
                               std::vector<scs_int>& scs_i,
                               std::vector<scs_int>& scs_p);

        /*!
         * \brief Initialize the SCS solver.
         */
        void initializeSCS();

        /*!
         * \brief Solve the optimization problem.
         * \return The status of the solver.
         */
        QuantLib::Integer solve();

        /*!
         * \brief Solve the optimization problem with warm start.
         * \param warm_start Warm start parameter.
         * \return The status of the solver.
         */
        QuantLib::Integer solve(int warm_start);

        /*!
         * \brief Update the right-hand side and linear form vectors.
         * \param b The new right-hand side vector.
         * \param c The new linear form vector.
         * \return The status of the update.
         */
        QuantLib::Integer update(Eigen::VectorXd& b, Eigen::VectorXd& c) const;

        /*!
         * \brief Get the solution vector x.
         * \return The solution vector x.
         */
        Eigen::VectorXd solution_x() const;

        /*!
         * \brief Get the solution vector y.
         * \return The solution vector y.
         */
        Eigen::VectorXd solution_y() const;

        /*!
         * \brief Get the solution vector s.
         * \return The solution vector s.
         */
        Eigen::VectorXd solution_s() const;

      private:
        // SCS data structures
        std::vector<scs_float> scs_P_x_; /*!< SCS values for matrix P. */
        std::vector<scs_int> scs_P_i_;   /*!< SCS row indices for matrix P. */
        std::vector<scs_int> scs_P_p_;   /*!< SCS column pointers for matrix P. */
        std::vector<scs_float> scs_A_x_; /*!< SCS values for matrix A. */
        std::vector<scs_int> scs_A_i_;   /*!< SCS row indices for matrix A. */
        std::vector<scs_int> scs_A_p_;   /*!< SCS column pointers for matrix A. */
        Eigen::VectorXd scs_b_;          /*!< SCS right-hand side vector. */
        Eigen::VectorXd scs_c_;          /*!< SCS linear form vector. */

        scs_int nVariables_;    /*!< Number of variables. */
        scs_int nConstraints_;  /*!< Number of constraints. */
        scs_int nEqualities_;   /*!< Number of equality constraints. */
        scs_int nInequalities_; /*!< Number of inequality constraints. */

        // Add tolerances to constructor
        double epsAbsolute_ = 1e-12;
        double epsRelative_ = 1e-12;
        double epsInfeasible_ = 1e-13;

        //bool is_updated_b_but_data_not_reset_; /*!< Flag indicating if b is updated but data not
        //                                          reset. */
        //bool is_updated_c_but_data_not_reset_; /*!< Flag indicating if c is updated but data not
        //                                          reset. */

        ScsData scs_data_;         /*!< SCS data structure. */
        ScsCone scs_cone_;         /*!< SCS cone structure. */
        ScsSettings scs_settings_; /*!< SCS settings structure. */
        ScsInfo scs_info_;         /*!< SCS info structure. */
        ScsWork* scs_work_;        /*!< SCS work structure. */
        ScsSolution scs_sol_;      /*!< SCS solution structure. */

        /*!
         * \brief Create an SCS matrix from given data.
         * \param rows Number of rows.
         * \param cols Number of columns.
         * \param nnz Number of non-zero elements.
         * \param x Values of the matrix.
         * \param i Row indices of the matrix.
         * \param p Column pointers of the matrix.
         * \return Pointer to the created SCS matrix.
         */
        ScsMatrix* createSCSMatrix(
            scs_int rows, scs_int cols, scs_int nnz, scs_float* x, scs_int* i, scs_int* p);

        /*!
         * \brief Create an SCS vector from an Eigen vector.
         * \param src The Eigen vector to convert.
         * \return Pointer to the created SCS vector.
         */
        scs_float* createCVector(const Eigen::VectorXd& src);
    };

}

#endif // scssolver_hpp
