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

/*! \file splineconstraints.hpp
    \brief Spline constraints and SCS solver integration
*/

#ifndef spline_constraints_hpp
#define spline_constraints_hpp

#include "scssolver.hpp"
#include <ql/types.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <stack>
#include <tuple>
#include <vector>

namespace PermutationIterators {

    /*!
     * \brief Custom iterator that applies a permutation to triplets.
     */
    class PermutationIteratorTriplet {
      public:
        using Triplet = Eigen::Triplet<double>;
        using Iterator = std::vector<Triplet>::const_iterator;

        /*!
         * \brief Constructor for PermutationIteratorTriplet.
         * \param it Iterator to the triplets.
         * \param permutation Vector of permutation indices.
         */
        PermutationIteratorTriplet(Iterator it, std::vector<int> permutation)
        : it_(it), permutation_(std::move(permutation)) {}

        /*!
         * \brief Increment the iterator.
         * \return Reference to the incremented iterator.
         */
        PermutationIteratorTriplet& operator++() {
            ++it_;
            return *this;
        }

        /*!
         * \brief Compare two iterators for inequality.
         * \param other The other iterator to compare.
         * \return True if the iterators are not equal, false otherwise.
         */
        bool operator!=(const PermutationIteratorTriplet& other) const { return it_ != other.it_; }

        /*!
         * \brief Dereference the iterator to get the permuted triplet.
         * \return The permuted triplet.
         */
        Triplet operator*() const { return {permutation_[it_->row()], it_->col(), it_->value()}; }

      private:
        Iterator it_;
        std::vector<int> permutation_;
    };

    /*!
     * \brief Custom iterator that applies a permutation to vectors.
     */
    class PermutationIteratorVector {
      public:
        using Iterator = std::vector<double>::const_iterator;

        /*!
         * \brief Constructor for PermutationIteratorVector.
         * \param it Iterator to the vector.
         * \param permutation Vector of permutation indices.
         */
        PermutationIteratorVector(Iterator it, std::vector<int> permutation)
        : it_(it), permutation_(std::move(permutation)) {}

        /*!
         * \brief Increment the iterator.
         * \return Reference to the incremented iterator.
         */
        PermutationIteratorVector& operator++() {
            ++it_;
            return *this;
        }

        /*!
         * \brief Compare two iterators for inequality.
         * \param other The other iterator to compare.
         * \return True if the iterators are not equal, false otherwise.
         */
        bool operator!=(const PermutationIteratorVector& other) const { return it_ != other.it_; }

      private:
        Iterator it_;
        std::vector<int> permutation_;
    };
}

namespace QuantLib {

    /*!
     * \brief Class for managing spline constraints and objective forms
     */
    class SplineConstraints {
      public:
        /*!
         * \brief Enum representing the type of constraint.
         */
        enum class ConstraintType { Equal, LessEqual }; // NOLINT(performance-enum-size)

        /*!
         * \brief Default destructor for SplineConstraints.
         */
        SplineConstraints() = default;

        /*!
         * \brief Constructor for SplineConstraints with specified parameters.
         * \param numVariables Number of variables.
         * \param P_quadForm Quadratic form matrix.
         * \param A_constraints Constraint matrix.
         * \param b_rhs Right-hand side vector.
         * \param c_linearForm Linear form vector.
         * \param constraintTypes Vector of constraint types.
         * \param fitData Whether to fit data or set as constraints
         * \param epsAbsolute Termination criterion for absolute error, see e.g.
         * https://www.cvxgrp.org/scs/algorithm/index.html#termination
         * \param epsRelative Termination criterion for relative error, see e.g.
         * https://www.cvxgrp.org/scs/algorithm/index.html#termination
         * \param epsInfeasible Termination criterion for infeasibility i.e. primal or dual
         * unbounded, see e.g. https://www.cvxgrp.org/scs/algorithm/index.html#termination
         */
        SplineConstraints(Size numVariables,
                          const std::vector<std::vector<double>>& P_quadForm = {},
                          const std::vector<std::vector<double>>& A_constraints = {},
                          const std::vector<double>& b_rhs = {},
                          const std::vector<double>& c_linearForm = {},
                          const std::vector<ConstraintType>& constraintTypes = {},
                          bool fitData = false,
                          double epsAbsolute = 1e-12,
                          double epsRelative = 1e-12,
                          double epsInfeasible = 1e-13);

        /*!
         * \brief Constructor for SplineConstraints with specified Eigen matrices and vectors.
         * \param numVariables Number of variables.
         * \param P_quadForm Quadratic form matrix.
         * \param A_triplets Triplets for the constraint matrix.
         * \param b_rhs Right-hand side vector.
         * \param c_linearForm Linear form vector.
         * \param constraintTypes Vector of constraint types.
         * \param fitData Whether to fit data or set as constraints
         * \param epsAbsolute Termination criterion for absolute error, see e.g.
         * https://www.cvxgrp.org/scs/algorithm/index.html#termination
         * \param epsRelative Termination criterion for relative error, see e.g.
         * https://www.cvxgrp.org/scs/algorithm/index.html#termination
         * \param epsInfeasible Termination criterion for infeasibility i.e. primal or dual
         * unbounded, see e.g. https://www.cvxgrp.org/scs/algorithm/index.html#termination
         */
        SplineConstraints(Size numVariables,
                          const Eigen::SparseMatrix<double>& P_quadForm,
                          const std::vector<Eigen::Triplet<double>>& A_triplets = {},
                          const std::vector<double>& b_rhs = {},
                          const std::vector<double>& c_linearForm = {},
                          const std::vector<ConstraintType>& constraintTypes = {},
                          bool fitData = false,
                          double epsAbsolute = 1e-12,
                          double epsRelative = 1e-12,
                          double epsInfeasible = 1e-13);

        ///*!
        // * \brief Default destructor for SplineConstraints.
        // */
        //~SplineConstraints();

        ///*!
        // * \brief Copy constructor for SplineConstraints.
        // * \param other The other SplineConstraints to copy.
        // */
        // SplineConstraints(const SplineConstraints& other);

        /*!
         * \brief Update the ordering of constraints.
         */
        void updateOrdering();

        /*!
         * \brief Add a linear constraint to the system.
         * \param constraint The constraint vector.
         * \param rhs The right-hand side value.
         * \param constraintType The type of the constraint (default is Equal).
         */
        void addLinearConstraint(const Eigen::VectorXd& constraint,
                                 Real rhs,
                                 ConstraintType constraintType = ConstraintType::Equal);

        /*!
         * \brief Solve the system of constraints.
         * \return The status of the solver.
         */
        int solve();

        /*!
         * \brief Get the solution vector.
         * \return The solution vector.
         */
        Eigen::VectorXd getSolution() const;

        /*!
         * \brief Get the number of variables.
         * \return The number of variables.
         */
        Size getNumVariables() const;
        // ReSharper disable once CppInconsistentNaming
        Integer get_num_variables() const { return static_cast<Integer>(getNumVariables()); }

        /*!
         * \brief Get the number of constraints.
         * \return The number of constraints.
         */
        Size getNConstraints() const;
        // ReSharper disable once CppInconsistentNaming
        Integer get_num_constraints() const { return static_cast<Integer>(getNConstraints()); }

        /*!
         * \brief Get the number of parameters.
         * \return The number of parameters.
         */
        Size getNParameters() const;
        // ReSharper disable once CppInconsistentNaming
        Integer get_num_parameters() const { return static_cast<Integer>(getNParameters()); }

        /*!
         * \brief Update the right-hand side vector.
         * \param parameters The new parameters.
         */
        void update_b(const std::vector<Real>& parameters);

        /*!
         * \brief Update the linear form vector.
         * \param parameters The new parameters.
         */
        void update_c(const std::vector<Real>& parameters);

        /*!
         * \brief Add new parameters to the system.
         * \param nNewParameters Number of new parameters.
         * \param B_new New B matrix.
         * \param C_new New C matrix.
         */
        void addParameters(Size nNewParameters,
                           const Eigen::SparseMatrix<Real>& B_new,
                           const Eigen::SparseMatrix<Real>& C_new);

        /*!
         * \brief Push the current state onto the stack.
         */
        void push();

        /*!
         * \brief Pop the state from the stack.
         */
        void pop();

        /*!
         * \brief Get a slice of matrix A, this is a hack to solve overdetermined problems.
         */
        Eigen::SparseMatrix<Real> getSliceOfA(Integer firstRow, Integer lastRow) const;

        // ReSharper disable once CppInconsistentNaming
        [[nodiscard]] std::vector<std::vector<Real>> get_p_matrix() const {
            // Initialize the output vector
            std::vector<std::vector<Real>> result(this->P_.rows(),
                                                  std::vector<Real>(this->P_.cols(), 0.0));

            // Convert the sparse matrix to a dense std::vector<std::vector<Real>>
            for (int k = 0; k < this->P_.outerSize(); ++k) {
                for (Eigen::SparseMatrix<Real>::InnerIterator it(this->P_, k); it; ++it) {
                    result[it.row()][it.col()] = it.value();
                }
            }

            return result;
        }

        // ReSharper disable once CppInconsistentNaming
        [[nodiscard]] std::vector<std::vector<Real>> get_a_matrix() const {
            // Initialize the output vector
            std::vector<std::vector<Real>> result(this->A_.rows(),
                                                  std::vector<Real>(this->A_.cols(), 0.0));

            // Convert the sparse matrix to a dense std::vector<std::vector<Real>>
            for (int k = 0; k < this->A_.outerSize(); ++k) {
                for (Eigen::SparseMatrix<Real>::InnerIterator it(this->A_, k); it; ++it) {
                    result[it.row()][it.col()] = it.value();
                }
            }

            return result;
        }

        // ReSharper disable once CppInconsistentNaming
        [[nodiscard]] std::vector<Real> get_b_vector() const { return this->b_list_; }

        // ReSharper disable once CppInconsistentNaming
        [[nodiscard]] std::vector<Real> get_c_vector() const { return this->c_list_; }

        // ReSharper disable once CppInconsistentNaming
        [[nodiscard]] std::vector<Real> get_parameters() const { return this->parameters_list_; }

        // ReSharper disable once CppInconsistentNaming
        [[nodiscard]] std::vector<ConstraintType> get_constraint_types() const {
            return this->constraintTypes_;
        }
        bool fitData_ = false; // TODO: Hack, should not be public like this

        void setP(Eigen::SparseMatrix<Real>& P) {
            P_ = P;
            scsDataIsUpToDate_ = false;
        }
        void setCVector(Eigen::VectorXd& c) {
            c_ = c;
            c_list_ = std::vector<Real>(c_.data(), c_.data() + c_.size());

            scsDataIsUpToDate_ = false;
        }

      private:
        Size numVariables_, numConstraints_, numEqualities_, numInequalities_, numParameters_;
        Eigen::SparseMatrix<Real> P_; /*!< Primary quadratic form matrix. */
        std::vector<Eigen::Triplet<Real>>
            A_triplets_;              /*!< Primary triplets for the constraint matrix. */
        Eigen::SparseMatrix<Real> A_; /*!< Secondary constraint matrix, needs to be updated. */
        Eigen::VectorXd c_;           /*!< Primary linear form vector. */
        std::vector<Real> c_list_;    /*!< Secondary linear form list, needs to be updated. */
        std::vector<Real> b_list_;    /*!< Primary right-hand side list. */
        Eigen::VectorXd b_;           /*!< Secondary right-hand side vector, needs to be updated. */
        std::vector<ConstraintType> constraintTypes_;  /*!< Primary vector of constraint types. */
        std::vector<Real> parameters_list_;            /*!< Primary list of parameters. */
        Eigen::VectorXd parameters_;                   /*!< Secondary vector of parameters. */
        std::vector<Eigen::Triplet<Real>> B_triplets_; /*!< Primary triplets for the B matrix. */
        Eigen::SparseMatrix<Real> B_; /*!< Secondary B matrix, needs to be updated. */
        std::vector<Eigen::Triplet<Real>> C_triplets_; /*!< Primary triplets for the C matrix. */
        Eigen::SparseMatrix<Real> C_; /*!< Secondary C matrix, needs to be updated. */
        std::stack<std::tuple<Size, Size, Size, Size, Size, Size, Size>>
            constraintStack_; /*!< Stack for storing constraint states. */

        SCS::SCSSolver* scsData_ = nullptr; /*!< Pointer to the SCS solver data. */

        bool scsDataIsUpToDate_ = false; /*!< True if the SCS data structures are up to date. */
        bool isOrdered_ =
            false; /*!< True if the constraints have been ordered according to constraint type. */
        bool isSolved_ = false; /*!< True if the system has been solved at least once, so a solution
                                   is available. */
        bool hasParameters_ =
            false; /*!< True if the constraints have parameters (usually interpolation nodes). */
        std::vector<int> permutation_; /*!< Permutations to order constraints. */
        int warmStart_ = 0;            /*!< Warm start parameter. */

        double epsAbsolute_ = 1e-12;
        double epsRelative_ = 1e-12;
        double epsInfeasible_ = 1e-13;

        /*!
         * \brief Reorder constraints by their types.
         */
        void reorderByConstraints();

        /*!
         * \brief Update parameter values from an Eigen vector.
         * \param parameters The new parameter values.
         */
        void updateParameterValues(const Eigen::VectorXd& parameters);

        /*!
         * \brief Update parameter values from a vector.
         * \param parameters The new parameter values.
         */
        void updateParameterValues(const std::vector<double>& parameters);

        /*!
         * \brief Update the SCS data structures.
         */
        void updateScsData();
    };
}

#endif // spline_constraints_hpp
