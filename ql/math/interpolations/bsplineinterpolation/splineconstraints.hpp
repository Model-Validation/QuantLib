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

#ifndef spline_constraints_hpp
#define spline_constraints_hpp

#include <stack>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <ql/types.hpp>
#include <scs.h>
#include <scs_types.h>
#include <tuple>


namespace SCS {
    class SCSSolver {
      public:
        SCSSolver();
        SCSSolver(const Eigen::SparseMatrix<double>& P,
                  const Eigen::SparseMatrix<double>& A,
                  const Eigen::VectorXd& b,
                  const Eigen::VectorXd& c,
                  QuantLib::Integer nEqualities,
                  QuantLib::Integer nInequalities);

        ~SCSSolver();

        // Method to convert Eigen structures to SCS format
        void convertEigenToCSC(const Eigen::SparseMatrix<double>& mat,
                               std::vector<scs_float>& scs_x,
                               std::vector<scs_int>& scs_i,
                               std::vector<scs_int>& scs_p);

        void initializeSCS();

        QuantLib::Integer solve();
        QuantLib::Integer solve(int warm_start);

        QuantLib::Integer update(Eigen::VectorXd& b, Eigen::VectorXd& c);

        Eigen::VectorXd solution_x();
        Eigen::VectorXd solution_y();
        Eigen::VectorXd solution_s();

      private:
        // SCS data structures
        std::vector<scs_float> scs_P_x_;
        std::vector<scs_int> scs_P_i_;
        std::vector<scs_int> scs_P_p_;
        std::vector<scs_float> scs_A_x_;
        std::vector<scs_int> scs_A_i_;
        std::vector<scs_int> scs_A_p_;
        Eigen::VectorXd scs_b_;
        Eigen::VectorXd scs_c_;

        scs_int nVariables_, nConstraints_, nEqualities_, nInequalities_;

        bool is_updated_b_but_data_not_reset_, is_updated_c_but_data_not_reset_;

        ScsData scs_data_;
        ScsCone scs_cone_;
        ScsSettings scs_settings_;
        ScsInfo scs_info_;
        ScsWork* scs_work_;
        ScsSolution scs_sol_;

        ScsMatrix* createSCSMatrix(
            scs_int rows, scs_int cols, scs_int nnz, scs_float* x, scs_int* i, scs_int* p);
        scs_float* SCSSolver::createCVector(const Eigen::VectorXd& src);
    };
}

namespace PermutationIterators {
    // Custom iterator that applies the permutation
    class PermutationIteratorTriplet {
      public:
        using Triplet = Eigen::Triplet<double>;
        using Iterator = std::vector<Triplet>::const_iterator;

        PermutationIteratorTriplet(Iterator it, const std::vector<int>& permutation)
        : it_(it), permutation_(permutation) {}

        PermutationIteratorTriplet& operator++() {
            ++it_;
            return *this;
        }

        bool operator!=(const PermutationIteratorTriplet& other) const { return it_ != other.it_; }

        Triplet operator*() const {
            return Triplet(permutation_[it_->row()], it_->col(), it_->value());
        }

      private:
        Iterator it_;
        const std::vector<int>& permutation_;
    };


    class PermutationIteratorVector {
      public:
        using Iterator = std::vector<double>::const_iterator;

        PermutationIteratorVector(Iterator it, const std::vector<int>& permutation)
        : it_(it), permutation_(permutation) {}

        PermutationIteratorVector& operator++() {
            ++it_;
            return *this;
        }

        bool operator!=(const PermutationIteratorVector& other) const { return it_ != other.it_; }

        //Triplet operator*() const {
        //    return T(permutation_[it_->row()], it_->col(), it_->value());
        //}

      private:
        Iterator it_;
        const std::vector<int>& permutation_;
    };
}

namespace QuantLib {

    class SplineConstraints {
      public:
        enum class ConstraintType { Equal, LessEqual };

        SplineConstraints(Size nVariables = 0,
                          const std::vector<std::vector<double>>& P = {},
                          const std::vector<std::vector<double>>& A = {},
                          const std::vector<double>& b = {},
                          const std::vector<double>& c = {},
                          const std::vector<ConstraintType>& constraintTypes = {}
                          //Size nParameters = 0,
                          //const std::vector<double>& parameters = {},
                          //const std::vector<std::vector<double>>& B = {},
                          //const std::vector<std::vector<double>>& C = {}
        );

        ~SplineConstraints() = default;

        SplineConstraints(const SplineConstraints& other);

        void updateOrdering();


        void addLinearConstraint(const Eigen::VectorXd& constraint,
                                 double rhs,
                                 ConstraintType constraintType = ConstraintType::Equal);

        //void addObjectiveFunction(const std::vector<std::vector<double>>& P,
        //                          const std::vector<double>& c = {});
        int solve();
        Eigen::VectorXd getSolution();
        Size getNVariables() const;
        Size getNConstraints() const;
        Size getNParameters() const;
        void update_b(std::vector<double> parameters);
        void update_c(std::vector<double> parameters);

        // Public methods
        //void setParameterMatrixB(Size nParameters = 0,
        //                         const std::vector<std::vector<double>>& B_parameterMatrix = {});
        void addParameters(Size nNewParameters,
                           Eigen::SparseMatrix<double>& B_new,
                           Eigen::SparseMatrix<double>& C_new);

        void push();
        void pop();

      private:
        Size nVariables_, nConstraints_, nEqualities_, nInequalities_, nParameters_;
        Eigen::SparseMatrix<double> P_; // Primary
        std::vector<Eigen::Triplet<double>> A_triplets_; // Primary
        Eigen::SparseMatrix<double> A_; // Secondary, needs to be updated
        Eigen::VectorXd c_; // Primary
        std::vector<double> c_list_; // Secondary, needs to be updated
        std::vector<double> b_list_; // Primary
        Eigen::VectorXd b_;          // Secondary, needs to be updated
        std::vector<ConstraintType> constraintTypes_; // Primary
        std::vector<double> parameters_list_; // Primary
        Eigen::VectorXd parameters_; // Secondary
        std::vector<Eigen::Triplet<double>> B_triplets_; // Primary
        Eigen::SparseMatrix<double> B_; // Secondary, needs to be updated
        std::vector<Eigen::Triplet<double>> C_triplets_; // Primary
        Eigen::SparseMatrix<double> C_; // Secondary, needs to be updated
        std::stack<std::tuple<Size, Size, Size, Size, Size>> constraintStack_;

       
        SCS::SCSSolver* scsData_ = nullptr; 

        // State information
        bool scsDataIsUpToDate_ = false; // True if the SCS data structures are up to date
        bool isOrdered_ = false; // True if the constrained have been ordered according to constraint type
        bool isSolved_ = false; // True if the system has been solved at least once, so a solution is available
        bool hasParameters_ = false; // True if the constraints have parameters (usually interpolation nodes)
        std::vector<int> permutation_; // Permutations to order constraints
        int warmStart_ = 0;

        // methods
        void reorderByConstraints();
        void updateParameterValues(const Eigen::VectorXd& parameters);
        void updateParameterValues(const std::vector<double> parameters);
        void updateScsData();
    };
}

#endif // spline_constraints_hpp