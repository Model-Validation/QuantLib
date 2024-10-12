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

#include <vector>
#include <Eigen/Sparse>
#include <ql/types.hpp>

using namespace QuantLib;

namespace EigenUtilities {
    std::vector<Eigen::Triplet<double>>
    convertToTriplets(const std::vector<std::vector<double>>& matrix) {
        std::vector<Eigen::Triplet<double>> triplets;

        for (Size i = 0; i < matrix.size(); ++i) {
            for (Size j = 0; j < matrix[i].size(); ++j) {
                if (matrix[i][j] != 0.0) {
                    triplets.emplace_back(i, j, matrix[i][j]);
                }
            }
        }
        return triplets;
    }

    std::vector<Eigen::Triplet<double>>
    convertToTriplets(const Eigen::SparseMatrix<double>& matrix) {
        std::vector<Eigen::Triplet<double>> triplets(matrix.nonZeros()); // Preallocate space

        Size index = 0; // Keep track of the current index

        for (int k = 0; k < matrix.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, k); it; ++it) {
                triplets[index++] = Eigen::Triplet<double>(static_cast<Size>(it.row()),
                                                           static_cast<Size>(it.col()), it.value());
            }
        }
        return triplets;
    }

    Eigen::SparseMatrix<double>
    convertToEigenSparseMatrix(const std::vector<std::vector<double>>& vecMatrix) {
        Size rows = vecMatrix.size();
        Size cols = vecMatrix[0].size();
        Eigen::SparseMatrix<double> sparseMatrix(rows, cols);

        for (Size i = 0; i < rows; ++i) {
            for (Size j = 0; j < cols; ++j) {
                if (vecMatrix[i][j] != 0.0) {
                    sparseMatrix.insert(i, j) = vecMatrix[i][j];
                }
            }
        }

        return sparseMatrix;
    }

    Eigen::VectorXd convertToEigenVector(const std::vector<double>& vec) {
        Eigen::VectorXd eigenVector(vec.size());
        for (Size i = 0; i < vec.size(); ++i) {
            eigenVector[i] = vec[i];
        }
        return eigenVector;
    }
}