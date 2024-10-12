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

#ifndef eigen_utilities_hpp
#define eigen_utilities_hpp

#include <Eigen/Sparse>
#include <vector>

namespace EigenUtilities {

    /*!
     * \brief Convert a 2D vector to a list of Eigen triplets.
     * \param matrix 2D vector representing the matrix.
     * \return Vector of Eigen triplets.
     */
    std::vector<Eigen::Triplet<double>>
    convertToTriplets(const std::vector<std::vector<double>>& matrix);

    /*!
     * \brief Convert an Eigen sparse matrix to a list of Eigen triplets.
     * \param matrix Eigen sparse matrix.
     * \return Vector of Eigen triplets.
     */
    std::vector<Eigen::Triplet<double>>
    convertToTriplets(const Eigen::SparseMatrix<double>& matrix);

    /*!
     * \brief Convert a 2D vector to an Eigen sparse matrix.
     * \param vecMatrix 2D vector representing the matrix.
     * \return Eigen sparse matrix.
     */
    Eigen::SparseMatrix<double>
    convertToEigenSparseMatrix(const std::vector<std::vector<double>>& vecMatrix);

    /*!
     * \brief Convert a 1D vector to an Eigen vector.
     * \param vec 1D vector.
     * \return Eigen vector.
     */
    Eigen::VectorXd convertToEigenVector(const std::vector<double>& vec);
}

#endif // eigen_utilities_hpp
