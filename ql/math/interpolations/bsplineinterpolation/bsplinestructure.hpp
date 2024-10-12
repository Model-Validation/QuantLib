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

/*! \file bsplinestructure.hpp
    \brief Structure to manage composite B-spline segments
*/


#ifndef composite_spline_structure_hpp
#define composite_spline_structure_hpp

#include "bsplineevaluator.hpp"
#include "splineconstraints.hpp"
#include "splinesegment.hpp"
#include "splinesegment.hpp"
#include <ql/shared_ptr.hpp>
#include <ql/types.hpp>
#include <vector>

namespace QuantLib {

    /*!
     * \brief Class for managing a composite B-spline structure.
     */
    class BSplineStructure {
      public:
        /*!
         * @brief Constructor for BSplineStructure with specified spline segments and
         * constraints.
         * @param splineSegments Vector of shared pointers to B-spline structures.
         * @param splineConstraints Shared pointer to spline constraints.
         */
        BSplineStructure(
            const std::vector<ext::shared_ptr<BSplineSegment>>& splineSegments,
            const ext::shared_ptr<SplineConstraints>& splineConstraints);

        /*!
         * \brief Default constructor for BSplineStructure.
         */
        BSplineStructure() = default;

        /*!
         * \brief Get the range of the spline.
         * \return A pair representing the range of the spline.
         */
        std::pair<Real, Real> range() const;

        /*!
         * \brief Set the constraints for the spline.
         * \param splineConstraints Shared pointer to spline constraints.
         */
        void setConstraints(const ext::shared_ptr<SplineConstraints>& splineConstraints);

        /*!
         * \brief Transforms values according to the transforms applicable for the spline.
         * \param abscissae The abscissa values for the points
         * \param values The values to be transformed
         */
        std::vector<Real> transform(const std::vector<Real>& abscissae,
                                    const std::vector<Real>& values) const;

        /*!
         * \brief Add interpolation nodes to the spline.
         * \param interpolationNodes Vector of interpolation nodes.
         */
        void addInterpolationNodes(const std::vector<Real>& interpolationNodes) const;

        /*!
         * \brief Interpolate the spline at given nodes and values.
         * \param interpolationNodes Vector of interpolation nodes.
         * \param values Vector of values at the nodes.
         * \return A vector of interpolated values.
         */
        Eigen::VectorXd interpolate(const std::vector<Real>& interpolationNodes,
                                    const std::vector<Real>& values) const;

        /*!
         * \brief Solve the spline for given parameters.
         * \param parameters Vector of parameters.
         * \return A vector of solutions.
         */
        Eigen::VectorXd solve(const std::vector<Real>& parameters) const;

        /*!
         * \brief Get the solution of the spline.
         * \return A vector representing the solution.
         */
        Eigen::VectorXd getSolution() const;

        /*!
         * \brief Get the number of variables in the spline.
         * \return The number of variables.
         */
        Size getNumVariables() const { return splineConstraints_->getNumVariables(); }

        /*!
         * \brief Evaluate the spline structure at a given point.
         * \param x The point at which to evaluate.
         * \return Vector of evaluated values over spline basis functions.
         */
        Eigen::VectorXd evaluateAll(Real x) const;

        /*!
         * \brief Evaluate the spline structure at a given point.
         * \param coefficients Coefficients of the spline
         * \param x The point at which to evaluate.
         * \param side The side for evaluation.
         * \return Vector of evaluated values.
         */
        Real value(const Eigen::VectorXd& coefficients, Real x, BSplineSegment::Side side=BSplineSegment::Side::Right) const;

        /*!
         * \brief Get the spline segments.
         * \return Vector of shared pointers to the spline segments.
         */
        const std::vector<ext::shared_ptr<BSplineSegment>>& getSplineSegments() const {
            return splineSegments_;
        }

      private:
        std::vector<ext::shared_ptr<BSplineSegment>> splineSegments_;
        ext::shared_ptr<SplineConstraints> splineConstraints_;

        BSplineEvaluator spline_;
    };
}
#endif // composite_spline_structure_hpp
