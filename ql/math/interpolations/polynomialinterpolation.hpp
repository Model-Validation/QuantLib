/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2023 Skandinaviska Enskilda Banken AB (publ)

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

/*! \file polynomialinterpolation.hpp
    \brief polynomial interpolation between discrete points
*/

#ifndef quantlib_polynomial_interpolation_hpp
#define quantlib_polynomial_interpolation_hpp

#include <ql/math/interpolation.hpp>
#include <ql/math/matrix.hpp>
#include <ql/math/polynomialmathfunction.hpp>
#include <ql/math/linearleastsquaresregression.hpp>
#include <ql/methods/finitedifferences/tridiagonaloperator.hpp>
#include <vector>

namespace QuantLib {
    namespace detail {
        template <class I1, class I2>
        class PolynomialInterpolationImpl;
    }

    //! %Polynomial interpolation between discrete points
    /*! \ingroup interpolations
        \warning See the Interpolation class for information about the
                 required lifetime of the underlying data.
    */
    class PolynomialInterpolation : public Interpolation {
      public:
        /*! \pre the \f$ x \f$ values must be sorted. */
        template <class I1, class I2>
        PolynomialInterpolation(const I1& xBegin, const I1& xEnd, const I2& yBegin) {
            impl_ = ext::shared_ptr<Interpolation::Impl>(
                new detail::PolynomialInterpolationImpl<I1, I2>(xBegin, xEnd, yBegin));
            impl_->update();
        }
    };

    //! %Polynomial-interpolation factory and traits
    /*! \ingroup interpolations */
    class Polynomial1D {
      public:
        template <class I1, class I2>
        Interpolation interpolate(const I1& xBegin, const I1& xEnd, const I2& yBegin) const {
            return PolynomialInterpolation(xBegin, xEnd, yBegin);
        }
        static const bool global = false;
        static const Size requiredPoints = 1;
    };

    namespace detail {

        template <class I1, class I2>
        class PolynomialInterpolationImpl : public Interpolation::templateImpl<I1, I2> {
          public:
            PolynomialInterpolationImpl(const I1& xBegin, const I1& xEnd, const I2& yBegin)
            : Interpolation::templateImpl<I1, I2>(xBegin, xEnd, yBegin, Polynomial1D::requiredPoints),
              polynomial_(std::vector<Real>{1}) {}
            void update() override {
                std::vector<ext::function<Real(Real)>> v{};
                for (int i = 0; i < xValues().size(); ++i) {
                    v.push_back([i](Real x) -> Real { return std::pow(x, i); });
                }
                auto fit = LinearRegression(xValues(), yValues(), v);
                auto coeffs = fit.coefficients();

                polynomial_ = PolynomialFunction(std::vector<Real>(coeffs.begin(), coeffs.end()));
            }
            Real value(Real x) const override { return polynomial_(x); }
            Real primitive(Real x) const override { return polynomial_.primitive(x); }
            Real derivative(Real x) const override { return polynomial_.derivative(x); }
            Real secondDerivative(Real x) const override {
                QL_FAIL(
                    "secondDerivative has not yet been implemented for polynomial interpolation");
            }

          private:
            PolynomialFunction polynomial_;
        };

    }

}
#endif