/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2002, 2003, 2008, 2009 Ferdinando Ametrano
 Copyright (C) 2004, 2007, 2008 StatPro Italia srl
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

/*! \file ratetimeinterpolation.hpp
    \brief ratetime-linear and ratetime-cubic interpolation between discrete points
*/

#ifndef quantlib_ratetime_interpolation_hpp
#define quantlib_ratetime_interpolation_hpp

#include <ql/errors.hpp>
#include <ql/types.hpp>
#include <ql/math/interpolations/cubicinterpolation.hpp>
#include <ql/math/interpolations/linearinterpolation.hpp>
#include <ql/math/interpolations/mixedinterpolation.hpp>
#include <ql/math/interpolations/bsplineinterpolation/bsplineinterpolation.hpp>
#include <vector>
#include <utility>

namespace QuantLib {

    namespace detail {
        template <class I1, class I2, class I>
        class RateTimeInterpolationImpl;
        template <class I1, class I2, class IN1, class IN2>
        class RateTimeMixedInterpolationImpl;
    }

    //! %ratetime-linear interpolation between discrete points
    /*! \ingroup interpolations
        \warning See the Interpolation class for information about the
                 required lifetime of the underlying data.
    */
    class RateTimeLinearInterpolation : public Interpolation {
    public:
        /*! \pre the \f$ x \f$ values must be sorted. */
        template <class I1, class I2>
        RateTimeLinearInterpolation(const I1& xBegin, const I1& xEnd, const I2& yBegin) {
            impl_ = ext::shared_ptr<Interpolation::Impl>(
                new detail::RateTimeInterpolationImpl<I1, I2, Linear>(xBegin, xEnd, yBegin));
            impl_->update();
        }
    };

    //! ratetime-linear interpolation factory and traits
    /*! \ingroup interpolations */
    class RateTimeLinear {
    public:
        template <class I1, class I2>
        // ReSharper disable once CppMemberFunctionMayBeStatic
        Interpolation interpolate(const I1& xBegin, const I1& xEnd, const I2& yBegin) const {
            return static_cast<Interpolation>(RateTimeLinearInterpolation(xBegin, xEnd, yBegin));
        }

        static constexpr bool global = false;
        static constexpr Size requiredPoints = 2;
    };

    //! %ratetime-linear interpolation between discrete points
    /*! \ingroup interpolations
        \warning See the Interpolation class for information about the
                 required lifetime of the underlying data.
    */
    class RateTimeBSplineInterpolation : public Interpolation {
    public:
        /*! \pre the \f$ x \f$ values must be sorted. */
        template <class I1, class I2>
        RateTimeBSplineInterpolation(const I1& xBegin,
                                     const I1& xEnd,
                                     const I2& yBegin,
                                     ext::shared_ptr<BSplineStructure> splineStructure = {}) {
            impl_ = ext::shared_ptr<Interpolation::Impl>(
                new detail::RateTimeInterpolationImpl<I1, I2, BSplineModel>(
                    xBegin, xEnd, yBegin, BSplineModel(splineStructure)));
            impl_->update();
        }
    };

    class RateTimeBSpline {
    public:
        RateTimeBSpline(ext::shared_ptr<BSplineStructure>& splineStructure)
            :
            splineStructure_(std::move(splineStructure)) {
        };

        template <class I1, class I2>
        Interpolation interpolate(const I1& xBegin, const I1& xEnd, const I2& yBegin) const {
            return static_cast<Interpolation>(RateTimeBSplineInterpolation(xBegin, xEnd,
                yBegin, splineStructure_));
        }

        static constexpr bool global = false;
        static constexpr Size requiredPoints = 2;

        // Getter functions for the start and end points of splineSegment_
        double getStartPoint() const {
            return splineStructure_->range().first;
        }

        double getEndPoint() const { return splineStructure_->range().second; }

    private:
        ext::shared_ptr<BSplineStructure> splineStructure_;
    };

    //! %ratetime-cubic interpolation between discrete points
    /*! \ingroup interpolations */
    class RateTimeCubicInterpolation : public Interpolation {
    public:
        /*! \pre the \f$ x \f$ values must be sorted. */
        template <class I1, class I2>
        RateTimeCubicInterpolation(const I1& xBegin,
                                   const I1& xEnd,
                                   const I2& yBegin,
                                   CubicInterpolation::DerivativeApprox da,
                                   bool monotonic,
                                   CubicInterpolation::BoundaryCondition leftC,
                                   Real leftConditionValue,
                                   CubicInterpolation::BoundaryCondition rightC,
                                   Real rightConditionValue) {
            impl_ = ext::shared_ptr<Interpolation::Impl>(
                new detail::RateTimeInterpolationImpl<I1, I2, Cubic>(
                    xBegin, xEnd, yBegin,
                    Cubic(da, monotonic, leftC, leftConditionValue, rightC, rightConditionValue)));
            impl_->update();
        }
    };

    //! ratetime-cubic interpolation factory and traits
    /*! \ingroup interpolations */
    class RateTimeCubic {
    public:
        RateTimeCubic(CubicInterpolation::DerivativeApprox da,
                      bool monotonic = true,
                      CubicInterpolation::BoundaryCondition leftCondition =
                          CubicInterpolation::SecondDerivative,
                      Real leftConditionValue = 0.0,
                      CubicInterpolation::BoundaryCondition rightCondition =
                          CubicInterpolation::SecondDerivative,
                      Real rightConditionValue = 0.0)
            : da_(da), monotonic_(monotonic), leftType_(leftCondition), rightType_(rightCondition),
              leftValue_(leftConditionValue), rightValue_(rightConditionValue) {
        }

        template <class I1, class I2>
        Interpolation interpolate(const I1& xBegin, const I1& xEnd, const I2& yBegin) const {
            return static_cast<Interpolation>(RateTimeCubicInterpolation(
                xBegin, xEnd, yBegin, da_, monotonic_, leftType_,
                leftValue_, rightType_, rightValue_));
        }

        static constexpr bool global = true;
        static constexpr Size requiredPoints = 2;

    private:
        CubicInterpolation::DerivativeApprox da_;
        bool monotonic_;
        CubicInterpolation::BoundaryCondition leftType_, rightType_;
        Real leftValue_, rightValue_;
    };

    // convenience classes

    class DefaultRateTimeCubic : public RateTimeCubic {
    public:
        DefaultRateTimeCubic() : RateTimeCubic(CubicInterpolation::Kruger) {
        }
    };

    class MonotonicRateTimeCubic : public RateTimeCubic {
    public:
        MonotonicRateTimeCubic()
            : RateTimeCubic(CubicInterpolation::Spline,
                            true,
                            CubicInterpolation::SecondDerivative,
                            0.0,
                            CubicInterpolation::SecondDerivative,
                            0.0) {
        }
    };

    class KrugerRateTime : public RateTimeCubic {
    public:
        KrugerRateTime()
            : RateTimeCubic(CubicInterpolation::Kruger,
                            false,
                            CubicInterpolation::SecondDerivative,
                            0.0,
                            CubicInterpolation::SecondDerivative,
                            0.0) {
        }
    };


    class RateTimeCubicNaturalSpline : public RateTimeCubicInterpolation {
    public:
        /*! \pre the \f$ x \f$ values must be sorted. */
        template <class I1, class I2>
        RateTimeCubicNaturalSpline(const I1& xBegin, const I1& xEnd, const I2& yBegin)
            : RateTimeCubicInterpolation(xBegin,
                                         xEnd,
                                         yBegin,
                                         CubicInterpolation::Spline,
                                         false,
                                         CubicInterpolation::SecondDerivative,
                                         0.0,
                                         CubicInterpolation::SecondDerivative,
                                         0.0) {
        }
    };

    class MonotonicRateTimeCubicNaturalSpline : public RateTimeCubicInterpolation {
    public:
        /*! \pre the \f$ x \f$ values must be sorted. */
        template <class I1, class I2>
        MonotonicRateTimeCubicNaturalSpline(const I1& xBegin, const I1& xEnd, const I2& yBegin)
            : RateTimeCubicInterpolation(xBegin,
                                         xEnd,
                                         yBegin,
                                         CubicInterpolation::Spline,
                                         true,
                                         CubicInterpolation::SecondDerivative,
                                         0.0,
                                         CubicInterpolation::SecondDerivative,
                                         0.0) {
        }
    };

    class KrugerRateTimeCubic : public RateTimeCubicInterpolation {
    public:
        /*! \pre the \f$ x \f$ values must be sorted. */
        template <class I1, class I2>
        KrugerRateTimeCubic(const I1& xBegin, const I1& xEnd, const I2& yBegin)
            : RateTimeCubicInterpolation(xBegin,
                                         xEnd,
                                         yBegin,
                                         CubicInterpolation::Kruger,
                                         false,
                                         CubicInterpolation::SecondDerivative,
                                         0.0,
                                         CubicInterpolation::SecondDerivative,
                                         0.0) {
        }
    };

    class HarmonicRateTimeCubic : public RateTimeCubicInterpolation {
    public:
        /*! \pre the \f$ x \f$ values must be sorted. */
        template <class I1, class I2>
        HarmonicRateTimeCubic(const I1& xBegin, const I1& xEnd, const I2& yBegin)
            : RateTimeCubicInterpolation(xBegin,
                                         xEnd,
                                         yBegin,
                                         CubicInterpolation::Harmonic,
                                         false,
                                         CubicInterpolation::SecondDerivative,
                                         0.0,
                                         CubicInterpolation::SecondDerivative,
                                         0.0) {
        }
    };

    class FritschButlandRateTimeCubic : public RateTimeCubicInterpolation {
    public:
        /*! \pre the \f$ x \f$ values must be sorted. */
        template <class I1, class I2>
        FritschButlandRateTimeCubic(const I1& xBegin, const I1& xEnd, const I2& yBegin)
            : RateTimeCubicInterpolation(xBegin,
                                         xEnd,
                                         yBegin,
                                         CubicInterpolation::FritschButland,
                                         false,
                                         CubicInterpolation::SecondDerivative,
                                         0.0,
                                         CubicInterpolation::SecondDerivative,
                                         0.0) {
        }
    };

    class RateTimeParabolic : public RateTimeCubicInterpolation {
    public:
        /*! \pre the \f$ x \f$ values must be sorted. */
        template <class I1, class I2>
        RateTimeParabolic(const I1& xBegin, const I1& xEnd, const I2& yBegin)
            : RateTimeCubicInterpolation(xBegin,
                                         xEnd,
                                         yBegin,
                                         CubicInterpolation::Parabolic,
                                         false,
                                         CubicInterpolation::SecondDerivative,
                                         0.0,
                                         CubicInterpolation::SecondDerivative,
                                         0.0) {
        }
    };

    class MonotonicRateTimeParabolic : public RateTimeCubicInterpolation {
    public:
        /*! \pre the \f$ x \f$ values must be sorted. */
        template <class I1, class I2>
        MonotonicRateTimeParabolic(const I1& xBegin, const I1& xEnd, const I2& yBegin)
            : RateTimeCubicInterpolation(xBegin,
                                         xEnd,
                                         yBegin,
                                         CubicInterpolation::Parabolic,
                                         true,
                                         CubicInterpolation::SecondDerivative,
                                         0.0,
                                         CubicInterpolation::SecondDerivative,
                                         0.0) {
        }
    };

    //! %ratetime-mixedlinearcubic interpolation between discrete points
    /*! \ingroup interpolations */
    class RateTimeMixedLinearCubicInterpolation : public Interpolation {
    public:
        /*! \pre the \f$ x \f$ values must be sorted. */
        template <class I1, class I2>
        RateTimeMixedLinearCubicInterpolation(const I1& xBegin,
                                              const I1& xEnd,
                                              const I2& yBegin,
                                              const Size n,
                                              MixedInterpolation::Behavior behavior,
                                              CubicInterpolation::DerivativeApprox da,
                                              bool monotonic,
                                              CubicInterpolation::BoundaryCondition leftC,
                                              Real leftConditionValue,
                                              CubicInterpolation::BoundaryCondition rightC,
                                              Real rightConditionValue) {
            impl_ = ext::shared_ptr<Interpolation::Impl>(
                new detail::RateTimeInterpolationImpl<I1, I2, MixedLinearCubic>(
                    xBegin, xEnd, yBegin,
                    MixedLinearCubic(n, behavior, da, monotonic, leftC, leftConditionValue, rightC,
                                     rightConditionValue)));
            impl_->update();
        }
    };

    //! ratetime-cubic interpolation factory and traits
    /*! \ingroup interpolations */
    class RateTimeMixedLinearCubic {
    public:
        RateTimeMixedLinearCubic(const Size n,
                                 MixedInterpolation::Behavior behavior,
                                 CubicInterpolation::DerivativeApprox da,
                                 bool monotonic = true,
                                 CubicInterpolation::BoundaryCondition leftCondition =
                                     CubicInterpolation::SecondDerivative,
                                 Real leftConditionValue = 0.0,
                                 CubicInterpolation::BoundaryCondition rightCondition =
                                     CubicInterpolation::SecondDerivative,
                                 Real rightConditionValue = 0.0)
            : n_(n), behavior_(behavior), da_(da), monotonic_(monotonic), leftType_(leftCondition),
              rightType_(rightCondition), leftValue_(leftConditionValue),
              rightValue_(rightConditionValue) {
        }

        template <class I1, class I2>
        Interpolation interpolate(const I1& xBegin, const I1& xEnd, const I2& yBegin) const {
            return static_cast<Interpolation>(RateTimeMixedLinearCubicInterpolation(
                xBegin, xEnd, yBegin, n_, behavior_, da_,
                monotonic_, leftType_, leftValue_,
                rightType_, rightValue_));
        }

        static constexpr bool global = true;
        static constexpr Size requiredPoints = 3;

    private:
        Size n_;
        MixedInterpolation::Behavior behavior_;
        CubicInterpolation::DerivativeApprox da_;
        bool monotonic_;
        CubicInterpolation::BoundaryCondition leftType_, rightType_;
        Real leftValue_, rightValue_;
    };

    // convenience classes

    class DefaultRateTimeMixedLinearCubic : public RateTimeMixedLinearCubic {
    public:
        explicit DefaultRateTimeMixedLinearCubic(
            const Size n,
            MixedInterpolation::Behavior behavior = MixedInterpolation::ShareRanges)
            : RateTimeMixedLinearCubic(n, behavior, CubicInterpolation::Kruger) {
        }
    };

    class MonotonicRateTimeMixedLinearCubic : public RateTimeMixedLinearCubic {
    public:
        explicit MonotonicRateTimeMixedLinearCubic(
            const Size n,
            MixedInterpolation::Behavior behavior = MixedInterpolation::ShareRanges)
            : RateTimeMixedLinearCubic(n,
                                       behavior,
                                       CubicInterpolation::Spline,
                                       true,
                                       CubicInterpolation::SecondDerivative,
                                       0.0,
                                       CubicInterpolation::SecondDerivative,
                                       0.0) {
        }
    };

    class KrugerRateTimeMixedLinearCubic : public RateTimeMixedLinearCubic {
    public:
        explicit KrugerRateTimeMixedLinearCubic(
            const Size n,
            MixedInterpolation::Behavior behavior = MixedInterpolation::ShareRanges)
            : RateTimeMixedLinearCubic(n,
                                       behavior,
                                       CubicInterpolation::Kruger,
                                       false,
                                       CubicInterpolation::SecondDerivative,
                                       0.0,
                                       CubicInterpolation::SecondDerivative,
                                       0.0) {
        }
    };


    class RateTimeMixedLinearCubicNaturalSpline : public RateTimeMixedLinearCubicInterpolation {
    public:
        /*! \pre the \f$ x \f$ values must be sorted. */
        template <class I1, class I2>
        RateTimeMixedLinearCubicNaturalSpline(
            const I1& xBegin,
            const I1& xEnd,
            const I2& yBegin,
            const Size n,
            MixedInterpolation::Behavior behavior = MixedInterpolation::ShareRanges)
            : RateTimeMixedLinearCubicInterpolation(xBegin,
                                                    xEnd,
                                                    yBegin,
                                                    n,
                                                    behavior,
                                                    CubicInterpolation::Spline,
                                                    false,
                                                    CubicInterpolation::SecondDerivative,
                                                    0.0,
                                                    CubicInterpolation::SecondDerivative,
                                                    0.0) {
        }
    };

    /*! \ingroup interpolations
        \warning See the Interpolation class for information about the
                 required lifetime of the underlying data.
    */
    class MixedRateTimeLinearCubicInterpolation : public Interpolation {
    public:
        /*! \pre the \f$ x \f$ values must be sorted. */
        template <class I1, class I2>
        MixedRateTimeLinearCubicInterpolation(const I1& xBegin,
                                              const I1& xEnd,
                                              const I2& yBegin,
                                              Size n,
                                              MixedInterpolation::Behavior behavior,
                                              CubicInterpolation::DerivativeApprox da,
                                              bool monotonic,
                                              CubicInterpolation::BoundaryCondition leftC,
                                              Real leftConditionValue,
                                              CubicInterpolation::BoundaryCondition rightC,
                                              Real rightConditionValue) {
            impl_ = ext::shared_ptr<Interpolation::Impl>(
                new detail::MixedInterpolationImpl<I1, I2, RateTimeLinear, Cubic>(
                    xBegin, xEnd, yBegin, n, behavior, RateTimeLinear(),
                    Cubic(da, monotonic, leftC, leftConditionValue, rightC, rightConditionValue)));
            impl_->update();
        }
    };

    //! mixed linear/cubic interpolation factory and traits
    /*! \ingroup interpolations */
    class MixedRateTimeLinearCubic {
    public:
        MixedRateTimeLinearCubic(Size n,
                                 MixedInterpolation::Behavior behavior,
                                 CubicInterpolation::DerivativeApprox da,
                                 bool monotonic = true,
                                 CubicInterpolation::BoundaryCondition leftCondition =
                                     CubicInterpolation::SecondDerivative,
                                 Real leftConditionValue = 0.0,
                                 CubicInterpolation::BoundaryCondition rightCondition =
                                     CubicInterpolation::SecondDerivative,
                                 Real rightConditionValue = 0.0)
            : n_(n), behavior_(behavior), da_(da), monotonic_(monotonic), leftType_(leftCondition),
              rightType_(rightCondition), leftValue_(leftConditionValue),
              rightValue_(rightConditionValue) {
        }

        template <class I1, class I2>
        MixedRateTimeLinearCubicInterpolation
        interpolate(const I1& xBegin, const I1& xEnd, const I2& yBegin) const {
            return MixedRateTimeLinearCubicInterpolation(xBegin, xEnd, yBegin, n_, behavior_, da_,
                                                         monotonic_, leftType_, leftValue_,
                                                         rightType_, rightValue_);
        }

        // fix below
        // ReSharper disable once CppVariableCanBeMadeConstexpr
        static const bool global = true;
        // ReSharper disable once CppVariableCanBeMadeConstexpr
        static const Size requiredPoints = 3;

    private:
        Size n_;
        MixedInterpolation::Behavior behavior_;
        CubicInterpolation::DerivativeApprox da_;
        bool monotonic_;
        CubicInterpolation::BoundaryCondition leftType_, rightType_;
        Real leftValue_, rightValue_;
    };

    class MixedRateTimeLinearParabolic : public MixedRateTimeLinearCubic {
    public:
        MixedRateTimeLinearParabolic(Size n)
            : MixedRateTimeLinearCubic(n,
                                       MixedInterpolation::SplitRanges,
                                       CubicInterpolation::Parabolic,
                                       false) {
        }
    };

    ///*! \ingroup interpolations
    //    \warning See the Interpolation class for information about the
    //             required lifetime of the underlying data.
    //*/
    //class MixedRateTimeBSplineBSplineInterpolation : public Interpolation {
    //  public:
    //    /*! \pre the \f$ x \f$ values must be sorted. */
    //    template <class I1, class I2>
    //    MixedRateTimeBSplineBSplineInterpolation(const I1& xBegin,
    //                                          const I1& xEnd,
    //                                          const I2& yBegin,
    //                                          ext::shared_ptr<BSplineSegment>& splineSegment1,
    //                                          ext::shared_ptr<BSplineSegment>& splineSegment2
    //    ) {
    //        impl_ = ext::shared_ptr<Interpolation::Impl>(
    //            new detail::MixedBSplineInterpolationImpl<I1, I2, RateTimeBSpline, BSplineModel>(
    //                xBegin, xEnd, yBegin, RateTimeBSpline(splineSegment1), BSplineModel(splineSegment2)
    //            )
    //        );
    //        impl_->update();
    //    }
    //};

    ////! mixed rate time b-spline and b-spline interpolation factory and traits
    ///*! \ingroup interpolations */

    //class MixedRateTimeBSplineBSpline {
    //  public:
    //    MixedRateTimeBSplineBSpline(ext::shared_ptr<BSplineSegment>& splineSegment1,
    //                                ext::shared_ptr<BSplineSegment>& splineSegment2)
    //    : splineSegment1_(splineSegment1), splineSegment2_(splineSegment2) {}

    //    template <class I1, class I2>
    //    Interpolation interpolate(const I1& xBegin, const I1& xEnd, const I2& yBegin) const {
    //        return MixedRateTimeBSplineBSplineInterpolation(xBegin, xEnd, yBegin, splineSegment1_, splineSegment2_);
    //    }
    //    // fix below
    //    static const bool global = true;
    //    static const Size requiredPoints = 2;

    //  private:
    //    ext::shared_ptr<BSplineSegment>& splineSegment1_;
    //    ext::shared_ptr<BSplineSegment>& splineSegment2_;
    //};

    // convenience classes
    namespace detail {
        template <class I1, class I2, class Interpolator>
        class RateTimeInterpolationImpl : public Interpolation::templateImpl<I1, I2> {
        public:
            RateTimeInterpolationImpl(const I1& xBegin,
                                      const I1& xEnd,
                                      const I2& yBegin,
                                      const Interpolator& factory = Interpolator())
                : Interpolation::templateImpl<I1, I2>(
                      xBegin, xEnd, yBegin, Interpolator::requiredPoints),
                  ratetimeY_(std::distance(xBegin, xEnd)) {
                interpolation_ =
                    factory.interpolate(this->xBegin_, this->xEnd_, ratetimeY_.begin());
            }

            void update() override {
                for (Size i = 0; i < ratetimeY_.size(); ++i) {
                    ratetimeY_[i] = this->yBegin_[i] * this->xBegin_[i];
                }
                interpolation_.update();
            }

            Real value(Real x) const override {
                if (x == 0.0) {
                    return interpolation_.derivative(0.0, true);
                } else {
                    return interpolation_(x, true) / x;
                }
            }

            Real derivative(Real x) const override {
                if (x == 0.0) {
                    return interpolation_.secondDerivative(0.0, true) / 2.0;
                } else {
                    return (interpolation_.derivative(x, true) - value(x)) / x;
                }
            }

            Real primitive(Real) const override {
                QL_FAIL("RateTimeInterpolation primitive not implemented");
            }

            Real secondDerivative(Real x) const override {
                if (x == 0.0) {
                    Real dt = 1.0e-6;
                    return (derivative(dt) - derivative(0.0)) / dt;
                } else {
                    return (interpolation_.derivative(x, true) - 2.0 * derivative(x)) / x;
                }
            }

        private:
            std::vector<Real> ratetimeY_;
            Interpolation interpolation_;
        };

    }
}

#endif