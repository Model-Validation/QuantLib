/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2025 AcadiaSoft Inc.

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

#include <ql/math/comparison.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/pricingengines/bacheliercalculator.hpp>

namespace QuantLib {

    class BachelierCalculator::Calculator : public AcyclicVisitor,
                                            public Visitor<Payoff>,
                                            public Visitor<PlainVanillaPayoff>,
                                            public Visitor<CashOrNothingPayoff>,
                                            public Visitor<AssetOrNothingPayoff>,
                                            public Visitor<GapPayoff> {
      private:
        BachelierCalculator& bachelier_;

      public:
        explicit Calculator(BachelierCalculator& bachelier) : bachelier_(bachelier) {}
        void visit(Payoff&) override;
        void visit(PlainVanillaPayoff&) override;
        void visit(CashOrNothingPayoff&) override;
        void visit(AssetOrNothingPayoff&) override;
        void visit(GapPayoff&) override;
    };


    BachelierCalculator::BachelierCalculator(const ext::shared_ptr<StrikedTypePayoff>& p,
                                             Real forward,
                                             Real stdDev,
                                             Real discount)
    : strike_(p->strike()), forward_(forward), stdDev_(stdDev), discount_(discount),
      variance_(stdDev * stdDev) {
        initialize(p);
    }

    BachelierCalculator::BachelierCalculator(
        Option::Type optionType, Real strike, Real forward, Real stdDev, Real discount)
    : strike_(strike), forward_(forward), stdDev_(stdDev), discount_(discount),
      variance_(stdDev * stdDev) {
        initialize(ext::shared_ptr<StrikedTypePayoff>(new PlainVanillaPayoff(optionType, strike)));
    }

    void BachelierCalculator::initialize(const ext::shared_ptr<StrikedTypePayoff>& p) {
        QL_REQUIRE(stdDev_ >= 0.0, "stdDev (" << stdDev_ << ") must be non-negative");
        QL_REQUIRE(discount_ > 0.0, "discount (" << discount_ << ") must be positive");

        if (stdDev_ >= QL_EPSILON) {
            d_ = (forward_ - strike_) / stdDev_;
            auto f = CumulativeNormalDistribution();
            cum_d_ = f(d_);
            n_d_ = f.derivative(d_);

        } else {
            if (close(forward_, strike_)) {
                d_ = 0.0;
                cum_d_ = 0.5;
                n_d_ = M_SQRT_2 * M_1_SQRTPI;
            } else if (forward_ > strike_) {
                d_ = QL_MAX_REAL;
                cum_d_ = 1.0;
                n_d_ = 0.0;
            } else {
                d_ = QL_MIN_REAL;
                cum_d_ = 0.0;
                n_d_ = 0.0;
            }
        }

        x_ = forward_ - strike_;
        DxDforward_ = 1.0;
        DxDstrike_ = -1.0;

        y_ = stdDev_;
        DyDforward_ = 0.0;
        DyDstrike_ = 0.0;


        // this part is always executed.
        // in case of plain-vanilla payoffs, it is also the only part
        // which is executed.
        switch (p->optionType()) {
            case Option::Call:
                alpha_ = cum_d_;                     //  N(d)
                DalphaDd_ = n_d_;                    //  n(d)
                D2alphaD2d_ = -d_ * n_d_;            // -d * n(d)
                beta_ = n_d_;                        // n(d)
                DbetaDd_ = -d_ * n_d_;               // -d * n(d)
                D2betaD2d_ = (d_ * d_ - 1.0) * n_d_; // (d^2 -1) * n(d)
                break;
            case Option::Put:
                alpha_ = -1.0 + cum_d_;              // -N(-d)
                DalphaDd_ = n_d_;                    //  n(d)
                D2alphaD2d_ = -d_ * n_d_;            // -d * n(d)
                beta_ = n_d_;                        //  n(d)
                DbetaDd_ = -d_ * n_d_;               // -d * n(d)
                D2betaD2d_ = (d_ * d_ - 1.0) * n_d_; // (d^2 -1) * n(d)
                break;
            default:
                QL_FAIL("invalid option type");
        }

        // now dispatch on type.

        Calculator calc(*this);
        p->accept(calc);
    }

    void BachelierCalculator::Calculator::visit(Payoff& p) {
        QL_FAIL("unsupported payoff type: " << p.name());
    }

    void BachelierCalculator::Calculator::visit(PlainVanillaPayoff&) {}

    void BachelierCalculator::Calculator::visit(CashOrNothingPayoff& payoff) {
        bachelier_.x_ = payoff.cashPayoff();
        bachelier_.DxDforward_ = 0.0;
        bachelier_.DxDstrike_ = 0.0;
        bachelier_.y_ = bachelier_.beta_ = bachelier_.DbetaDd_ = bachelier_.DyDforward_ = 0.0;
    }

    void BachelierCalculator::Calculator::visit(AssetOrNothingPayoff& payoff) {
        // NPV = discount * forward * N(d) + stdDev * n(d) for calls, for puts  - forward * N(-d) +
        // stdDev * n(d)
        bachelier_.x_ = bachelier_.forward_;
        bachelier_.DxDforward_ = 1.0;
        bachelier_.DxDstrike_ = 0.0;
    }

    void BachelierCalculator::Calculator::visit(GapPayoff& payoff) {
        bachelier_.x_ = bachelier_.forward_ - payoff.secondStrike();
    }

    Real BachelierCalculator::value() const {
        Real result = discount_ * (x_ * alpha_ + y_ * beta_);
        return result;
    }

    Real BachelierCalculator::delta(Real spot) const {

        QL_REQUIRE(spot > 0.0, "positive spot value required: " << spot << " not allowed");
        Real DforwardDs = forward_ / spot;
        Real temp2 = deltaForward() * DforwardDs;
        return discount_ * temp2;
    }

    Real BachelierCalculator::deltaForward() const {
        Real DdDForward = 1.0 / stdDev_;
        Real DalphaDForward = DalphaDd_ * DdDForward;
        Real DbetaDForward = DbetaDd_ * DdDForward;
        Real temp2 =
            DalphaDForward * x_ + alpha_ * DxDforward_ + DbetaDForward * y_ + beta_ * DyDforward_;
        return discount_ * temp2;
    }

    Real BachelierCalculator::elasticity(Real spot) const {
        Real val = value();
        Real del = delta(spot);
        if (val > QL_EPSILON)
            return del / val * spot;
        else if (std::fabs(del) < QL_EPSILON)
            return 0.0;
        else if (del > 0.0)
            return QL_MAX_REAL;
        else
            return QL_MIN_REAL;
    }

    Real BachelierCalculator::elasticityForward() const {
        Real val = value();
        Real del = deltaForward();
        if (val > QL_EPSILON)
            return del / val * forward_;
        else if (std::fabs(del) < QL_EPSILON)
            return 0.0;
        else if (del > 0.0)
            return QL_MAX_REAL;
        else
            return QL_MIN_REAL;
    }

    Real BachelierCalculator::gamma(Real spot) const {
        QL_REQUIRE(spot > 0.0, "positive spot value required: " << spot << " not allowed");
        Real DforwardDs = forward_ / spot;
        Real temp2 = gammaForward() * DforwardDs * DforwardDs;
        return discount_ * temp2;
    }

    Real BachelierCalculator::gammaForward() const {
        Real DdDForward = 1.0 / stdDev_;
        Real temp2 = 2.0 * DalphaDd_ * DdDForward * DxDforward_ +
                     D2alphaD2d_ * DdDForward * DxDforward_ + DyDforward_ * D2betaD2d_ * DdDForward;
        return discount_ * temp2;
    }

    Real BachelierCalculator::theta(Real spot, Time maturity) const {

        QL_REQUIRE(maturity >= 0.0, "maturity (" << maturity << ") must be non-negative");
        if (close(maturity, 0.0))
            return 0.0;
        return -(std::log(discount_) * value() + std::log(forward_ / spot) * spot * delta(spot) +
                 0.5 * variance_ * gamma(spot)) /
               maturity;
    }

    Real BachelierCalculator::vega(Time maturity) const {
        QL_REQUIRE(maturity >= 0.0, "negative maturity not allowed");
        Real sigma = stdDev_ / std::sqrt(maturity);
        Real DdDsigma = -d_ / sigma;
        // actually DalphaDsigma / SQRT(T)
        Real DalphaDsigma = DalphaDd_ * DdDsigma;
        Real DbetaDsigma = DbetaDd_ * DdDsigma;
        Real DyDsigma = std::sqrt(maturity);
        Real temp2 = DalphaDsigma * x_ + DbetaDsigma * y_ + DyDsigma * beta_;
        return discount_ * temp2;
    }

    Real BachelierCalculator::rho(Time maturity) const {
        QL_REQUIRE(maturity >= 0.0, "negative maturity not allowed");
        Real DForwardDrho = maturity * forward_;
        return deltaForward() * DForwardDrho + maturity * value();
    }

    Real BachelierCalculator::dividendRho(Time maturity) const {
        QL_REQUIRE(maturity >= 0.0, "negative maturity not allowed");
        // actually DalphaDq / T
        Real DForwardDrho = -maturity * forward_;
        return deltaForward() * DForwardDrho - maturity * value();
    }

    Real BachelierCalculator::strikeSensitivity() const {
        Real DdDstrike = -stdDev_; // actually -1/stdDev
        Real DalphaDstrike = -DalphaDd_ / DdDstrike;
        Real DbetaDstrike = -DbetaDd_ / DdDstrike;

        Real temp2 =
            DalphaDstrike * x_ + alpha_ * DxDstrike_ + DbetaDstrike * y_ + beta_ * DyDstrike_;
        return discount_ * temp2;
    }

    Real BachelierCalculator::strikeGamma() const {
        Real DdDStrike = -stdDev_;
        Real DalphaDstrike = -DalphaDd_ / DdDStrike;
        Real DbetaDstrike = -DbetaDd_ / DdDStrike;

        Real temp2 = 2.0 * DalphaDstrike * DxDstrike_ + D2alphaD2d_ / DdDStrike * DxDstrike_ +
                     DyDstrike_ * D2betaD2d_ / DdDStrike;

        return discount_ * temp2;
    }
}
