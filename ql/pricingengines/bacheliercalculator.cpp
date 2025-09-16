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

#include <ql/pricingengines/bacheliercalculator.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/math/comparison.hpp>

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
    : strike_(p->strike()), forward_(forward), stdDev_(stdDev),
      discount_(discount), variance_(stdDev*stdDev) {
        initialize(p);
    }

    BachelierCalculator::BachelierCalculator(Option::Type optionType,
                                     Real strike,
                                     Real forward,
                                     Real stdDev,
                                     Real discount)
    : strike_(strike), forward_(forward), stdDev_(stdDev),
      discount_(discount), variance_(stdDev*stdDev) {
        initialize(ext::shared_ptr<StrikedTypePayoff>(new
            PlainVanillaPayoff(optionType, strike)));
    }

    void BachelierCalculator::initialize(const ext::shared_ptr<StrikedTypePayoff>& p) {
        QL_REQUIRE(stdDev_>=0.0,
                   "stdDev (" << stdDev_ << ") must be non-negative");
        QL_REQUIRE(discount_>0.0,
                   "discount (" << discount_ << ") must be positive");

        if (stdDev_>=QL_EPSILON) {
            d_ = (forward_ - strike_) / stdDev_;
            auto f = CumulativeNormalDistribution();
            cum_d_ = f(d_);
            n_d_ = f.derivative(d_);

        } else {
            if (close(forward_, strike_)) {
                d_ = 0.0;
                cum_d_ = 0.5;
                n_d_ = M_SQRT_2 * M_1_SQRTPI;
            } else if (forward_>strike_) {
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

        // the following one will probably disappear as soon as
        // super-share will be properly handled
        DxDs_ = 0.0;

        // this part is always executed.
        // in case of plain-vanilla payoffs, it is also the only part
        // which is executed.
        switch (p->optionType()) {
            case Option::Call:
                alpha_ = cum_d_;       //  N(d)
                DalphaDd_ = n_d_;      //  n(d)
                beta_ = n_d_;          // n(d)
                DbetaDd_ = -d_ * n_d_; // -d * n(d)
                break;
            case Option::Put:
                alpha_ = -1.0 + cum_d_; // -N(-d)
                DalphaDd_ = n_d_;       //  n(d)
                beta_ = n_d_;           //  n(d)
                DbetaDd_ = -d_ * n_d_;  // -d * n(d)
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
        bachelier_.y_ = bachelier_.beta_ = bachelier_.DbetaDd_ = bachelier_.DyDforward_ = 0.0;
        bachelier_.x_ = payoff.cashPayoff();
        bachelier_.DxDstrike_ = 0.0;
        switch (payoff.optionType()) {
            case Option::Call:
                bachelier_.alpha_ = bachelier_.cum_d_;
                bachelier_.DalphaDd_ = bachelier_.n_d_;
                break;
            case Option::Put:
                bachelier_.alpha_ = -1.0 + bachelier_.cum_d_;
                bachelier_.DalphaDd_ = bachelier_.n_d_;
                break;
            default:
                QL_FAIL("invalid option type");
        }
    }

    void BachelierCalculator::Calculator::visit(AssetOrNothingPayoff& payoff) {
        // NPV = discount * forward * N(d) + stdDev * n(d) for calls, for puts  - forward * N(-d) +
        // stdDev * n(d)
        bachelier_.x_ = bachelier_.forward_;
        bachelier_.DxDforward_ = 1.0;
        bachelier_.y_ = bachelier_.stdDev_;
        bachelier_.DyDforward_ = 0.0;

        switch (payoff.optionType()) {
            case Option::Call:
                bachelier_.alpha_ = bachelier_.cum_d_;                  //  N(d)
                bachelier_.DalphaDd_ = bachelier_.n_d_;                 //  n(d)
                bachelier_.beta_ = bachelier_.n_d_;                     // n(d)
                bachelier_.DbetaDd_ = -bachelier_.d_ * bachelier_.n_d_; // -d * n(d)
                break;
            case Option::Put:
                bachelier_.alpha_ = -1.0 + bachelier_.cum_d_;           // -N(-d)
                bachelier_.DalphaDd_ = bachelier_.n_d_;                 //  n(d)
                bachelier_.beta_ = bachelier_.n_d_;                     //  n(d)
                bachelier_.DbetaDd_ = -bachelier_.d_ * bachelier_.n_d_; // -d * n(d)
                break;
            default:
                QL_FAIL("invalid option type");
        }
    }

    void BachelierCalculator::Calculator::visit(GapPayoff& payoff) {
        bachelier_.x_ =  bachelier_.forward_ - payoff.secondStrike();
        bachelier_.DxDforward_ = 1.0;
        bachelier_.DxDstrike_ = 0.0;
        bachelier_.y_ = bachelier_.stdDev_;
        bachelier_.DyDforward_ = 0.0;

        switch (payoff.optionType()) {
            case Option::Call:
                bachelier_.alpha_ = bachelier_.cum_d_;                  //  N(d)
                bachelier_.DalphaDd_ = bachelier_.n_d_;                 //  n(d)
                bachelier_.beta_ = bachelier_.n_d_;                     // n(d)
                bachelier_.DbetaDd_ = -bachelier_.d_ * bachelier_.n_d_; // -d * n(d)
                break;
            case Option::Put:
                bachelier_.alpha_ = -1.0 + bachelier_.cum_d_;           // -N(-d)
                bachelier_.DalphaDd_ = bachelier_.n_d_;                 //  n(d)
                bachelier_.beta_ = bachelier_.n_d_;                     //  n(d)
                bachelier_.DbetaDd_ = -bachelier_.d_ * bachelier_.n_d_; // -d * n(d)
                break;
            default:
                QL_FAIL("invalid option type");
        }
    }

    Real BachelierCalculator::value() const {
        Real result = discount_ * (x_ * alpha_ + y_ * beta_);
        return result;
    }

    Real BachelierCalculator::delta(Real spot) const {

        QL_REQUIRE(spot > 0.0, "positive spot value required: " << spot << " not allowed");

        Real DforwardDs = forward_ / spot;
        Real DdDForward = 1.0 / stdDev_;
        Real DalphaDs = DalphaDd_ * DdDForward * DforwardDs;

        Real DbetaDs = DbetaDd_ * DdDForward * DforwardDs;

        Real temp2 = DalphaDs * x_ + alpha_ * DxDforward_ * DforwardDs + DbetaDs * y_ +
                     beta_ * DyDforward_ * DforwardDs;

        return discount_ * temp2;
    }

    Real BachelierCalculator::deltaForward() const {

        Real DdDForward = 1.0 / stdDev_;
        Real DalphaDs = DalphaDd_ * DdDForward;

        Real DbetaDs = DbetaDd_ * DdDForward;
        Real temp2 = DalphaDs * x_ + alpha_ * DxDforward_ + DbetaDs * y_ + beta_ * DyDforward_;

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
        return 0.0;
        QL_REQUIRE(spot > 0.0, "positive spot value required: " << spot << " not allowed");

        Real DforwardDs = forward_ / spot;

        Real temp = stdDev_ * spot;
        Real DalphaDs = DalphaDd_ / temp;
        Real DbetaDs = DbetaDd_ / temp;

        Real D2alphaDs2 = -DalphaDs / spot * (1 + d1_ / stdDev_);
        Real D2betaDs2 = -DbetaDs / spot * (1 + d2_ / stdDev_);

        Real temp2 = D2alphaDs2 * forward_ + 2.0 * DalphaDs * DforwardDs + D2betaDs2 * x_ +
                     2.0 * DbetaDs * DxDs_;

        return discount_ * temp2;
    }

    Real BachelierCalculator::gammaForward() const {
        return 0.0;
        Real DdDForward = 1.0 / stdDev_;
        Real temp2 = 2.0 * DalphaDd_ * DdDForward * DxDforward_ - x_ * d_ * DalphaDd_ * DdDForward;


        return discount_ * temp2;
    }

    Real BachelierCalculator::theta(Real spot, Time maturity) const {
        return 0.0;
        QL_REQUIRE(maturity >= 0.0, "maturity (" << maturity << ") must be non-negative");
        if (close(maturity, 0.0))
            return 0.0;
        return -(std::log(discount_) * value() + std::log(forward_ / spot) * spot * delta(spot) +
                 0.5 * variance_ * spot * spot * gamma(spot)) /
               maturity;
    }

    Real BachelierCalculator::vega(Time maturity) const {
        return 0.0;
        QL_REQUIRE(maturity >= 0.0, "negative maturity not allowed");

        Real temp = std::log(strike_ / forward_) / variance_;
        // actually DalphaDsigma / SQRT(T)
        Real DalphaDsigma = DalphaDd_ * (temp + 0.5);
        Real DbetaDsigma = DbetaDd_ * (temp - 0.5);

        Real temp2 = DalphaDsigma * forward_ + DbetaDsigma * x_;

        return discount_ * std::sqrt(maturity) * temp2;
    }

    Real BachelierCalculator::rho(Time maturity) const {
        return 0.0;
        QL_REQUIRE(maturity >= 0.0, "negative maturity not allowed");

        // actually DalphaDr / T
        Real DalphaDr = DalphaDd_ / stdDev_;
        Real DbetaDr = DbetaDd_ / stdDev_;
        Real temp = DalphaDr * forward_ + alpha_ * forward_ + DbetaDr * x_;

        return maturity * (discount_ * temp - value());
    }

    Real BachelierCalculator::dividendRho(Time maturity) const {
        return 0.0;
        QL_REQUIRE(maturity >= 0.0, "negative maturity not allowed");

        // actually DalphaDq / T
        Real DalphaDq = -DalphaDd_ / stdDev_;
        Real DbetaDq = -DbetaDd_ / stdDev_;

        Real temp = DalphaDq * forward_ - alpha_ * forward_ + DbetaDq * x_;

        return maturity * discount_ * temp;
    }

    Real BachelierCalculator::strikeSensitivity() const {
        return 0.0;
        Real temp = stdDev_ * strike_;
        Real DalphaDstrike = -DalphaDd_ / temp;
        Real DbetaDstrike = -DbetaDd_ / temp;

        Real temp2 = DalphaDstrike * forward_ + DbetaDstrike * x_ + beta_ * DxDstrike_;

        return discount_ * temp2;
    }

    Real BachelierCalculator::strikeGamma() const {
        return 0.0;
        Real temp = stdDev_ * strike_;
        Real DalphaDstrike = -DalphaDd_ / temp;
        Real DbetaDstrike = -DbetaDd_ / temp;

        Real D2alphaD2strike = -DalphaDstrike / strike_ * (1 - d1_ / stdDev_);
        Real D2betaD2strike = -DbetaDstrike / strike_ * (1 - d2_ / stdDev_);

        Real temp2 =
            D2alphaD2strike * forward_ + D2betaD2strike * x_ + 2.0 * DbetaDstrike * DxDstrike_;

        return discount_ * temp2;
    }
}
