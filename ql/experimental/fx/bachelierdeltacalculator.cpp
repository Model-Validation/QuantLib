/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2025 AcadiaSoft, Inc.

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

#include <ql/experimental/fx/bachelierdeltacalculator.hpp>

namespace QuantLib {

    BachelierDeltaCalculator::BachelierDeltaCalculator(
                        Option::Type ot,
                        DeltaVolQuote::DeltaType dt,
                        Real spot,
                        DiscountFactor dDiscount,   // domestic discount
                        DiscountFactor fDiscount,   // foreign  discount
                        Real stdDev):
    dt_(dt), ot_(ot),
    dDiscount_(dDiscount), fDiscount_(fDiscount),
    stdDev_(stdDev), spot_(spot),
    forward_(spot*fDiscount/dDiscount), phi_(Integer(ot)) {
    }


    Real BachelierDeltaCalculator::deltaFromStrike(Real strike) const {

        QL_REQUIRE(strike >= 0.0, "positive strike value required: " << strike << " not allowed");

        Real res = 0.0;

        switch (dt_) {
            case DeltaVolQuote::Spot:
                res = phi_ * fDiscount_ * cumD(strike);
                break;

            case DeltaVolQuote::Fwd:
                res = phi_ * cumD(strike);
                break;

            case DeltaVolQuote::PaSpot:
                QL_FAIL("Premium adjusted spot delta not implemented yet.");
                // its not stable around F = 0
                //res = fDiscount_ * (phi_ * strike * cumD(strike) - stdDev_ *  nD(strike)) / forward_;
                break;

            case DeltaVolQuote::PaFwd:
                // its not stable around F = 0
                QL_FAIL("Premium adjusted forward delta not implemented yet.");
                //res = fDiscount_ * (phi_ * strike * cumD(strike) - stdDev_ * nD(strike)) / forward_;
                break;

            default:
                QL_FAIL("invalid delta type");
        }
        return res;
    }

    Real BachelierDeltaCalculator::strikeFromDelta(Real delta) const {
        return(strikeFromDelta(delta, dt_));
    }

    Real BachelierDeltaCalculator::strikeFromDelta(Real delta,
                                               DeltaVolQuote::DeltaType dt)
                                                                        const{
        Real res=0.0;
        Real arg=0.0;
        InverseCumulativeNormal f;

        QL_REQUIRE(delta*phi_>=0.0, "Option type and delta are incoherent.");

        switch (dt) {
            case DeltaVolQuote::Spot:
                QL_REQUIRE(std::fabs(delta) <= fDiscount_, "Spot delta out of range.");
                arg = -phi_ * f(phi_ * delta / fDiscount_) * stdDev_;
                res = forward_ + arg;
                break;

            case DeltaVolQuote::Fwd:
                QL_REQUIRE(std::fabs(delta) <= 1.0, "Forward delta out of range.");
                arg = -phi_ * f(phi_ * delta) * stdDev_;
                res = forward_ + arg;
                break;

            case DeltaVolQuote::PaSpot:
            case DeltaVolQuote::PaFwd: 
                QL_FAIL("Premium adjusted delta not implemented yet.");
                break;
        }
        return res;
        
    }

    Real BachelierDeltaCalculator::atmStrike(DeltaVolQuote::AtmType atmT) const {

        Real res=0.0;

        switch(atmT) {
          case DeltaVolQuote::AtmSpot:
            res=spot_;
            break;

          case DeltaVolQuote::AtmDeltaNeutral:
            res = forward_; // call delta + put delta = 0 -> N(d) + (N(d) - 1) = 0 -> N(d) = 0.5 -> d=0 -> K=F
            break;

          case DeltaVolQuote::AtmFwd:
            res=forward_; 
            break;

          case DeltaVolQuote::AtmGammaMax: case DeltaVolQuote::AtmVegaMax:
            res=forward_; ////n(d) max at 0 -> K = F
            break;

          case DeltaVolQuote::AtmPutCall50:
            QL_REQUIRE(dt_==DeltaVolQuote::Fwd,
                       "|PutDelta|=CallDelta=0.50 only possible for forward delta.");
            res=forward_;
            break;

          default:
            QL_FAIL("invalid atm type");
        }

        return res;
    }


    Real BachelierDeltaCalculator::cumD(Real strike) const {

        Real d=0.0;
        Real cum_d1_pos_ = 1.0; // N(d1)
        Real cum_d1_neg_ = 0.0; // N(-d1)

        CumulativeNormalDistribution f;

        if (stdDev_>=QL_EPSILON) {
            d = (forward_ - strike) / stdDev_;
            return f(phi_*d);
        } else {
            if (forward_<strike) {
                cum_d1_pos_ = 0.0;
                cum_d1_neg_ = 1.0;
            } else if(forward_==strike){
                return 0.5;
            }
        }

        if (phi_>0) { // if Call
            return cum_d1_pos_;
        } else {
            return cum_d1_neg_;
        }
    }


    Real BachelierDeltaCalculator::nD(Real strike) const {

        Real d=0.0;
        Real n_d = 0.0; // n(d)
        if (stdDev_>=QL_EPSILON){
            d = (forward_ - strike) / stdDev_;
            CumulativeNormalDistribution f;
            n_d = f.derivative(d);
        }
        return n_d;
    }

    void BachelierDeltaCalculator::setDeltaType(DeltaVolQuote::DeltaType dt){
        dt_=dt;
    }

    void BachelierDeltaCalculator::setOptionType(Option::Type ot){
        ot_=ot;
        phi_=Integer(ot_);
    }


    // helper classes

    BachelierDeltaPremiumAdjustedSolverClass::BachelierDeltaPremiumAdjustedSolverClass(
                        Option::Type ot,
                        DeltaVolQuote::DeltaType dt,
                        Real spot,
                        DiscountFactor dDiscount,   // domestic discount
                        DiscountFactor fDiscount,   // foreign  discount
                        Real stdDev,
                        Real delta):
    bdc_(ot,dt,spot,dDiscount,fDiscount,stdDev), delta_(delta) {}


    Real BachelierDeltaPremiumAdjustedSolverClass::operator()(Real strike) const {
        return bdc_.deltaFromStrike(strike)-delta_;
    }


}
