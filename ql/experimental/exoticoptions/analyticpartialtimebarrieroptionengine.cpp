/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2014 Master IMAFA - Polytech'Nice Sophia - Université de Nice Sophia Antipolis
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

#include <ql/exercise.hpp>
#include <ql/experimental/exoticoptions/analyticpartialtimebarrieroptionengine.hpp>
#include <ql/math/distributions/bivariatenormaldistribution.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
#include <ql/quotes/simplequote.hpp>
#include <utility>

namespace QuantLib {

    AnalyticPartialTimeBarrierOptionEngine::AnalyticPartialTimeBarrierOptionEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process)
    : process_(std::move(process)) {
        registerWith(process_);
    }

    void AnalyticPartialTimeBarrierOptionEngine::calculate() const {
        ext::shared_ptr<PlainVanillaPayoff> payoff =
            ext::dynamic_pointer_cast<PlainVanillaPayoff>(arguments_.payoff);
        QL_REQUIRE(payoff, "non-plain payoff given");
        QL_REQUIRE(payoff->strike()>0.0,
                   "strike must be positive");

        Real spot = process_->x0();
        QL_REQUIRE(spot > 0.0, "negative or null underlying given");

        results_.additionalResults["residualTime"] = residualTime();
        results_.additionalResults["coverEventTime"] = coverEventTime();
        results_.additionalResults["volatility(t1)"] = volatility(coverEventTime());
        results_.additionalResults["volatility(T2)"] = volatility(residualTime());
        results_.additionalResults["barrier"] = barrier();
        results_.additionalResults["strike"] = strike();
        results_.additionalResults["riskFreeRate"] = riskFreeRate();
        results_.additionalResults["dividendYield"] = dividendYield();

        PartialBarrier::Type barrierType = arguments_.barrierType;
        PartialBarrier::Range barrierRange = arguments_.barrierRange;

        switch (payoff->optionType()) {
          //Call Option
          case Option::Call:
            switch (barrierType) {
              case PartialBarrier::DownOut:
                switch (barrierRange) {
                  case PartialBarrier::Start:
                    results_.value = CA(1);
                    break;
                  case PartialBarrier::EndB1:
                    results_.value = CoB1();
                    break;
                  case PartialBarrier::EndB2:
                    results_.value = CoB2(PartialBarrier::DownOut);
                    break;
                  default:
                    QL_FAIL("invalid barrier range");
                }
                break;

              case PartialBarrier::DownIn:
                switch (barrierRange) {
                  case PartialBarrier::Start:
                    results_.value = CIA(1);
                    break;
                  case PartialBarrier::EndB1:
                    results_.value = underlyingOption()->NPV() - CoB1();
                    break;
                  case PartialBarrier::EndB2:
                    results_.value = underlyingOption()->NPV() - CoB2(PartialBarrier::DownOut);
                    break;
                  default:
                    QL_FAIL("invalid barrier range");
                }
                break;

              case PartialBarrier::UpOut:
                switch (barrierRange) {
                  case PartialBarrier::Start:
                    results_.value = CA(-1);
                    break;
                  case PartialBarrier::EndB1:
                    results_.value = CoB1();
                    break;
                  case PartialBarrier::EndB2:
                    results_.value = CoB2(PartialBarrier::UpOut);
                    break;
                  default:
                    QL_FAIL("invalid barrier range");
                }
                break;

              case PartialBarrier::UpIn:
                switch (barrierRange) {
                  case PartialBarrier::Start:
                    results_.value = CIA(-1);
                    break;
                  case PartialBarrier::EndB1:
                    results_.value = underlyingOption()->NPV() - CoB1();
                    break;
                  case PartialBarrier::EndB2:
                    results_.value = underlyingOption()->NPV() - CoB2(PartialBarrier::UpOut);
                    break;
                  default:
                    QL_FAIL("invalid barrier range");
                }
                break;
              default:
                QL_FAIL("unknown barrier type");
            }
            break;

          case Option::Put: {
            /*
            We use two approaches here to properly value put options. 
            Specifically, out-style options are directly transformed into call options, 
            while in-style options are valued as a vanilla put minus a corresponding out-style option.
            */
            ext::shared_ptr<EuropeanExercise> exercise =
                ext::dynamic_pointer_cast<EuropeanExercise>(arguments_.exercise);

            // For the put-call tranformation, we invert the barrier
            // direction. Up becomes down, down becomes up
            if (barrierType == PartialBarrier::DownIn || barrierType == PartialBarrier::UpIn) {
              auto vanillaPutOption = underlyingOption();
              PartialBarrier::Type barrierTypePut =
                  (barrierType == PartialBarrier::DownIn ? PartialBarrier::DownOut :
                                                            PartialBarrier::UpOut);
              ext::shared_ptr<PlainVanillaPayoff> payoffPut=
                  ext::make_shared<PlainVanillaPayoff>(Option::Put, strike());

              PartialTimeBarrierOption ptBarrierPutOut(
                  barrierTypePut, barrierRange, barrier(), rebate(),
                  arguments_.coverEventDate, payoffPut, exercise);
              ptBarrierPutOut.setPricingEngine(
                  ext::make_shared<AnalyticPartialTimeBarrierOptionEngine>(process_));

              results_.value = vanillaPutOption->NPV() - ptBarrierPutOut.NPV();
              for (auto res : ptBarrierPutOut.additionalResults()) {
                results_.additionalResults["PUT_" + res.first] = res.second;
              }
            } else if (barrierType == PartialBarrier::DownOut ||
                       barrierType == PartialBarrier::UpOut) {
              PartialBarrier::Type barrierTypeCall =
                  (barrierType == PartialBarrier::DownOut ? PartialBarrier::UpOut :
                                                            PartialBarrier::DownOut);
              ext::shared_ptr<PlainVanillaPayoff> payoffCall =
                  ext::make_shared<PlainVanillaPayoff>(Option::Call, underlying());

              PartialTimeBarrierOption ptBarrierCall(
                  barrierTypeCall, barrierRange, underlying() * strike() / barrier(), rebate(),
                  arguments_.coverEventDate, payoffCall, exercise);

              Handle<Quote> strikeAsSpot = Handle<Quote>(ext::make_shared<SimpleQuote>(strike()));
              // New process with flipped rfr/div term structures
              ext::shared_ptr<GeneralizedBlackScholesProcess> gbsp =
                  ext::make_shared<GeneralizedBlackScholesProcess>(
                      strikeAsSpot, process_->riskFreeRate(), process_->dividendYield(),
                      process_->blackVolatility());
              ptBarrierCall.setPricingEngine(
                  ext::make_shared<AnalyticPartialTimeBarrierOptionEngine>(gbsp));

              results_.value = ptBarrierCall.NPV();
              for (auto res : ptBarrierCall.additionalResults()) {
                results_.additionalResults["CALL_" + res.first] = res.second;
              }
            } else {
              QL_FAIL("unknown barrier type");
            }

            break;
          }
          default:
            QL_FAIL("unknown option type");
        }
    }

    Real AnalyticPartialTimeBarrierOptionEngine::CoB2(
                                      PartialBarrier::Type barrierType) const {
        Real result = 0.0;
        Real b = riskFreeRate()-dividendYield();
        if (strike()<barrier()){
            switch (barrierType) {
              case PartialBarrier::DownOut:
                result = underlying()*std::exp((b-riskFreeRate())*residualTime());
                result *= (M(g1(),e1(),rho())-HS(underlying(),barrier(),2*(mu()+1))*M(g3(),-e3(),-rho()));
                result -= strike()*std::exp(-riskFreeRate()*residualTime())*(M(g2(),e2(),rho())-HS(underlying(),barrier(),2*mu())*M(g4(),-e4(),-rho()));
                return result;

              case PartialBarrier::UpOut:
                result = underlying()*std::exp((b-riskFreeRate())*residualTime());
                result *= (M(-g1(),-e1(),rho())-HS(underlying(),barrier(),2*(mu()+1))*M(-g3(),e3(),-rho()));
                result -= strike()*std::exp(-riskFreeRate()*residualTime())*(M(-g2(),-e2(),rho())-HS(underlying(),barrier(),2*mu())*M(-g4(),e4(),-rho()));
                result -= underlying()*std::exp((b-riskFreeRate())*residualTime())*(M(-d1(),-e1(),rho())-HS(underlying(),barrier(),2*(mu()+1))*M(e3(),-f1(),-rho()));
                result += strike()*std::exp(-riskFreeRate()*residualTime())*(M(-d2(),-e2(),rho())-HS(underlying(),barrier(),2*mu())*M(e4(),-f2(),-rho()));
                return result;

              default:
                QL_FAIL("invalid barrier type");
            }
        } else {
            switch (barrierType) {
              case PartialBarrier::DownOut:
                return CoB1();
              case PartialBarrier::UpOut:
                return 0.0;
              default:
                QL_FAIL("invalid barrier type");
            }
        }
    }

    Real AnalyticPartialTimeBarrierOptionEngine::CoB1() const {
        Real result = 0.0;
        Real b = riskFreeRate()-dividendYield();
        if (strike()>barrier()) {
            result = underlying()*std::exp((b-riskFreeRate())*residualTime());
            result *= (M(d1(),e1(),rho())-HS(underlying(),barrier(),2*(mu()+1))*M(f1(),-e3(),-rho()));
            result -= (strike()*std::exp(-riskFreeRate()*residualTime()))*(M(d2(),e2(),rho())-HS(underlying(),barrier(),2*mu())*M(f2(),-e4(),-rho()));
            return result;
        } else {
            Real S1 = underlying()*std::exp((b-riskFreeRate())*residualTime());
            Real X1 = (strike()*std::exp(-riskFreeRate()*residualTime()));
            Real HS1 = HS(underlying(),barrier(),2*(mu()+1));
            Real HS2 = HS(underlying(), barrier(), 2 * mu());
            result = S1;
            result *= (M(-g1(),-e1(),rho())-HS1*M(-g3(),e3(),-rho()));
            result -= X1*(M(-g2(), -e2(), rho()) - HS2*M(-g4(), e4(), -rho()));
            result -= S1*(M(-d1(), -e1(), rho()) - HS1*M(-f1(), e3(), -rho()));
            result += X1*(M(-d2(), -e2(), rho()) - HS2*M(-f2(), e4(), -rho()));
            result += S1*(M(g1(), e1(), rho()) - HS1*M(g3(), -e3(), -rho()));
            result -= X1*(M(g2(), e2(), rho()) - HS2*M(g4(), -e4(), -rho()));
            return result;
        }
    }

    // eta = -1: Up-and-In Call
    // eta =  1: Down-and-In Call
    Real AnalyticPartialTimeBarrierOptionEngine::CIA(Integer eta) const {
        return underlyingOption()->NPV() - CA(eta);
    }

    Real AnalyticPartialTimeBarrierOptionEngine::CA(Integer eta) const {
        //Partial-Time-Start- OUT  Call Option calculation
        Real b = riskFreeRate()-dividendYield();
        Real result;
        result = underlying()*std::exp((b-riskFreeRate())*residualTime());
        result *= (M(d1(),eta*e1(),eta*rho())-HS(underlying(),barrier(),2*(mu()+1))*M(f1(),eta*e3(),eta*rho()));
        result -= (strike()*std::exp(-riskFreeRate()*residualTime())*(M(d2(),eta*e2(),eta*rho())-HS(underlying(),barrier(),2*mu())*M(f2(),eta*e4(),eta*rho())));
        return result;
    }

    Real AnalyticPartialTimeBarrierOptionEngine::underlying() const {
        return process_->x0();
    }

    Real AnalyticPartialTimeBarrierOptionEngine::strike() const {
        ext::shared_ptr<PlainVanillaPayoff> payoff =
            ext::dynamic_pointer_cast<PlainVanillaPayoff>(arguments_.payoff);
        QL_REQUIRE(payoff, "non-plain payoff given");
        return payoff->strike();
    }

    Time AnalyticPartialTimeBarrierOptionEngine::residualTime() const {
        return process_->time(arguments_.exercise->lastDate());
    }

    Time AnalyticPartialTimeBarrierOptionEngine::coverEventTime() const {
        return process_->time(arguments_.coverEventDate);
    }

    Volatility AnalyticPartialTimeBarrierOptionEngine::volatility(Time t) const {
        return process_->blackVolatility()->blackVol(t, strike());
    }

    Real AnalyticPartialTimeBarrierOptionEngine::stdDeviation() const {
        Time T = residualTime();
        return volatility(T) * std::sqrt(T);
    }

    Real AnalyticPartialTimeBarrierOptionEngine::barrier() const {
        return arguments_.barrier;
    }

    Real AnalyticPartialTimeBarrierOptionEngine::rebate() const {
        return arguments_.rebate;
    }

    Rate AnalyticPartialTimeBarrierOptionEngine::riskFreeRate() const {
        return process_->riskFreeRate()->zeroRate(residualTime(), Continuous,
                                                  NoFrequency);
    }

    DiscountFactor AnalyticPartialTimeBarrierOptionEngine::riskFreeDiscount() const {
        return process_->riskFreeRate()->discount(residualTime());
    }

    Rate AnalyticPartialTimeBarrierOptionEngine::dividendYield() const {
        return process_->dividendYield()->zeroRate(residualTime(), Continuous,
                                                   NoFrequency);
    }

    DiscountFactor AnalyticPartialTimeBarrierOptionEngine::dividendDiscount() const {
        return process_->dividendYield()->discount(residualTime());
    }


    Real AnalyticPartialTimeBarrierOptionEngine::f1() const {
        Real S = underlying();
        Real T = residualTime();
        Real sigma = volatility(T);
        return (std::log(S / strike()) + 2 * std::log(barrier() / S) + ((riskFreeRate()-dividendYield()) + (std::pow(sigma, 2) / 2))*T) / (sigma*std::sqrt(T));
    }

    Real AnalyticPartialTimeBarrierOptionEngine::f2() const {
        Time T = residualTime();
        return f1() - volatility(T)*std::sqrt(T);
    }

    Real AnalyticPartialTimeBarrierOptionEngine::M(Real a,Real b,Real rho) const {
        BivariateCumulativeNormalDistributionDr78 CmlNormDist(rho);
        return CmlNormDist(a,b);
    }

    Real AnalyticPartialTimeBarrierOptionEngine::rho() const {
        return std::sqrt(coverEventTime()/residualTime());
    }

    Rate AnalyticPartialTimeBarrierOptionEngine::mu() const {
        Volatility vol = volatility(coverEventTime());
        return ((riskFreeRate() - dividendYield()) - (vol * vol) / 2) / (vol * vol);
    }

    Real AnalyticPartialTimeBarrierOptionEngine::d1() const {
        Real b = riskFreeRate()-dividendYield();
        Time T2 = residualTime();
        Volatility vol = volatility(T2);
        return (std::log(underlying()/strike())+(b+vol*vol/2)*T2)/(std::sqrt(T2)*vol);
    }

    Real AnalyticPartialTimeBarrierOptionEngine::d2() const {
        Time T2 = residualTime();
        Volatility vol = volatility(T2);
        return d1() - vol*std::sqrt(T2);
    }

    Real AnalyticPartialTimeBarrierOptionEngine::e1() const {
        Real b = riskFreeRate()-dividendYield();
        Time T1 = coverEventTime();
        Volatility vol = volatility(T1);
        return (std::log(underlying()/barrier())+(b+vol*vol/2)*T1)/(std::sqrt(T1)*vol);
    }

    Real AnalyticPartialTimeBarrierOptionEngine::e2() const {
        Time T1 = coverEventTime();
        Volatility vol = volatility(T1);
        return e1() - vol*std::sqrt(T1);
    }

    Real AnalyticPartialTimeBarrierOptionEngine::e3() const {
        Time T1 = coverEventTime();
        Real vol = volatility(T1);
        return e1()+(2*std::log(barrier()/underlying()) /(vol*std::sqrt(T1)));
    }

    Real AnalyticPartialTimeBarrierOptionEngine::e4() const {
        Time t = coverEventTime();
        return e3()-volatility(t)*std::sqrt(t);
    }

    Real AnalyticPartialTimeBarrierOptionEngine::g1() const {
        Real b = riskFreeRate()-dividendYield();
        Time T2 = residualTime();
        Volatility vol = volatility(T2);
        return (std::log(underlying()/barrier())+(b+vol*vol/2)*T2)/(std::sqrt(T2)*vol);
    }

    Real AnalyticPartialTimeBarrierOptionEngine::g2() const {
        Time T2 = residualTime();
        Volatility vol = volatility(T2);
        return g1() - vol*std::sqrt(T2);
    }

    Real AnalyticPartialTimeBarrierOptionEngine::g3() const {
        Time T2 = residualTime();
        Real vol = volatility(T2);
        return g1()+(2*std::log(barrier()/underlying()) /(vol*std::sqrt(T2)));
    }

    Real AnalyticPartialTimeBarrierOptionEngine::g4() const {
        Time T2 = residualTime();
        Real vol = volatility(T2);
        return g3()-vol*std::sqrt(T2);
    }

    Real AnalyticPartialTimeBarrierOptionEngine::HS(Real S, Real H, Real power) const {
        return std::pow((H/S),power);
    }

    const ext::shared_ptr<VanillaOption>
    AnalyticPartialTimeBarrierOptionEngine::underlyingOption() const {
        ext::shared_ptr<EuropeanExercise> exercise =
            ext::dynamic_pointer_cast<EuropeanExercise>(arguments_.exercise);

        ext::shared_ptr<PlainVanillaPayoff> payoff =
            ext::dynamic_pointer_cast<PlainVanillaPayoff>(arguments_.payoff);

        ext::shared_ptr<VanillaOption> europeanOption= ext::make_shared<VanillaOption>(payoff, exercise);

        europeanOption->setPricingEngine(ext::make_shared<AnalyticEuropeanEngine>(process_));
        return europeanOption;
    }


}

