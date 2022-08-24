/*
 Copyright (C) 2022 Skandinaviska Enskilda Banken AB (publ)

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
#include <ql/pricingengines/asian/discretearithmeticasianlevyengine.hpp>
#include <ql/pricingengines/blackformula.hpp>

using namespace QuantLib;

void DiscreteArithmeticAsianLevyEngine::calculate() const {

    // Enforce a few required things
    QL_REQUIRE(arguments_.exercise->type() == Exercise::European, "not a European Option");
    QL_REQUIRE(arguments_.averageType == Average::Type::Arithmetic,
               "must be Arithmetic Average::Type");

    // Calculate the accrued portion
    Real runningAccumulator = arguments_.runningAccumulator;
    Size pastFixings = arguments_.pastFixings;
    Real accruedAverage = 0;
    Size m = arguments_.fixingDates.size() + pastFixings;
    if (pastFixings != 0) {
        accruedAverage = runningAccumulator / m;
    }

    // Populate some additional results that don't change
    Real discount = process_->riskFreeRate()->discount(arguments_.exercise->lastDate());
    results_.additionalResults["discount"] = discount;
    results_.additionalResults["accrued"] = accruedAverage;

    ext::shared_ptr<PlainVanillaPayoff> payoff =
        ext::dynamic_pointer_cast<PlainVanillaPayoff>(arguments_.payoff);
    QL_REQUIRE(payoff, "non-plain payoff given");

    // TODO: If not model dependent, return early.
    /*if (pastFixings > 0) {

        if (accruedAverage > 1.0 * arguments_.fixingDates.size() / pastFixings * payoff->strike()) {
            if (payoff->optionType() == Option::Type::Call) {
                results_.value = 1010101;
            } else if (payoff->optionType() == Option::Type::Put) {
                results_.value = 0;
                return;
            } else {
                QL_FAIL("unexpected option type " << payoff->optionType());
            }
        }
    }*/

    // We will read the volatility off the surface at the effective strike
    Real effectiveStrike = payoff->strike() - accruedAverage;
    results_.additionalResults["strike"] = payoff->strike();
    results_.additionalResults["effective_strike"] = effectiveStrike;
    // We should only get this far when the effectiveStrike > 0 but will check anyway
    QL_REQUIRE(effectiveStrike > 0.0, "expected effectiveStrike to be positive");

    // Expected value of the non-accrued portion of the average prices
    // In general, m will equal n below if there is no accrued. If accrued, m > n.
    Real EA = 0.0;
    std::vector<Real> forwards;
    std::vector<Time> times;
    for (const auto& fd : arguments_.fixingDates) {
        Real spot = process_->stateVariable()->value();
        DiscountFactor dividendDiscount = process_->dividendYield()->discount(fd);
        DiscountFactor riskFreeDiscountForFwdEstimation = process_->riskFreeRate()->discount(fd);

        forwards.push_back(spot * dividendDiscount / riskFreeDiscountForFwdEstimation);
        times.push_back(process_->blackVolatility()->timeFromReference(fd));

        EA += forwards.back();
    }
    EA /= m;

    // Calculate the expected value of A^2

    Size n = forwards.size();
    Time t1 = times.front();
    Real tn = times.back();
    Real vol = process_->blackVolatility()->blackVol(tn, effectiveStrike);

    Real h = (tn - t1) / (n - 1);
    Real exp_1 = (1 - exp(pow(vol, 2) * h * n)) / (1 - exp(pow(vol, 2) * h));
    Real exp_2 = 2 / (1 - exp(pow(vol, 2) * h));
    // Real EA2 = (pow(EA, 2) * exp(pow(vol, 2) * t1)) / pow(n, 2) * (exp_1 + exp_2 * (n - exp_1));

    Real EA2 = 0.0;
    for (Size i = 0; i < n; ++i) {
        EA2 += forwards[i] * forwards[i] * exp(vol * vol * times[i]);
        for (Size j = 0; j < i; ++j) {
            EA2 += 2 * forwards[i] * forwards[j] * exp(vol * vol * times[j]);
        }
    }
    EA2 /= m * m;

    // Asian volatility
    Real sigma = sqrt(log(EA2 / (EA * EA)) / tn);

    // Populate results
    results_.value =
        blackFormula(payoff->optionType(), effectiveStrike, EA, sigma * sqrt(tn), discount);
    results_.additionalResults["forward"] = EA;
    results_.additionalResults["exp_A_2"] = EA2;
    results_.additionalResults["t1"] = t1;
    results_.additionalResults["tte"] = tn;
    results_.additionalResults["vol"] = vol;
    results_.additionalResults["sigma"] = sigma;
    results_.additionalResults["times"] = times;
    results_.additionalResults["forwards"] = forwards;
}
