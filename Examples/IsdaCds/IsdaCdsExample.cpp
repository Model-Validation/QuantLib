/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2024 QuantLib Contributors

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.
*/

/*! \file IsdaCdsExample.cpp
    \brief Example demonstrating ISDA CDS interface usage
*/

#include <ql/qldefines.hpp>
#if !defined(BOOST_ALL_NO_LIB) && defined(BOOST_MSVC)
#  include <ql/auto_link.hpp>
#endif

#include <ql/experimental/isda/isdacdsinterface.hpp>
#include <ql/instruments/creditdefaultswap.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/credit/flathazardrate.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/settings.hpp>

#include <iostream>
#include <iomanip>

using namespace QuantLib;
using namespace std;

int main() {
    try {
        
        cout << "\n=== ISDA CDS Interface Example ===" << endl;
        
        // Set up market data
        Calendar calendar = TARGET();
        Date today(15, May, 2024);
        Settings::instance().evaluationDate() = today;
        
        // Discount curve (flat 2% rate)
        Handle<Quote> riskFreeRate(ext::make_shared<SimpleQuote>(0.02));
        Handle<YieldTermStructure> discountCurve(
            ext::make_shared<FlatForward>(today, riskFreeRate, Actual365Fixed()));
        
        // Credit curve (flat 300bp hazard rate)
        Handle<Quote> hazardRate(ext::make_shared<SimpleQuote>(0.03));
        Handle<DefaultProbabilityTermStructure> creditCurve(
            ext::make_shared<FlatHazardRate>(today, hazardRate, Actual365Fixed()));
        
        // CDS parameters
        Real recoveryRate = 0.4;  // 40% recovery
        Real notional = 10000000; // $10MM
        Rate spread = 0.01;       // 100bp spread
        
        // Create CDS schedule (quarterly payments, 5 years)
        Date startDate = calendar.advance(today, 1, Days);
        Date endDate = calendar.advance(today, 5, Years);
        
        Schedule schedule(startDate, endDate, 3*Months, calendar,
                         Following, Following, DateGeneration::Forward, false);
        
        // Create CDS instrument
        CreditDefaultSwap cds(Protection::Buyer, notional, spread, schedule,
                             Following, Actual360(), true, true, startDate);
        
        // Create ISDA interface
        IsdaCdsInterface isdaInterface;
        
        cout << "\n--- Market Data ---" << endl;
        cout << "Valuation Date: " << today << endl;
        cout << "Risk-free Rate: " << io::rate(riskFreeRate->value()) << endl;
        cout << "Hazard Rate: " << io::rate(hazardRate->value()) << endl;
        cout << "Recovery Rate: " << io::rate(recoveryRate) << endl;
        cout << "CDS Spread: " << io::rate(spread) << endl;
        cout << "Notional: $" << notional << endl;
        cout << "Maturity: " << endDate << endl;
        
        // Price CDS using ISDA methodology
        cout << "\n--- ISDA CDS Pricing ---" << endl;
        
        Real isdaPrice = isdaInterface.priceCds(
            cds, discountCurve, creditCurve, recoveryRate, today);
        
        cout << "ISDA CDS Price: $" << fixed << setprecision(2) << isdaPrice << endl;
        
        // Calculate par spread
        Real parSpread = isdaInterface.calculateParSpread(
            schedule, discountCurve, creditCurve, recoveryRate, today, startDate);
        
        cout << "ISDA Par Spread: " << io::rate(parSpread) << endl;
        
        // Bootstrap default curve from CDS quotes
        cout << "\n--- Curve Bootstrapping ---" << endl;
        
        vector<Real> cdsQuotes = {0.0050, 0.0075, 0.0100, 0.0125, 0.0150};
        vector<Period> tenors = {1*Years, 2*Years, 3*Years, 4*Years, 5*Years};
        
        cout << "Bootstrapping from CDS quotes:" << endl;
        for (size_t i = 0; i < cdsQuotes.size(); ++i) {
            cout << "  " << tenors[i] << ": " << io::rate(cdsQuotes[i]) << endl;
        }
        
        TCurve* bootstrappedCurve = isdaInterface.bootstrapDefaultCurve(
            discountCurve, cdsQuotes, tenors, recoveryRate, today);
        
        if (bootstrappedCurve) {
            cout << "Successfully bootstrapped ISDA spread curve" << endl;
            cout << "Curve has " << bootstrappedCurve->fNumItems << " points" << endl;
            
            // Clean up
            IsdaCdsInterface::freeTCurve(bootstrappedCurve);
        }
        
        // Demonstrate curve conversion
        cout << "\n--- Curve Conversion ---" << endl;
        
        TCurve* isdaDiscountCurve = isdaInterface.convertYieldCurve(
            discountCurve, today);
        
        if (isdaDiscountCurve) {
            cout << "Successfully converted QuantLib yield curve to ISDA format" << endl;
            cout << "ISDA curve has " << isdaDiscountCurve->fNumItems << " points" << endl;
            
            // Clean up
            IsdaCdsInterface::freeTCurve(isdaDiscountCurve);
        }
        
        TCurve* isdaCreditCurve = isdaInterface.convertDefaultCurve(
            creditCurve, today, recoveryRate);
        
        if (isdaCreditCurve) {
            cout << "Successfully converted QuantLib default curve to ISDA format" << endl;
            cout << "ISDA curve has " << isdaCreditCurve->fNumItems << " points" << endl;
            
            // Clean up
            IsdaCdsInterface::freeTCurve(isdaCreditCurve);
        }
        
        cout << "\n=== Example completed successfully ===" << endl;
        
        return 0;
        
    } catch (exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
}