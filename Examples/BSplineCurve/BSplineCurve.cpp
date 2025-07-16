/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*!
 Copyright (C) 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2004 Ferdinando Ametrano
 Copyright (C) 2005, 2006, 2007 StatPro Italia srl
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

// ReSharper disable CppClangTidyMiscUseAnonymousNamespace
#if !defined(BOOST_ALL_NO_LIB) && defined(BOOST_MSVC)
#endif

#undef QL_USE_STD_SHARED_PTR

#include "BSplineCurve.h"
#include <chrono>
#include <cmath>
#include <exception>
#include <iomanip>
#include <ios>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <Eigen/Dense>
//#include <boost/accumulators/statistics/variance.hpp>
//#include <boost/accumulators/statistics/stats.hpp>
//#include <boost/accumulators/statistics/mean.hpp>
//#include <boost/accumulators/framework/extractor.hpp>
//#include <boost/accumulators/framework/accumulator_set.hpp>
#include "ql/termstructures/globalbootstrap.hpp"
#include <ql/compounding.hpp>
#include <ql/currencies/europe.hpp>
#include <ql/handle.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/math/interpolations/bsplineinterpolation/bsplineevaluator.hpp>
#include <ql/math/interpolations/bsplineinterpolation/bsplineinterpolation.hpp>
#include <ql/math/interpolations/bsplineinterpolation/splineconstraints.hpp>
#include <ql/math/interpolations/bsplineinterpolation/splinesegment.hpp>
#include <ql/math/interpolations/spreadedinterpolation.hpp>
#include <ql/quote.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/shared_ptr.hpp>
#include <ql/termstructures/yield/bootstraptraits.hpp>
#include <ql/termstructures/yield/oisratehelper.hpp>
#include <ql/termstructures/yield/piecewiseyieldcurve.hpp>
#include <ql/termstructures/yield/ratehelpers.hpp>
#include <ql/termstructures/yield/zerocurve.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/time/calendar.hpp>
#include <ql/time/calendars/sweden.hpp>
#include <ql/time/date.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/frequency.hpp>
#include <ql/time/period.hpp>
#include <ql/time/timeunit.hpp>
#include <ql/types.hpp>
#include <Examples/RateTimeCurve/printinformation.hpp>
#include <ql/math/interpolations/ratetimeinterpolation.hpp>
#include <ql/math/interpolations/bsplineinterpolation/eigenutilities.hpp>
#include <ql/math/interpolations/bsplineinterpolation/bsplinestructure.hpp>

using namespace QuantLib;
//using namespace boost::accumulators;
//  TODO Lets ReSharper relax about the type mismatch between QuantLib::Size and Eigen::Index, should revisit this
//  ReSharper disable CppClangTidyBugproneNarrowingConversions


// Define a small epsilon for floating point comparison
constexpr double global_epsilon = 1e-6;

// Function to compare vectors
bool compareVectors(const Eigen::VectorXd& vec1,
                    const Eigen::VectorXd& vec2,
                    const double epsilon = global_epsilon) {
    if (vec1.size() != vec2.size())
        return false;

    for (int i = 0; i < vec1.size(); ++i) {
        if (std::fabs(vec1[i] - vec2[i]) > epsilon)
            return false;
    }

    return true;
}

// Function to compare scalars
bool compareScalars(double val1, double val2, const double epsilon = global_epsilon) {
    return std::fabs(val1 - val2) < epsilon;
}

// Function to display the progress bar
void printProgressBar(size_t progress, size_t total) {
    constexpr int barWidth = 70;
    const double progressRatio = static_cast<double>(progress) / static_cast<double>(total);
    const int pos = static_cast<int>(barWidth * progressRatio);

    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << static_cast<int>(progressRatio * 100.0) << " %\r";
    std::cout.flush();
}

static std::vector<ext::shared_ptr<RateHelper>>
createRateHelpers(const ext::shared_ptr<OvernightIndex>& overnightIndex, const Calendar& calendar) {
    // Define the necessary variables
    std::vector<ext::shared_ptr<RateHelper>> rateHelpers;

    // DepositRateHelper
    Rate depositRate = 0.03746;
    rateHelpers.push_back(ext::make_shared<DepositRateHelper>(depositRate, overnightIndex));

    // DatedOISRateHelper
    std::vector<std::tuple<Rate, Date, Date>> oisData = {
        //{0.0375, Date(2, July, 2024), Date(3, July, 2024)},
        //{0.0371, Date(20, August, 2024), Date(21, August, 2024)},
        //{0.0354, Date(1, October, 2024), Date(2, October, 2024)},
        //{0.0339, Date(12, November, 2024), Date(13, November, 2024)},
        {0.0325, Date(7, January, 2025), Date(8, January, 2025)},
        //{0.0311, Date(4, February, 2025), Date(5, February, 2025)},
        //{0.0298, Date(1, April, 2025), Date(2, April, 2025)},
        //{0.0286, Date(13, May, 2025), Date(14, May, 2025)},
        {0.0276, Date(1, July, 2025), Date(2, July, 2025)},
        //{0.0266, Date(26, August, 2025), Date(27, August, 2025)},
        //{0.0257, Date(30, September, 2025), Date(1, October, 2025)},
        //{0.0250, Date(11, November, 2025), Date(12, November, 2025)},
        //{0.0244, Date(5, January, 2026), Date(7, January, 2026)},
        //{0.0239, Date(10, February, 2026), Date(11, February, 2026)},
        //{0.0235, Date(7, April, 2026), Date(8, April, 2026)},
        //{0.0232, Date(12, May, 2026), Date(13, May, 2026)}
        };

    for (const auto& ois : oisData) {
        Rate rate = std::get<0>(ois);
        Date startDate = std::get<1>(ois);
        Date endDate = std::get<2>(ois);
        rateHelpers.push_back(ext::make_shared<DatedOISRateHelper>(
            startDate, endDate, Handle<Quote>(ext::make_shared<SimpleQuote>(rate)),
            overnightIndex));
    }

    // Additional OISRateHelpers
    std::vector<std::pair<Period, Rate>> oisRateData = {
        {2 * Years, 0.02870671}, {3 * Years, 0.02661334}, {4 * Years, 0.02546776}};
        //{5 * Years, 0.02477894},  {6 * Years, 0.02451402},  {7 * Years, 0.02444801},
        //{8 * Years, 0.02449672},  {9 * Years, 0.02459572},  {10 * Years, 0.02467563},
        //{12 * Years, 0.02483502}, {15 * Years, 0.02489491}, {20 * Years, 0.02429415},
        //{25 * Years, 0.02333397}, {30 * Years, 0.02235453}};

    for (const auto& oisRate : oisRateData) {
        Period tenor = oisRate.first;
        Rate rate = oisRate.second;
        rateHelpers.push_back(ext::make_shared<OISRateHelper>(
            2, tenor, Handle<Quote>(ext::make_shared<SimpleQuote>(rate)), overnightIndex));
    }

    return rateHelpers;
}

static std::vector<ext::shared_ptr<RateHelper>>
createRateHelpers3(const ext::shared_ptr<OvernightIndex>& overnightIndex, const Calendar& calendar) {
        // Define the necessary variables
        std::vector<ext::shared_ptr<RateHelper>> rateHelpers;

        Handle<YieldTermStructure> discountingCurve = {};

        // Add DepositRateHelper
        Rate depositRate = 0.035; // Example rate for the deposit helper
        rateHelpers.push_back(ext::make_shared<DepositRateHelper>(
            Handle<Quote>(ext::make_shared<SimpleQuote>(depositRate)),
            overnightIndex));

        // Define DatedOISRateHelper data
        struct DatedOISData {
            Rate rate;
            Date startDate;
            Date endDate;
        };

        // Populate with data based on your specifications
        std::vector<DatedOISData> datedOISData = {
            {0.035, Date(1, October, 2024), Date(2, October, 2024)},
            {0.0318, Date(12, November, 2024), Date(13, November, 2024)},
            {0.0284, Date(7, January, 2025), Date(8, January, 2025)},
            {0.0257, Date(4, February, 2025), Date(5, February, 2025)},
            {0.023, Date(25, March, 2025), Date(26, March, 2025)},
            {0.0206, Date(13, May, 2025), Date(14, May, 2025)},
            {0.019, Date(24, June, 2025), Date(25, June, 2025)},
            {0.0178, Date(26, August, 2025), Date(27, August, 2025)},
            {0.017, Date(30, September, 2025), Date(1, October, 2025)},
            {0.0166, Date(11, November, 2025), Date(12, November, 2025)},
            {0.0166, Date(5, January, 2026), Date(7, January, 2026)},
            {0.0167, Date(10, February, 2026), Date(11, February, 2026)},
            {0.0168, Date(7, April, 2026), Date(8, April, 2026)},
            {0.0169, Date(12, May, 2026), Date(13, May, 2026)},
            {0.0171, Date(7, July, 2026), Date(8, July, 2026)}
        };

        // Common settings for DatedOISRateHelpers
        bool telescopicValueDates = false;
        BusinessDayConvention paymentConvention = Following;
        Integer paymentLag = 0;
        Frequency paymentFrequency = Annual;
        Calendar paymentCalendar = Sweden();
        Real overnightSpread = 0.0;
        RateAveraging::Type averagingMethod = RateAveraging::Compound;

        // Add DatedOISRateHelpers to the rateHelpers vector
        for (const auto& data : datedOISData) {
            rateHelpers.push_back(ext::make_shared<DatedOISRateHelper>(
                data.startDate,                                          // Start date
                data.endDate,                                            // End date
                Handle<Quote>(ext::make_shared<SimpleQuote>(data.rate)), // Rate
                overnightIndex,                                          // Overnight index
                discountingCurve,                                        // Discounting curve
                telescopicValueDates,                                    // Telescopic value dates
                averagingMethod,                                         // Averaging method
                paymentLag,                                              // Payment lag
                paymentConvention,                                       // Payment convention
                paymentFrequency,                                        // Payment frequency
                paymentCalendar,                                         // Payment calendar
                overnightSpread                                         // Overnight spread
            )
            );                                       // Averaging method
        }


        // List of tenors and corresponding rates
        std::vector<std::pair<Period, Rate>> oisData = {
            {Period(2, Years), 0.0201805668},  {Period(3, Years), 0.0193503335},
            {Period(4, Years), 0.0192698979},  {Period(5, Years), 0.0194097386},
            {Period(6, Years), 0.019701},      {Period(7, Years), 0.020023},
            {Period(8, Years), 0.020364},      {Period(9, Years), 0.020715},
            {Period(10, Years), 0.0210547023}, {Period(12, Years), 0.021685},
            {Period(15, Years), 0.022175},     {Period(20, Years), 0.021925},
            {Period(25, Years), 0.021105},     {Period(30, Years), 0.020285}
        };


        // Common settings for the OISRateHelper
        paymentLag = 2;
        Period forwardStart = 0 * Days;
        Pillar::Choice pillar = Pillar::MaturityDate;
        auto customPillarDate = Date();
        paymentConvention = ModifiedFollowing;

        // Loop through each pair of tenor and rate to create OISRateHelper instances
        for (const auto& [tenor, rate] : oisData) {
            rateHelpers.push_back(ext::make_shared<OISRateHelper>(
                2,                                                  // Settlement days
                tenor,                                              // Tenor
                Handle<Quote>(ext::make_shared<SimpleQuote>(rate)), // Rate
                overnightIndex,                                     // Overnight index
                discountingCurve,                                   // Discounting curve
                telescopicValueDates,                               // Telescopic value dates
                paymentLag,                                         // Payment lag
                paymentConvention,                                  // Payment convention
                paymentFrequency,                                   // Payment frequency
                paymentCalendar,                                    // Payment calendar
                forwardStart,                                       // Forward start
                overnightSpread,                                    // Overnight spread
                pillar,                                             // Pillar choice
                customPillarDate,                                   // Custom pillar date
                averagingMethod                                     // Averaging method
                ));
        }

        return rateHelpers;
    }

static std::vector<ext::shared_ptr<RateHelper>>
createRateHelpers2(const ext::shared_ptr<OvernightIndex>& overnightIndex, const Calendar& calendar) {
    // Define the necessary variables
    std::vector<ext::shared_ptr<RateHelper>> rateHelpers;

    // DepositRateHelper
    Rate depositRate = 0.035;
    rateHelpers.push_back(ext::make_shared<DepositRateHelper>(depositRate, overnightIndex));

    // DatedOISRateHelper
    const std::vector<std::tuple<Rate, Date, Date>> oisData = {
        {0.035, Date(15, January, 2025), Date(16, January, 2025)},
        {0.030, Date(15, July, 2025), Date(16, July, 2025)},
    };

    for (const auto& ois : oisData) {
        Rate rate = std::get<0>(ois);
        Date startDate = std::get<1>(ois);
        Date endDate = std::get<2>(ois);
        rateHelpers.push_back(ext::make_shared<DatedOISRateHelper>(
            startDate, endDate, Handle<Quote>(ext::make_shared<SimpleQuote>(rate)),
            overnightIndex));
    }

    // Additional OISRateHelpers
    std::vector<std::pair<Period, Rate>> oisRateData = {{2 * Years + 1 * Months, 0.03},
                                                        {3 * Years + 1 * Months, 0.025},
                                                        {4 * Years, 0.025}};

    for (const auto& oisRate : oisRateData) {
        Period tenor = oisRate.first;
        Rate rate = oisRate.second;
        rateHelpers.push_back(ext::make_shared<OISRateHelper>(
            2, tenor, Handle<Quote>(ext::make_shared<SimpleQuote>(rate)), overnightIndex));
    }

    return rateHelpers;
}

void testBSplineEvaluator(const std::vector<double>& knots,
                          size_t degree,
                          const Eigen::VectorXd& expectedResult,
                          double xEval1,
                          double xEval2,
                          const Eigen::VectorXd& coefficients,
                          double expectedValue) {
    const auto bSpline = std::make_unique<BSplineEvaluator>(knots, degree);

    // Just for the printing below
    Eigen::VectorXd eigenKnots = Eigen::Map<const Eigen::VectorXd>(knots.data(), knots.size());

    std::cout << '\n' << std::string(80, '#') << '\n';
    std::cout << "Test polynomial of degree " << degree << " with knots:";
    std::cout << std::fixed << std::setprecision(0) << eigenKnots.transpose() << '\n';
    std::cout << std::string(80, '#') << '\n' << '\n';

    // Evaluate the B-spline basis functions at xEval1
    Eigen::VectorXd result = bSpline->evaluateAll(xEval1);

    // Check if the result matches the expected result
    if (compareVectors(result, expectedResult)) {
        std::cout << "B-spline basis values at x = " << std::setprecision(5) << xEval1
                  << " MATCH the expected result.\n";
    } else {
        std::cout << "B-spline basis values at x = " << std::setprecision(5) << xEval1
                  << " do NOT MATCH the expected result.\n";
    }
    std::cout << '\n'
        << std::setprecision(1) << "Difference is: " << std::scientific
              << (expectedResult - result).transpose() << '\n'
        << '\n';

    // Print the result
    std::cout << std::fixed << std::setprecision(5) << "B-spline basis values at x = " << xEval1
              << ":\n"
              << result.transpose() << '\n';

    // Value of a spline with given coefficients at another evaluation point
    double value = bSpline->value(coefficients, xEval2);

    // Check if the computed value matches the expected value
    if (compareScalars(value, expectedValue)) {
        std::cout << "B-spline value at x = " << xEval2 << " MATCHES the expected value.\n";
    } else {
        std::cout << "B-spline value at x = " << xEval2 << " does NOT MATCH the expected value.\n";
    }
    std::cout << "Difference is: " << std::scientific << std::setprecision(1)
              << (expectedValue - value) << '\n';

    // Print the result
    std::cout << std::fixed << std::setprecision(4) << "B-spline has value at x = " << xEval2
              << ": " << value << '\n'
              << '\n';
}


void testEvaluator1() {
    std::vector<double> knots = {0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4};
    size_t degree = 3;
    Eigen::VectorXd expectedResult(7);
    double xEval1 = 2.3;
    expectedResult << 0.0, 0.0, 343.0 / 6000.0, 3541.0 / 6000.0, 4151.0 / 12000.0, 27.0 / 4000.0,
        0.0;
    double xEval2 = 2.7;
    Eigen::VectorXd coefficients(7);
    coefficients << 0.1, 2.1, 1.6, 0.5, 3.3, 2.4, 1.1;
    double expectedValue = 268837.0 / 120000.0;

    testBSplineEvaluator(knots, degree, expectedResult, xEval1, xEval2, coefficients, expectedValue);
}

void testEvaluator2() {
    std::vector<double> knots = {0, 0, 0, 0, 0, 1, 2, 2, 2, 3, 3, 4, 5, 5, 5, 5, 5};
    size_t degree = 4;

    double xEval1 = 2.3;
    Eigen::VectorXd expectedResult(12);
    expectedResult << 0.0, 0.0, 0.0, 0.0, 2401.0 / 20000.0, 10633.0 / 20000.0, 12177.0 / 40000.0,
        1701.0 / 40000.0, 27.0 / 20000.0, 0.0, 0.0, 0.0;

    double xEval2 = 2.7;
    Eigen::VectorXd coefficients(12);
    coefficients << 0.1, 2.1, 1.6, 0.5, 3.3, 2.4, 1.1, 2.2, 4.5, 2.4, 1.1, 0.5;
    double expectedValue = 2074747.0 / 1200000.0;

    testBSplineEvaluator(knots, degree, expectedResult, xEval1, xEval2, coefficients, expectedValue);
}

/*
void performanceTest() {
    std::vector<double> knots;
    size_t degree;

    std::cout << std::string(80, '#') << std::endl;
    std::cout << "Calculating mean construction and evaluation times:" << std::endl;
    std::cout << std::string(80, '#') << std::endl << std::endl;

    const size_t numPoints = 30;
    degree = 3;
    const size_t numValues = 100'000;
    const size_t numIterations = 1'000;
    const size_t modulus = numIterations / 200;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 30);

    // Random evaluation points
    std::vector<double> values(numValues);
    for (size_t i = 0; i < numValues; ++i) {
        values[i] = dis(gen);
    }

    accumulator_set<double, stats<tag::mean, tag::variance>> constructionAcc;
    accumulator_set<double, stats<tag::mean, tag::variance>> evaluationAcc;

    for (size_t iter = 0; iter < numIterations; ++iter) {
        // Create random knot vectors, first create points, then a valid knot vector
        std::vector<double> points(numPoints);
        for (size_t i = 0; i < numPoints; ++i) {
            points[i] = dis(gen);
        }

        std::sort(points.begin(), points.end());

        knots = std::vector<double>(degree + 1, 0.0);
        knots.insert(knots.end(), points.begin(), points.end());
        knots.insert(knots.end(), degree + 1, 30.0);

        // Random coefficients
        Eigen::VectorXd coeffs = Eigen::VectorXd::Random(knots.size() - degree - 1);
        auto startConstruction = std::chrono::high_resolution_clock::now();
        BSplineEvaluator spline(knots, degree);
        auto endConstruction = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> constructionDuration = endConstruction - startConstruction;
        constructionAcc(constructionDuration.count());

        auto startEvaluation = std::chrono::high_resolution_clock::now();
        for (double x : values) {
            spline.value(coeffs, x);
        }
        auto endEvaluation = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> evaluationDuration = endEvaluation - startEvaluation;

        // Updating evaluation accumulator to record time per value
        evaluationAcc(evaluationDuration.count() / numValues);

        // Update the progress bar
        if ((iter + 1) % modulus == 0)
            printProgressBar(iter + 1, numIterations);
    }
    printProgressBar(numIterations, numIterations);

    std::cout << std::endl;

    double meanConstructionTime = mean(constructionAcc);
    double stdConstructionTime = std::sqrt(variance(constructionAcc));

    double meanEvaluationTime = mean(evaluationAcc);
    double stdEvaluationTime = std::sqrt(variance(evaluationAcc));

    std::cout << std::endl << std::setw(8) << std::fixed << std::setprecision(0);
    std::cout << "Average construction time: " << meanConstructionTime * 1e9 << " ns\n";
    std::cout << "Standard deviation of construction time: " << stdConstructionTime * 1e9
              << " ns\n";
    std::cout << "Average evaluation time per value: " << meanEvaluationTime * 1e9 << " ns\n";
    std::cout << "Standard deviation of evaluation time per value: " << stdEvaluationTime * 1e9
              << " ns\n";
}
*/

/*
spline_structure =
    (SplineStructure([ 0, 1, 2, 3, 4 ], 3, [ 0, 0, 0, 0, 1, 2, 2, 2, 2, 3, 4, 4, 4, 4 ])
         .derivative(4, 0, "equal", "left")
         .value(2, 6, "equal", "average")
         .value(3, 2, "equal")
         .value(4, 0, "equal", side = "left")
         .combine()
         .second_derivative(0, 0, "equal", "right")
         .spline_regularization()
    )
*/

using BSplineZeroCurve=PiecewiseYieldCurve<ZeroYield, BSplineModel>;
using BSplineRateTimeCurve = PiecewiseYieldCurve<RateTime, BSplineModel>;
using PiecewiseLinearZero = PiecewiseYieldCurve<ZeroYield, Linear>;
using MixedRateTimeBSplineBSplineZeroCurve = PiecewiseYieldCurve<ZeroYield, BSplineModel, IterativeBootstrap>;
using LinearSpread = PiecewiseYieldCurve<ZeroYield, SpreadedInterpolationModel<Linear> >;

int testInterpolation1() {
    std::vector<std::vector<Real>> P = {
        { 12.0, -16.5,  3.0,   1.5,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0},
        {-16.5,  24.0, -6.0,  -3.0,   1.5,   0.0,   0.0,  0.0,   0.0,   0.0},
        {  3.0,  -6.0,  6.0,  -6.0,   3.0,   0.0,   0.0,  0.0,   0.0,   0.0},
        {  1.5,  -3.0, -6.0,  24.0, -16.5,   0.0,   0.0,  0.0,   0.0,   0.0},
        {  0.0,   1.5,  3.0, -16.5,  12.0,   0.0,   0.0,  0.0,   0.0,   0.0},
        {  0.0,   0.0,  0.0,   0.0,   0.0,  12.0, -16.5,  3.0,   1.5,   0.0},
        {  0.0,   0.0,  0.0,   0.0,   0.0, -16.5,  24.0, -6.0,  -3.0,   1.5},
        {  0.0,   0.0,  0.0,   0.0,   0.0,   3.0,  -6.0,  6.0,  -6.0,   3.0},
        {  0.0,   0.0,  0.0,   0.0,   0.0,   1.5,  -3.0, -6.0,  24.0, -16.5},
        {  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   1.5,  3.0, -16.5,  12.0}
    };

    std::vector<std::vector<Real>> A = {
         { 0.00,  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, -3.00,  3.00 },
         { 0.00,  0.00, 0.00, 0.00, 0.50, 0.50, 0.00, 0.00,  0.00,  0.00 },
         { 0.00,  0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 0.50,  0.25, -1.00 },
         { 6.00, -9.00, 3.00, 0.00, 0.00, 0.00, 0.00, 0.00,  0.00,  0.00 }
    };

    std::vector<Real> b = {0, 6, 2, 0};
    std::vector<Real> c = {};

    auto splineConstraints = ext::make_shared<SplineConstraints>(P.size(), P, A, b);

    std::vector<Time> simpleKnots = {0, 1, 2, 3, 4};
    std::vector<Integer> knotIndices = {0, 0, 0, 0, 1, 2, 2, 2, 2, 3, 4, 4, 4, 4};

    Size degree = 3;
    auto firstSplineSegment = ext::make_shared<BSplineSegment>(simpleKnots, degree, knotIndices);
    std::vector<ext::shared_ptr<BSplineSegment>> splineSegments = {firstSplineSegment};
    auto splineStructure = ext::make_shared<BSplineStructure>(splineSegments, splineConstraints);

    std::vector xVec = {0.5, 1.5, 2.5, 3.5};
    std::vector yVec = {2.0, 0.0, 1.0, 4.0};

    Time xEval = 2.7;
    Real expectedValue = 9658418.0 / 3759325.0; // exact rational, double is 2.569189415653076 See BSplineCppWork.nb

    BSplineInterpolation bSpline(xVec, yVec, ext::make_shared<BSplineStructure>(*splineStructure));
    //BSplineInterpolation bSpline(xVec, yVec, splineStructure);

    std::cout << bSpline(xEval) - expectedValue << '\n';

    std::vector xVec2 = {0.0, 0.5, 1.5, 2.5, 3.5, 4.0};
    std::vector yVec2 = {4.0, 3.0, 2.0, 1.0, 3.0, 5.0};

    SpreadedInterpolation spreadedInterpolation(xVec2.begin(), xVec2.end(), yVec2.begin(),
        Linear(),
        ext::make_shared<Interpolation>(bSpline));

    std::cout << "Spreaded Interpolation Values:" << '\n';
    for (Size i = 0; i <= 40; i += static_cast<Size>(1)) {
        Real x = static_cast<Real>(i) / 10.0;
        std::cout << x << ",\t" << spreadedInterpolation(x) << '\n';
    }

    /* ### Setup phase ### */
    Calendar calendar = Sweden();

    // Bootstrap the piecewise linear zero curve
    DayCounter dayCounter = Actual360();

    // Define the SEK-STINATN index
    auto overnightIndex =
        ext::make_shared<OvernightIndex>("SEK-STINATN", 1, SEKCurrency(), calendar, dayCounter);

    // Create rate helpers
    std::vector<ext::shared_ptr<RateHelper>> rateHelpers =
        createRateHelpers(overnightIndex, calendar);

    Date referenceDate(20, June, 2024);
    Date settlementDate = calendar.advance(referenceDate, 2, Days);

    ///* ### Test PiecewiseYieldCurve<ZeroYield, BSplineModel> implementation ### */
    std::cout << '\n'
        << "### Test PiecewiseYieldCurve<ZeroYield, BSplineModel> implementation ###" << '\n'
        << '\n';

    try {
        //auto splineCurve = ext::make_shared<BSplineZeroCurve>(
        //    settlementDate, rateHelpers, dayCounter,
        //    BSplineModel(ext::make_shared<BSplineStructure>(*splineStructure))
        //);
        auto splineCurve = ext::make_shared<BSplineZeroCurve>(settlementDate, rateHelpers, dayCounter,
            BSplineModel(ext::make_shared<BSplineStructure>(*splineStructure)));

        Handle<YieldTermStructure> splineYts(splineCurve);

        auto start = std::chrono::high_resolution_clock::now();
        splineYts->update();
        //splineYts->zeroRate(settlementDate, dayCounter, Continuous, NoFrequency);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "Execution time: " << duration << " microseconds" << '\n';

        // Print the rate helper quotes for verification
        RateTimePrintInfo::printInformation(rateHelpers, splineCurve, settlementDate, splineYts,
                                            dayCounter);
        std::vector yVec3 = {0.04, 0.035, 0.03, 0.025, 0.027, 0.028};

        SpreadedInterpolation spreadedCurve(
            xVec2.begin(), xVec2.end(), yVec3.begin(), Linear(),
            ext::make_shared<Interpolation>(splineCurve->getInterpolation())
        );

        //std::cout << "Spreaded Interpolation Values:" << std::endl;
        //for (double x = 0.0; x <= 4.0; x += 0.1) {
        //    std::cout << std::fixed << std::setprecision(2) << x << ",\t" << std::setprecision(5)
        //              << spreadedCurve(x) << ",\t" << splineCurve->getInterpolation().operator()(x)
        //              << std::endl;
        //}

        // Create more rate helpers
        std::vector<ext::shared_ptr<RateHelper>> rateHelpers2 =
            createRateHelpers2(overnightIndex, calendar);

        auto spreadedSplineCurve = ext::make_shared<LinearSpread>(
            settlementDate, rateHelpers2, dayCounter,
            SpreadedInterpolationModel(
                Linear(), ext::make_shared<Interpolation>(splineCurve->getInterpolation())
            )
        );

        std::cout << "Spreaded Interpolation Values:" << '\n';
        for (Size i = 0; i <= 40; ++i) {
            Real x = static_cast<Real>(i) / 10.0;
            std::cout << std::fixed << std::setprecision(2) << x << ",\t" << std::setprecision(5)
                      << spreadedSplineCurve->zeroRate(x, Continuous, NoFrequency, true).rate()
                      << ",\t" << splineCurve->zeroRate(x, Continuous, NoFrequency, true).rate()
                      << '\n';
        }

    } catch (std::exception& e) {
        std::cerr << "BSplineZeroCurve error: " << e.what() << '\n';
        return 1;
    } catch (...) {
        std::cerr << "Unknown BSplineZeroCurve error" << '\n';
        return 1;
    }


    return 0;
}

void createFirstSegment(ext::shared_ptr<BSplineStructure>& splineStructure) {
    // clang-format off
    std::vector<std::vector<Real>> P = {
        {2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0}
    };

    //std::vector<Real> b = {0.03500000, 0.03500820, 0.03505760, 0.03523000, 0.03540230, 0.03563100, 0.03574490, 0.03597130, 0.03613960, 0.03633410, 0.03655330, 0.03668860, 0.03684880, 0.03705850, 0.03718710, 0.03738850, 0.03751150, 0.03781310};

    std::vector<std::vector<Real>> A = {};
    std::vector<Real> b = {};
    std::vector<Real> c = {};
    // clang-format on

    auto splineConstraints = ext::make_shared<SplineConstraints>(P.size(), P, A, b);

    std::vector<Time> simpleKnots = {0.0, 0.00547945, 0.11506849, 0.23013699, 0.38356164,
                                      0.46027397, 0.61369863, 0.72876712, 0.86301370, 1.01643836,
                                      1.11232877, 1.22739726, 1.38082192, 1.47671233, 1.63013699,
                                      1.72602740, 1.87945205, 2.00821918};
    std::vector<Integer> knotIndices = {0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 17};

    Size degree = 1;
    const auto splineSegment = ext::make_shared<BSplineSegment>(
        simpleKnots, degree, knotIndices, BSplineSegment::InterpolationSmoothness::Default,
        BSplineSegment::InterpolationTransform::RateTime);
    std::vector splineSegments = {splineSegment};

    splineStructure = ext::make_shared<BSplineStructure>(splineSegments, splineConstraints);
}


void createSecondStructure(ext::shared_ptr<BSplineStructure>& splineStructure) {
    // clang-format off
    // ReSharper disable once CppInconsistentNaming
   std::vector<std::vector<Real>> P = {
        {2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0}
    };

    // ReSharper disable once CppInconsistentNaming
    std::vector<std::vector<Real>> A = {
        {1.50549451, -3.00824176, 1.00274725, 0.99726776, -0.24829375, -0.24897401, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {-0.50274725, 0.0, 1.50274725, -1.49726776, 0.24829375, 0.24897401, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, -0.24965800, -0.24829375, 1.49658377, -1.49931601, 0.25034200, 0.25034200, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, -0.24965800, -0.25034200, 1.50000000, -1.50000000, 0.25068306, 0.24931694, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25068306, -0.25068306, 1.49863388, -1.49318429, 0.24694065, 0.24897588, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.24829375, -0.24694065, 1.49386641, -1.50206326, 0.25206091, 0.25137033, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25000047, -0.25206091, 1.50343312, -1.50068587, 0.24965706, 0.24965706, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25034294, -0.24965706, 1.50000000, -1.50000000, 0.33378871, 0.16621129, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.33378871, -0.33378871, 1.33242259, -0.83037620, 0.09908292, 0.06644811, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.19898022, -0.09908292, 0.69784409, -0.53383375, 0.08387482, 0.05017798, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.12506846, -0.08387482, 0.45894328, -0.32478106, 0.03738030, 0.03740077, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.06248289, -0.03738030, 0.29969894, -0.29980835, 0.04999998, 0.04997261, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.04997261, -0.04999998, 0.29980835, -0.29969894, 0.0, 0.09986319},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04997261, 0.04999998, -0.19989047, -0.19978106, 0.59934319, -0.29964425}
    };
    std::vector<Real> b = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    std::vector<Real> c = {};
    // clang-format on

    auto splineConstraints =
        ext::make_shared<SplineConstraints>(P.size(), P, A, b);

    std::vector<Time> simpleKnots = {
        2.00821918, 3.00547945,  4.00821918,  5.00821918,  6.00821918,  7.01369863,  8.01095890,
        9.01095890, 10.01095890, 12.01917808, 15.01369863, 20.01917808, 25.02191781, 30.02739726};
    std::vector<Integer> knotIndices = {0, 0, 0, 0, 1, 1, 2,  2,  3,  3,  4,  4,  5,  5,  6,  6,
                                        7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 13, 13};

    size_t degree = 3;
    const auto splineSegment = ext::make_shared<BSplineSegment>(simpleKnots, degree, knotIndices);
    std::vector splineSegments = {splineSegment};

    splineStructure = ext::make_shared<BSplineStructure>(splineSegments, splineConstraints);
}

void createCompositeStructure(ext::shared_ptr<BSplineStructure>& compositeSplineStructure) {
    constexpr Size nVariablesFirst = 18;
    constexpr Size nVariablesSecond = 28;
    constexpr Size nVariables = nVariablesFirst + nVariablesSecond;
    constexpr Size nConstraintsFirst = 0;
    constexpr Size nConstraintsSecond = 14;
    constexpr Size nConstraintsJoint = 1;
    constexpr Size nConstraints = nConstraintsFirst + nConstraintsSecond + nConstraintsJoint;

    Eigen::SparseMatrix<double> P(nVariables, nVariables);
    Eigen::SparseMatrix<double> P_half(nVariablesFirst - 1, nVariablesFirst);
    Eigen::MatrixXd P_first(nVariablesFirst, nVariablesFirst);
    Eigen::MatrixXd P_second(nVariablesSecond, nVariablesSecond);

    // Create a larger matrix to hold the direct sum
    Eigen::MatrixXd P_directSum = Eigen::MatrixXd::Zero(nVariables, nVariables);

    Eigen::MatrixXd A_first(nConstraintsFirst, nVariablesFirst),
        A_second(nConstraintsSecond, nVariablesSecond), A_joint(nConstraintsJoint, nVariables);

    // clang-format off

    std::vector<std::vector<Time>> simpleKnots = {
        {0.0, 0.01388889, 0.03611111, 0.15277778, 0.30833333, 0.38611111, 0.52222222, 0.65833333, 0.77500000, 0.95000000, 1.04722222, 1.16388889, 1.31944444, 1.41666667, 1.57222222, 1.66944444, 1.82500000, 2.03888889},
        {2.03888889, 3.05277778, 4.07500000, 5.08611111, 6.09722222, 7.11111111, 8.12777778, 9.14166667, 10.16111111, 12.18611111, 15.22777778, 20.30277778, 25.37500000, 30.44722223}
    };

    std::vector<std::vector<Integer>> knotIndices = {
        {0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 17},
        {0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 13, 13}
    };

    P_first <<
        144.00000000, -144.00000000, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        -144.00000000, 234.00000000, -90.00000000, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, -90.00000000, 107.14285714, -17.14285714, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, -17.14285714, 30.00000000, -12.85714286, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, -12.85714286, 38.57142857, -25.71428571, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, -25.71428571, 40.40816327, -14.69387755, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, -14.69387755, 29.38775510, -14.69387755, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -14.69387755, 31.83673469, -17.14285714, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -17.14285714, 28.57142857, -11.42857143, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -11.42857143, 32.00000000, -20.57142857, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -20.57142857, 37.71428571, -17.14285714, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -17.14285714, 30.00000000, -12.85714286, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12.85714286, 33.42857143, -20.57142857, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -20.57142857, 33.42857143, -12.85714286, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12.85714286, 33.42857143, -20.57142857, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -20.57142857, 33.42857143, -12.85714286, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -12.85714286, 22.20779221, -9.35064935,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -9.35064935, 9.35064935;

    P_second <<
        23.02714791, -34.54072187, 5.78034818, 5.73322578, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        -34.54072187, 69.08144374, -34.54072187, 0.00000000, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        5.78034818, -34.54072187, 45.86657451, -22.74633111, 2.80465495, 2.83547534, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        5.73322578, 0.00000000, -22.74633111, 45.12179813, -30.92105284, 2.81236005, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 2.80465495, -30.92105284, 45.18403229, -22.84028006, 2.88632283, 2.88632283, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 2.83547534, 2.81236005, -22.84028006, 46.18254430, -31.90814028, 2.91804066, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 2.88632283, -31.90814028, 46.43490793, -23.21745396, 2.90616279, 2.89820070, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 2.88632283, 2.91804066, -23.21745396, 46.43490793, -31.92001815, 2.89820070, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.90616279, -31.92001815, 46.37129846, -23.12212690, 2.88628491, 2.87839889, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.89820070, 2.89820070, -23.12212690, 46.11755722, -31.66234459, 2.87051287, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.88628491, -31.66234459, 45.99138012, -22.93286031, 2.85485916, 2.86268069, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.87839889, 2.87051287, -22.93286031, 45.74040444, -31.41131505, 2.85485916, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.85485916, -31.41131505, 45.74040444, -22.93286031, 2.88230958, 2.86660217, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.86268069, 2.85485916, -22.93286031, 45.99138012, -31.65051555, 2.87445588, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.88230958, -31.65051555, 45.92880728, -22.83925703, 3.77713505, 1.90152066, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.86660217, 2.87445588, -22.83925703, 45.42958482, -30.22254400, 1.89115815, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.77713505, -30.22254400, 34.05677947, -8.57259126, 0.57704864, 0.38417211, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.90152066, 1.89115815, -8.57259126, 8.63138819, -4.04487925, 0.19340352, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.57704864, -4.04487925, 4.81445214, -1.60261900, 0.16006413, 0.09593334, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.38417211, 0.19340352, -1.60261900, 2.13389818, -1.17272275, 0.06386795, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.16006413, -1.17272275, 1.36401156, -0.40875552, 0.02869344, 0.02870915, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.09593334, 0.06386795, -0.40875552, 0.48996967, -0.25822208, 0.01720663, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.02869344, -0.25822208, 0.36732628, -0.18376372, 0.02298304, 0.02298304, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.02870915, 0.01720663, -0.18376372, 0.36772872, -0.25287641, 0.02299563, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.02298304, -0.25287641, 0.36782938, -0.18391469, 0.00000000, 0.04597867,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.02298304, 0.02299563, -0.18391469, 0.36782938, -0.27587204, 0.04597867,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00000000, -0.27587204, 0.55174408, -0.27587204,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04597867, 0.04597867, -0.27587204, 0.18391469;

    A_second <<
        1.48147041, -2.95890411, 0.98630137, 0.97826087, -0.24223332, -0.24489522, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        -0.49516904, 0.0, 1.47743370, -1.46939320, 0.24223332, 0.24489522, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, -0.24422427, -0.24223332, 1.47006415, -1.48081427, 0.24860386, 0.24860386, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, -0.24590164, -0.24860386, 1.48351648, -1.48351648, 0.24759191, 0.24691358, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.24759191, -0.24759191, 1.48283815, -1.48012853, 0.24657395, 0.24590026, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.24725136, -0.24657395, 1.47877743, -1.47608262, 0.24522932, 0.24590118, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.24590118, -0.24522932, 1.47608262, -1.47877743, 0.24758728, 0.24623803, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.24691081, -0.24758728, 1.47810465, -1.47272971, 0.32533830, 0.16378485, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.32712098, -0.32533830, 1.30939358, -0.82229431, 0.09927040, 0.06608961, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.19718834, -0.09927040, 0.69119558, -0.52613554, 0.08215792, 0.04924079, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.12340592, -0.08215792, 0.45197041, -0.32024762, 0.03691042, 0.03693063, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.06158478, -0.03691042, 0.29559347, -0.29570138, 0.04930155, 0.04930155, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.04927457, -0.04930155, 0.29572837, -0.29572837, 0.0, 0.09857612,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04927457, 0.04930155, -0.19715225, -0.19715225, 0.59145674, -0.29572837;

    A_joint <<
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.49046322, 1.00000000, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    // clang-format on

    //P_half.reserve(2*(nVariablesFirst-1));
    //for (int i = 0; i < static_cast<int>(nVariablesFirst)-1; ++i) {
    //    P_half.insert(i, i) = simpleKnots[0][i + 1];     // Diagonal elements (tt[1] to tt[17])
    //    P_half.insert(i, i + 1) = -simpleKnots[0][i];    // Superdiagonal elements (-tt[0] to -tt[16])
    //}

    ////P.setIdentity();
    ////P = 2.0 * P;
    //P_first = 1.0 * (P_half.transpose() * P_half).toDense();
    P_directSum.topLeftCorner(P_first.rows(), P_first.cols()) = P_first;
    P_directSum.bottomRightCorner(P_second.rows(), P_second.cols()) = 2.0 * P_second;
    P = P_directSum.sparseView();


    // Create a larger matrix to hold the direct sum
    Eigen::MatrixXd A_directSum = Eigen::MatrixXd::Zero(nConstraintsFirst + nConstraintsSecond,
                                                      nVariablesFirst + nVariablesSecond);

    // Place matrix A in the top-left block
    //directSum.topLeftCorner(A_first.rows(), A_first.cols()) = A_first;
    A_directSum.bottomRightCorner(A_second.rows(), A_second.cols()) = A_second;

    Eigen::MatrixXd A(nConstraints, nVariables);

    // Vertical stacking of the constraints from both parts, and the constraint ensuring continuity across the boundary
    A << A_directSum,
         A_joint;

    std::vector<Eigen::Triplet<double>> A_triplets = EigenUtilities::convertToTriplets(A.sparseView());

    // clang-format off
    std::vector<Real> b = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    // clang-format on

    assert(b.size() == nConstraints);

    auto splineConstraints = ext::make_shared<SplineConstraints>(nVariables, P, A_triplets, b);

    std::vector<Size> degrees = {1, 3};

    using InterpolationTransform = BSplineSegment::InterpolationTransform;
    std::vector transforms = {InterpolationTransform::RateTime, InterpolationTransform::Default};

    std::vector<ext::shared_ptr<BSplineSegment>> splineSegments;

    for (Size i = 0; i < simpleKnots.size(); ++i) {
        splineSegments.push_back(ext::make_shared<BSplineSegment>(
            simpleKnots[i], static_cast<Integer>(degrees[i]), knotIndices[i],
            BSplineSegment::InterpolationSmoothness::Default,
            transforms[i]));
    }

    compositeSplineStructure = ext::make_shared<BSplineStructure>(splineSegments, splineConstraints);
}

int testInterpolation2() {
    ext::shared_ptr<BSplineStructure> splineStructureFirst, splineStructureSecond,
        splineStructureBoth;

    createFirstSegment(splineStructureFirst);

    std::vector<Time> xVecFirst = {0.0,        0.00547945, 0.11506849, 0.23013699, 0.38356164,
                              0.46027397, 0.61369863, 0.72876712, 0.86301370, 1.01643836,
                              1.11232877, 1.22739726, 1.38082192, 1.47671233, 1.63013699,
                              1.72602740, 1.87945205, 2.00821918};
    std::vector<Real> yVecFirst = {0.035,     0.0350082, 0.0350576, 0.03523,   0.0354023, 0.035631,
                              0.0357449, 0.0359713, 0.0361396, 0.0363341, 0.0365533, 0.0366886,
                              0.0368488, 0.0370585, 0.0371871, 0.0373885, 0.0375115, 0.0378131};

    // Time xEval = 2.7;
    // Real expectedValue =
    //     9658418.0 / 3759325.0; // exact rational, double is 2.569189415653076 See
    //     BSplineCppWork.nb

    //const BSplineInterpolation bSplineFirst(
    //    xVecFirst, yVecFirst, ext::make_shared<BSplineStructure>(*splineStructureFirst));
    const BSplineInterpolation bSplineFirst(xVecFirst, yVecFirst, splineStructureFirst);

    std::cout << bSplineFirst(1.7) << '\n';

    createSecondStructure(splineStructureSecond);

    std::vector xVecSecond = {2.00821918,  3.00547945,  4.00821918,  5.00821918,  6.00821918,
                              7.01369863,  8.01095890,  9.01095890,  10.01095890, 12.01917808,
                              15.01369863, 20.01917808, 25.02191781, 30.02739726};
    std::vector yVecSecond = {0.03497620, 0.03492030, 0.03481140, 0.03463480, 0.03437550,
                              0.03401640, 0.03354720, 0.03294820, 0.03220500, 0.03021500,
                              0.02610580, 0.02036040, 0.02226920, 0.02499780};

    //Time xEval = 2.7;
    //Real expectedValue =
    //    9658418.0 / 3759325.0; // exact rational, double is 2.569189415653076 See BSplineCppWork.nb

    //const BSplineInterpolation bSplineSecond(xVec, yVec, ext::make_shared<BSplineStructure>(*splineStructureSecond));
    const BSplineInterpolation bSplineSecond(xVecSecond, yVecSecond, splineStructureSecond);

    std::cout << bSplineSecond(2.7) << '\n';

    //std::vector<Time> xVec2 = {0.0, 0.5, 1.5, 2.5, 3.5, 4.0};
    //std::vector<Real> yVec2 = {4.0, 3.0, 2.0, 1.0, 3.0, 5.0};

    //SpreadedInterpolation spreadedInterpolation(xVec2.begin(), xVec2.end(), yVec2.begin(), Linear(),
    //                                            ext::make_shared<Interpolation>(bspline));

    //std::cout << "Spreaded Interpolation Values:" << std::endl;
    //for (double x = 0.0; x <= 4.0; x += 0.1) {
    //    std::cout << x << ",\t" << spreadedInterpolation(x) << std::endl;
    //}
    createCompositeStructure(splineStructureBoth);

    std::vector xVecBoth(xVecFirst);
    xVecBoth.insert(xVecBoth.end(), xVecSecond.begin(), xVecSecond.end());

    std::vector yVecBoth(yVecFirst);
    xVecBoth.insert(yVecBoth.end(), yVecSecond.begin(), yVecSecond.end());

    const BSplineInterpolation bSplineBoth(xVecBoth, xVecSecond, splineStructureBoth);
    std::cout << bSplineBoth(1.7)  << "," << bSplineBoth(2.7) << '\n';

    /* ### Setup phase ### */
    Calendar calendar = Sweden();

    // Bootstrap the piecewise linear zero curve
    DayCounter dayCounter = Actual360();

    // Define the SEK-STINATN index
    const auto overnightIndex =
        ext::make_shared<OvernightIndex>("SEK-STINATN", 1, SEKCurrency(), calendar, dayCounter);

    // Create rate helpers
    std::vector<ext::shared_ptr<RateHelper>> rateHelpers =
        createRateHelpers3(overnightIndex, calendar);

    const Date referenceDate(21, August, 2024);
    Date settlementDate = calendar.advance(referenceDate, 2, Days);

    ///* ### Test PiecewiseYieldCurve<ZeroYield, BSplineModel> implementation ### */
    std::cout << '\n'
        << "### Test PiecewiseYieldCurve<ZeroYield, BSplineModel> implementation ###"
              << '\n'
        << '\n';

    try {
        BSplineModel interpolatorFactory(splineStructureBoth);
        // MixedRateTimeBSplineBSpline interpolator_factory(ext::make_shared<BSplineSegment>(*splineStructureFirst),
        //                            ext::make_shared<BSplineSegment>(*splineStructureSecond), 18,
        //                            MixedInterpolation::Behavior::SplitRanges);

        auto splineCurve =
            ext::make_shared<MixedRateTimeBSplineBSplineZeroCurve>(settlementDate, rateHelpers, dayCounter, interpolatorFactory);

        //ext::shared_ptr<BSplineZeroCurve> splineCurve = ext::make_shared<BSplineZeroCurve>(
        //    settlementDate, rateHelpers, dayCounter,
        //    BSplineModel(ext::make_shared<BSplineSegment>(*compositeSplineStructure)));


        Handle<YieldTermStructure> splineYts(splineCurve);

        const auto start = std::chrono::high_resolution_clock::now();
        std::cout << splineYts->zeroRate(settlementDate, dayCounter, Continuous, NoFrequency).rate();
        const auto end = std::chrono::high_resolution_clock::now();
        const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "Execution time: " << duration << " microseconds" << '\n';

        // Print the rate helper quotes for verification
        RateTimePrintInfo::printInformation(rateHelpers, splineCurve, settlementDate, splineYts,
                                            dayCounter);

        std::cout << '\n' << "Interpolation Values:" << '\n';
        for (int x = 1; x <= 20; x += 1) {
            Real value = static_cast<Real>(x) / 2.0;
            std::cout << std::fixed << std::setprecision(2) << value << ",\t" << std::setprecision(6)
                      << splineCurve->zeroRate(value, Continuous, NoFrequency).rate()
                      << '\n';
        }

    } catch (std::exception& e) {
        std::cerr << "BSplineZeroCurve error: " << e.what() << '\n';
        return 1;
    } catch (...) {
        std::cerr << "Unknown BSplineZeroCurve error" << '\n';
        return 1;
    }


    return 0;
}

int testLinearInterpolation() {
    /* ### Setup phase ### */
    Calendar calendar = Sweden();

    // Bootstrap the piecewise linear zero curve
    DayCounter dayCounter = Actual360();

    // Define the SEK-STINATN index
    const auto overnightIndex =
        ext::make_shared<OvernightIndex>("SEK-STINATN", 1, SEKCurrency(), calendar, dayCounter);

    // Create rate helpers
    std::vector<ext::shared_ptr<RateHelper>> rateHelpers =
        createRateHelpers3(overnightIndex, calendar);

    const Date referenceDate(19, September, 2024);
    Date settlementDate = calendar.advance(referenceDate, 2, Days);

    ///* ### Test PiecewiseYieldCurve<ZeroYield, Linear> implementation ### */
    std::cout << '\n'
              << "### Test PiecewiseYieldCurve<ZeroYield, Linear> implementation ###" << '\n'
              << '\n';

    try {
        auto splineCurve =
            ext::make_shared<PiecewiseLinearZero>(referenceDate, rateHelpers, dayCounter);

        Handle<YieldTermStructure> splineYts(splineCurve);

        const auto start = std::chrono::high_resolution_clock::now();
        std::cout << splineYts->zeroRate(settlementDate, dayCounter, Continuous, NoFrequency).rate();
        const auto end = std::chrono::high_resolution_clock::now();
        const auto duration =
            std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "Execution time: " << duration << " microseconds" << '\n';

        // Print the rate helper quotes for verification
        RateTimePrintInfo::printInformation(rateHelpers, splineCurve, settlementDate, splineYts,
                                            dayCounter);

        std::cout << '\n' << "Interpolation Values:" << '\n';
        for (int x = 1; x <= 20; x += 1) {
            Real value = static_cast<Real>(x) / 2.0;
            std::cout << std::fixed << std::setprecision(2) << value << ",\t"
                      << std::setprecision(6)
                      << splineCurve->zeroRate(value, Continuous, NoFrequency).rate() << '\n';
        }

    } catch (std::exception& e) {
        std::cerr << "PiecewiseLinearZero error: " << e.what() << '\n';
        return 1;
    } catch (...) {
        std::cerr << "Unknown PiecewiseLinearZero error" << '\n';
        return 1;
    }


    return 0;
}

int testCompositeInterpolation1() {
    ext::shared_ptr<BSplineStructure> compositeSplineStructure;

    createCompositeStructure(compositeSplineStructure);

    // ### Setup phase ###
    Calendar calendar = Sweden();

    // Bootstrap the piecewise linear zero curve
    DayCounter dayCounter = Actual360();

    // Define the SEK-STINATN index
    const auto overnightIndex =
        ext::make_shared<OvernightIndex>("SEK-STINATN", 2, SEKCurrency(), calendar, dayCounter);

    // Create rate helpers
    std::vector<ext::shared_ptr<RateHelper>> rateHelpers =
        createRateHelpers3(overnightIndex, calendar);

    const Date referenceDate(19, September, 2024);
    const Date otherDate(13, November, 2024);
    Date settlementDate = calendar.advance(referenceDate, 2, Days);

    // ### Test PiecewiseYieldCurve<ZeroYield, BSplineModel> implementation ###
    std::cout << '\n'
        << "### Test PiecewiseYieldCurve<ZeroYield, BSplineModel> implementation ###"
              << '\n'
        << '\n';

    try {
        BSplineModel interpolatorFactory(compositeSplineStructure);

        auto splineCurve = ext::make_shared<MixedRateTimeBSplineBSplineZeroCurve>(
            referenceDate, rateHelpers, dayCounter, interpolatorFactory);

        Handle<YieldTermStructure> splineYts(splineCurve);

        const auto start = std::chrono::high_resolution_clock::now();
        // ReSharper disable once CppExpressionWithoutSideEffects
        splineYts->zeroRate(settlementDate, dayCounter, Continuous, NoFrequency).rate();
        const auto end = std::chrono::high_resolution_clock::now();
        const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "Execution time: " << duration << " milliseconds" << '\n';
        std::cout << splineYts->zeroRate(otherDate, dayCounter, Continuous, NoFrequency).rate() << std::endl;

        // Print the rate helper quotes for verification
        RateTimePrintInfo::printInformation(rateHelpers, splineCurve, settlementDate, splineYts,
                                            dayCounter);

        std::cout << '\n' << "Interpolation Values:" << '\n';
        for (int x = 1; x <= 20; x += 1) {
            Real value = static_cast<Real>(x) / 2.0;
            std::cout << std::fixed << std::setprecision(2) << value << ",\t"
                      << std::setprecision(6)
                      << splineCurve->zeroRate(value, Continuous, NoFrequency).rate() << '\n';
        }

    } catch (std::exception& e) {
        std::cerr << "BSplineZeroCurve error: " << e.what() << '\n';
        return 1;
    } catch (...) {
        std::cerr << "Unknown BSplineZeroCurve error" << '\n';
        return 1;
    }


    return 0;
}

//static void testBSplineInterpolation(const std::vector<double>& simpleKnots,
//                                     size_t degree,
//                                     std::vector<double> xVec,
//                                     std::vector<double> yVec,
//                                     double xEval,
//                                     double expectedValue) {
//    BSplineSegment compositeSplineStructure(simpleKnots, degree);
//}

int main() {
    const Calendar calendar = Sweden();
    Date todayDate(19, September, 2024);
    // must be a business day
    todayDate = calendar.adjust(todayDate);

    Settings::instance().evaluationDate() = todayDate;

    /* ### Setup phase ### */

    std::cout << '\n' << "### Test B-spline interpolation ###" << '\n' << '\n';

/*  ##################
        ##### Test 1 #####
        ################## */
    // testEvaluator1();

    /*  ##################
        ##### Test 2 #####
        ################## */
    // testEvaluator2();

    /*  ############################
        ##### Performance test #####
        ############################ */
    // performanceTest();

    /*  ############################
        ####### BSpline test #######
        ############################ */
    //testInterpolation1();
    //testInterpolation2();
    //testLinearInterpolation();
    testCompositeInterpolation1();

    //std::vector<Real> knots_ = {0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0};
    //std::vector<Real> coeffs = {0.23, -0.67, 0.44, -0.81, 0.16};
    //std::vector<Real> xx = {-1.0, -0.5, 0.0, 0.1, 0.5, 1.5, 2.5, 3.5, 3.9, 4.0, 4.5};

    //Size i;
    //Size degree_ = 1;

    //for (const auto& x : xx) {
    //    auto it = std::upper_bound(knots_.begin(), knots_.end(), x); // Should be binary search
    //    if (it != knots_.end()) {
    //        i = std::distance(knots_.begin(), it);
    //    } else if (x == knots_.back()) {
    //        std::cout << "x is the last knot" << std::endl;
    //        i = knots_.size() - degree_ - 1; // TODO Not -2, right?
    //        std::cout << std::distance(knots_.begin(), it)  - 1 << std::endl;
    //    } else {
    //        std::cout << "x is outside the range of the knots" << std::endl;
    //    }
    //    std::cout << "(x,i) is " << x << "," << i << std::endl;
    //}

    //BSplineEvaluator bspline(knots_, 1);
    //Eigen::VectorXd eigen_coeffs(Eigen::Map<Eigen::VectorXd>(coeffs.data(), coeffs.size()));
    //for (const auto& x : xx) {
    //    std::cout << bspline.value(eigen_coeffs, x) << std::endl;
    //}

    return 0;
}
