/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*!
 Copyright (C) 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2004 Ferdinando Ametrano
 Copyright (C) 2005, 2006, 2007 StatPro Italia srl
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

#if !defined(BOOST_ALL_NO_LIB) && defined(BOOST_MSVC)
#endif
#include "BSplineCurve.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <iomanip>
#include <ios>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/framework/extractor.hpp>
#include <boost/accumulators/framework/accumulator_set.hpp>
#include <ql/compounding.hpp>
#include <ql/currencies/europe.hpp>
#include <ql/handle.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/math/interpolations/bsplineinterpolation/bsplineevaluator.hpp>
#include <ql/math/interpolations/bsplineinterpolation/bsplineinterpolation.hpp>
#include <ql/math/interpolations/bsplineinterpolation/splineconstraints.hpp>
#include <ql/math/interpolations/bsplineinterpolation/splinestructure.hpp>
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


using namespace QuantLib;
using namespace boost::accumulators;


// Define a small epsilon for floating point comparison
const double EPSILON = 1e-6;

// Function to compare vectors
bool compareVectors(const Eigen::VectorXd& vec1,
                    const Eigen::VectorXd& vec2,
                    double epsilon = EPSILON) {
    if (vec1.size() != vec2.size())
        return false;

    for (int i = 0; i < vec1.size(); ++i) {
        if (std::fabs(vec1[i] - vec2[i]) > epsilon)
            return false;
    }

    return true;
}

// Function to compare scalars
bool compareScalars(double val1, double val2, double epsilon = EPSILON) {
    return std::fabs(val1 - val2) < epsilon;
}

// Function to display the progress bar
void printProgressBar(size_t progress, size_t total) {
    int barWidth = 70;
    double progressRatio = static_cast<double>(progress) / static_cast<double>(total);
    int pos = static_cast<int>(barWidth * progressRatio);

    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(progressRatio * 100.0) << " %\r";
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
createRateHelpers2(const ext::shared_ptr<OvernightIndex>& overnightIndex, const Calendar& calendar) {
    // Define the necessary variables
    std::vector<ext::shared_ptr<RateHelper>> rateHelpers;

    // DepositRateHelper
    Rate depositRate = 0.035;
    rateHelpers.push_back(ext::make_shared<DepositRateHelper>(depositRate, overnightIndex));

    // DatedOISRateHelper
    std::vector<std::tuple<Rate, Date, Date>> oisData = {
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
    auto bspline = std::make_unique<BSplineEvaluator>(knots, degree);
    Eigen::VectorXd result, eigenKnots;
    double value;

    // Just for the printing below
    eigenKnots = Eigen::Map<const Eigen::VectorXd>(knots.data(), knots.size());

    std::cout << std::endl << std::string(80, '#') << std::endl;
    std::cout << "Test polynomial of degree " << degree << " with knots:";
    std::cout << std::fixed << std::setprecision(0) << eigenKnots.transpose() << std::endl;
    std::cout << std::string(80, '#') << std::endl << std::endl;

    // Evaluate the B-spline basis functions at xEval1
    result = bspline->evaluateAll(xEval1);

    // Check if the result matches the expected result
    if (compareVectors(result, expectedResult)) {
        std::cout << "B-spline basis values at x = " << std::setprecision(5) << xEval1
                  << " MATCH the expected result.\n";
    } else {
        std::cout << "B-spline basis values at x = " << std::setprecision(5) << xEval1
                  << " do NOT MATCH the expected result.\n";
    }
    std::cout << std::endl
              << std::setprecision(1) << "Difference is: " << std::scientific
              << (expectedResult - result).transpose() << std::endl
              << std::endl;

    // Print the result
    std::cout << std::fixed << std::setprecision(5) << "B-spline basis values at x = " << xEval1
              << ":\n"
              << result.transpose() << std::endl;

    // Value of a spline with given coefficients at another evaluation point
    value = bspline->value(coefficients, xEval2);

    // Check if the computed value matches the expected value
    if (compareScalars(value, expectedValue)) {
        std::cout << "B-spline value at x = " << xEval2 << " MATCHES the expected value.\n";
    } else {
        std::cout << "B-spline value at x = " << xEval2 << " does NOT MATCH the expected value.\n";
    }
    std::cout << "Difference is: " << std::scientific << std::setprecision(1)
              << (expectedValue - value) << std::endl;

    // Print the result
    std::cout << std::fixed << std::setprecision(4) << "B-spline has value at x = " << xEval2
              << ": " << value << std::endl
              << std::endl;
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

    ext::shared_ptr<SplineConstraints> splineConstraints = ext::make_shared<SplineConstraints>(P.size(), P, A, b);

    std::vector<Time> simpleKnots = {0, 1, 2, 3, 4};
    std::vector<Integer> knotIndices = {0, 0, 0, 0, 1, 2, 2, 2, 2, 3, 4, 4, 4, 4};

    size_t degree = 3;
    ext::shared_ptr<BSplineStructure> splineStructure =
        ext::make_shared<BSplineStructure>(simpleKnots, degree, knotIndices, splineConstraints);

    std::vector<Time> xVec = {0.5, 1.5, 2.5, 3.5};
    std::vector<Real> yVec = {2.0, 0.0, 1.0, 4.0};

    Time xEval = 2.7;
    Real expectedValue = 9658418.0 / 3759325.0; // exact rational, double is 2.569189415653076 See BSplineCppWork.nb

    BSplineInterpolation bspline(xVec, yVec, ext::make_shared<BSplineStructure>(*splineStructure));

    std::cout << bspline(2.7) - expectedValue << std::endl;

    std::vector<Time> xVec2 = {0.0, 0.5, 1.5, 2.5, 3.5, 4.0};
    std::vector<Real> yVec2 = {4.0, 3.0, 2.0, 1.0, 3.0, 5.0};

    SpreadedInterpolation spreadedInterpolation(xVec2.begin(), xVec2.end(), yVec2.begin(), Linear(),
                                                ext::make_shared<Interpolation>(bspline));
    
    std::cout << "Spreaded Interpolation Values:" << std::endl;
    for (double x = 0.0; x <= 4.0; x += 0.1) {
        std::cout << x << ",\t" << spreadedInterpolation(x) << std::endl;
    }

    /* ### Setup phase ### */
    Calendar calendar = Sweden();

    // Bootstrap the piecewise linear zero curve
    DayCounter dayCounter = Actual360();

    // Define the SEK-STINATN index
    ext::shared_ptr<OvernightIndex> overnightIndex =
        ext::make_shared<OvernightIndex>("SEK-STINATN", 1, SEKCurrency(), calendar, dayCounter);

    // Create rate helpers
    std::vector<ext::shared_ptr<RateHelper>> rateHelpers =
        createRateHelpers(overnightIndex, calendar);

    Date referenceDate(20, June, 2024);
    Date settlementDate = calendar.advance(referenceDate, 2, Days);

    ///* ### Test PiecewiseYieldCurve<ZeroYield, BSplineModel> implementation ### */
    std::cout << std::endl
              << "### Test PiecewiseYieldCurve<ZeroYield, BSplineModel> implementation ###" << std::endl
              << std::endl;

    try {
        ext::shared_ptr<BSplineZeroCurve> splineCurve = ext::make_shared<BSplineZeroCurve>(
            settlementDate, rateHelpers, dayCounter,
            BSplineModel(ext::make_shared<BSplineStructure>(*splineStructure)));

        Handle<YieldTermStructure> splineYts(splineCurve);

        auto start = std::chrono::high_resolution_clock::now();
        splineYts->zeroRate(settlementDate, dayCounter, Continuous, NoFrequency);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "Execution time: " << duration << " microseconds" << std::endl;

        // Print the rate helper quotes for verification
        RateTimePrintInfo::printInformation(rateHelpers, splineCurve, settlementDate, splineYts,
                                            dayCounter);
        std::vector<Real> yVec3 = {0.04, 0.035, 0.03, 0.025, 0.027, 0.028};

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

        ext::shared_ptr<LinearSpread> spreadedSplineCurve = ext::make_shared<LinearSpread>(
            settlementDate, rateHelpers2, dayCounter,
            SpreadedInterpolationModel(
                Linear(), ext::make_shared<Interpolation>(splineCurve->getInterpolation())
            )
        );

        std::cout << "Spreaded Interpolation Values:" << std::endl;
        for (double x = 0.0; x <= 4.0; x += 0.1) {
            std::cout << std::fixed << std::setprecision(2) << x << ",\t" << std::setprecision(5)
                      << spreadedSplineCurve->zeroRate(x, Continuous, NoFrequency, true).rate()
                      << ",\t" << splineCurve->zeroRate(x, Continuous, NoFrequency, true).rate()
                      << std::endl;
        }

    } catch (std::exception& e) {
        std::cerr << "BSplineZeroCurve error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown BSplineZeroCurve error" << std::endl;
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
//    BSplineStructure splineStructure(simpleKnots, degree);
//}

int main() {
    /* ### Setup phase ### */

    std::cout << std::endl << "### Test B-spline interpolation ###" << std::endl << std::endl;

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
    testInterpolation1();


    return 0;
}
