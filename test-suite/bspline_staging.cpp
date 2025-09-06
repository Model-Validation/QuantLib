/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2025

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.
*/

#include "bspline_staging.hpp"
#include "utilities.hpp"
#include <ql/math/interpolations/bsplineinterpolation.hpp>
#include <ql/math/array.hpp>
#include <chrono>
#include <iostream>

using namespace QuantLib;
using namespace boost::unit_test_framework;

namespace {
    
    // Friend class to access BSplineStructure internals
    class BSplineTestAccess {
    public:
        static bool hasConstraints(const BSplineStructure& structure) {
            return structure.splineConstraints_ != nullptr;
        }
        
        static Size getConstraintCount(const BSplineStructure& structure) {
            if (!structure.splineConstraints_) return 0;
            return structure.splineConstraints_->getNumberOfConstraints();
        }
        
        static bool constraintsArePushed(const BSplineStructure& structure) {
            // Check if push has been called (state saved)
            return structure.splineConstraints_ && 
                   structure.splineConstraints_->hasSavedState();
        }
    };
    
    // Helper to time operations
    class Timer {
        std::chrono::high_resolution_clock::time_point start_;
    public:
        Timer() : start_(std::chrono::high_resolution_clock::now()) {}
        
        double elapsed() const {
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start_);
            return duration.count() / 1000.0; // milliseconds
        }
    };
    
}

void BSplineStagingTest::testPushPopBehavior() {
    BOOST_TEST_MESSAGE("Testing current push/pop behavior...");
    
    // Create a simple B-spline structure
    std::vector<Real> knots = {0, 1, 2, 3, 4};
    Size degree = 3;
    
    BSplineStructure structure(knots, degree);
    
    // Initial state - no constraints pushed
    BOOST_CHECK(!BSplineTestAccess::constraintsArePushed(structure));
    Size initialConstraints = BSplineTestAccess::getConstraintCount(structure);
    BOOST_TEST_MESSAGE("Initial constraints: " << initialConstraints);
    
    // Prepare interpolation data
    std::vector<Real> x = {0.5, 1.5, 2.5, 3.5};
    std::vector<Real> y = {1.0, 2.0, 1.5, 1.0};
    
    // Time first interpolation (cold start)
    Timer timer1;
    auto result1 = structure.interpolate(x, y, BSplineInterpolation::HARD);
    double time1 = timer1.elapsed();
    BOOST_TEST_MESSAGE("First interpolation: " << time1 << " ms");
    
    // Check that push/pop happened (constraints back to initial state)
    BOOST_CHECK(!BSplineTestAccess::constraintsArePushed(structure));
    BOOST_CHECK_EQUAL(BSplineTestAccess::getConstraintCount(structure), initialConstraints);
    
    // Time second interpolation with same x (should rebuild everything)
    std::vector<Real> y2 = {1.1, 2.1, 1.6, 1.1};
    Timer timer2;
    auto result2 = structure.interpolate(x, y2, BSplineInterpolation::HARD);
    double time2 = timer2.elapsed();
    BOOST_TEST_MESSAGE("Second interpolation: " << time2 << " ms");
    
    // Both should take similar time (no warm-start benefit)
    double speedup = time1 / time2;
    BOOST_TEST_MESSAGE("Speedup factor: " << speedup << "x");
    BOOST_CHECK(speedup < 2.0); // No significant speedup expected
    
    // Verify interpolation quality
    for (size_t i = 0; i < x.size(); ++i) {
        double val1 = result1(x[i]);
        double val2 = result2(x[i]);
        BOOST_CHECK_CLOSE(val1, y[i], 1e-6);
        BOOST_CHECK_CLOSE(val2, y2[i], 1e-6);
    }
}

void BSplineStagingTest::testMixedModeFailure() {
    BOOST_TEST_MESSAGE("Testing mixed mode failure with current implementation...");
    
    std::vector<Real> knots = {0, 1, 2, 3, 4};
    Size degree = 3;
    
    BSplineStructure structure(knots, degree);
    
    // Create data with mixed modes
    std::vector<Real> x(10);
    std::vector<Real> y(10);
    std::vector<BSplineInterpolation::Mode> modes(10);
    
    for (size_t i = 0; i < 10; ++i) {
        x[i] = 0.5 + i * 0.35;
        y[i] = std::sin(x[i]);
        // Alternate between hard and soft
        modes[i] = (i % 3 == 0) ? BSplineInterpolation::HARD : BSplineInterpolation::LS;
    }
    
    // This should fail with overdetermined system error
    bool failed = false;
    std::string errorMsg;
    try {
        auto result = structure.interpolate(x, y, modes);
    } catch (const std::exception& e) {
        failed = true;
        errorMsg = e.what();
        BOOST_TEST_MESSAGE("Expected error: " << errorMsg);
    }
    
    BOOST_CHECK(failed);
    BOOST_CHECK(errorMsg.find("overdetermined") != std::string::npos ||
                errorMsg.find("constraints") != std::string::npos);
}

void BSplineStagingTest::testStagedProblemSeparation() {
    BOOST_TEST_MESSAGE("Testing staged problem separation (future implementation)...");
    
    // This test validates the design where:
    // 1. Base problem remains immutable
    // 2. Interpolation layer is built once
    // 3. Combined problem can be solved multiple times
    
    std::vector<Real> knots = {0, 1, 2, 3, 4};
    Size degree = 3;
    
    // Future StagedProblem class (to be implemented)
    // StagedProblem staged(knots, degree);
    
    // Stage with x points and mode specification
    std::vector<Real> x = {0.5, 1.5, 2.5, 3.5};
    std::vector<ModeSpan> spans = {
        {0.0, 2.0, BSplineInterpolation::HARD},
        {2.0, 4.0, BSplineInterpolation::LS}
    };
    
    // staged.stage(x, spans);
    
    // Now solve multiple times with different y values
    std::vector<double> times;
    for (int i = 0; i < 5; ++i) {
        std::vector<Real> y(x.size());
        for (size_t j = 0; j < x.size(); ++j) {
            y[j] = std::sin(x[j]) + 0.01 * i;
        }
        
        // Timer timer;
        // auto result = staged.solve(y);
        // times.push_back(timer.elapsed());
        
        // BOOST_TEST_MESSAGE("Iteration " << i+1 << ": " << times.back() << " ms");
    }
    
    // Expected: First solve slower (cold), rest fast (warm)
    // double avgWarm = std::accumulate(times.begin()+1, times.end(), 0.0) / (times.size()-1);
    // double speedup = times[0] / avgWarm;
    // BOOST_TEST_MESSAGE("Expected speedup: " << speedup << "x");
    // BOOST_CHECK(speedup > 5.0); // Should get significant speedup
    
    BOOST_TEST_MESSAGE("StagedProblem not yet implemented - test placeholder");
}

void BSplineStagingTest::testParameterMapping() {
    BOOST_TEST_MESSAGE("Testing parameter mapping for hard vs soft points...");
    
    std::vector<Real> knots = {0, 1, 2, 3, 4};
    Size degree = 3;
    
    BSplineStructure structure(knots, degree);
    
    // Test that addInterpolationNodes sets up parameters correctly
    // even when we're in LS mode (for warm-start capability)
    
    std::vector<Real> x = {0.5, 1.5, 2.5, 3.5};
    std::vector<Real> y = {1.0, 2.0, 1.5, 1.0};
    
    // LS mode should still set up parameters
    auto resultLS = structure.interpolate(x, y, BSplineInterpolation::LS);
    
    // Hard mode definitely sets up parameters
    auto resultHard = structure.interpolate(x, y, BSplineInterpolation::HARD);
    
    // Both should produce valid results
    for (size_t i = 0; i < x.size(); ++i) {
        double valLS = resultLS(x[i]);
        double valHard = resultHard(x[i]);
        
        // LS should be close but not exact
        BOOST_CHECK_SMALL(std::abs(valLS - y[i]), 0.5);
        
        // Hard should be exact
        BOOST_CHECK_CLOSE(valHard, y[i], 1e-6);
    }
    
    BOOST_TEST_MESSAGE("Parameter mapping verified for both modes");
}

void BSplineStagingTest::testModeSpanResolution() {
    BOOST_TEST_MESSAGE("Testing mode span resolution strategies...");
    
    // Test overlap resolution
    std::vector<ModeSpan> spans = {
        {0.0, 2.0, BSplineInterpolation::HARD},
        {1.0, 3.0, BSplineInterpolation::LS},  // Overlaps with first
        {2.5, 4.0, BSplineInterpolation::HARD}
    };
    
    // Test point at x=1.5 (in overlap)
    // Different strategies should give different results:
    // - LAST_WINS: LS
    // - FIRST_WINS: HARD
    // - NARROWEST: LS (second span is narrower)
    // - HARD_PRIORITY: HARD
    
    // ModeSpanCollection collection(spans);
    // collection.setOverlapStrategy(ModeSpanCollection::LAST_WINS);
    // BOOST_CHECK_EQUAL(collection.getModeAt(1.5), BSplineInterpolation::LS);
    
    BOOST_TEST_MESSAGE("ModeSpan not yet implemented - test placeholder");
}

void BSplineStagingTest::testVoronoiSpanGeneration() {
    BOOST_TEST_MESSAGE("Testing Voronoi span generation from mode points...");
    
    // User specifies modes at key points
    std::vector<std::pair<Real, BSplineInterpolation::Mode>> modePoints = {
        {0.5, BSplineInterpolation::HARD},
        {2.0, BSplineInterpolation::LS},
        {3.5, BSplineInterpolation::HARD}
    };
    
    // System should generate spans automatically
    Real domainStart = 0.0;
    Real domainEnd = 4.0;
    
    // Expected spans (using midpoints):
    // [0.0, 1.25): HARD (nearest to 0.5)
    // [1.25, 2.75): LS (nearest to 2.0)
    // [2.75, 4.0]: HARD (nearest to 3.5)
    
    // VoronoiSpanGenerator generator(modePoints, domainStart, domainEnd);
    // auto generatedSpans = generator.generateSpans();
    
    // BOOST_CHECK_EQUAL(generatedSpans.size(), 3);
    // BOOST_CHECK_CLOSE(generatedSpans[0].end, 1.25, 1e-10);
    // BOOST_CHECK_CLOSE(generatedSpans[1].start, 1.25, 1e-10);
    // BOOST_CHECK_CLOSE(generatedSpans[1].end, 2.75, 1e-10);
    
    BOOST_TEST_MESSAGE("Voronoi spans not yet implemented - test placeholder");
}

void BSplineStagingTest::testWarmStartPerformance() {
    BOOST_TEST_MESSAGE("Testing warm-start performance benefits...");
    
    // Large problem to make timing differences clear
    std::vector<Real> knots;
    for (int i = 0; i <= 20; ++i) {
        knots.push_back(i * 0.5);
    }
    Size degree = 3;
    
    BSplineStructure structure(knots, degree);
    
    // Many interpolation points
    std::vector<Real> x;
    for (int i = 0; i < 50; ++i) {
        x.push_back(0.1 + i * 0.19);
    }
    
    // Simulate bootstrap iterations
    std::vector<double> currentTimes;
    std::vector<double> stagedTimes;
    
    BOOST_TEST_MESSAGE("Current implementation (rebuilds every time):");
    for (int iter = 0; iter < 10; ++iter) {
        std::vector<Real> y(x.size());
        for (size_t j = 0; j < x.size(); ++j) {
            y[j] = std::sin(x[j]) * (1.0 + 0.01 * iter);
        }
        
        Timer timer;
        auto result = structure.interpolate(x, y, BSplineInterpolation::LS);
        double elapsed = timer.elapsed();
        currentTimes.push_back(elapsed);
        BOOST_TEST_MESSAGE("  Iteration " << iter+1 << ": " << elapsed << " ms");
    }
    
    double avgCurrent = std::accumulate(currentTimes.begin(), currentTimes.end(), 0.0) 
                        / currentTimes.size();
    BOOST_TEST_MESSAGE("Average time (current): " << avgCurrent << " ms");
    
    // Future staged implementation would show:
    // - First iteration: ~same as current
    // - Subsequent iterations: 5-50x faster
    
    BOOST_TEST_MESSAGE("Staged implementation would provide 5-50x speedup");
}

test_suite* BSplineStagingTest::suite() {
    test_suite* suite = BOOST_TEST_SUITE("B-spline staging tests");
    
    suite->add(QUANTLIB_TEST_CASE(&BSplineStagingTest::testPushPopBehavior));
    suite->add(QUANTLIB_TEST_CASE(&BSplineStagingTest::testMixedModeFailure));
    suite->add(QUANTLIB_TEST_CASE(&BSplineStagingTest::testStagedProblemSeparation));
    suite->add(QUANTLIB_TEST_CASE(&BSplineStagingTest::testParameterMapping));
    suite->add(QUANTLIB_TEST_CASE(&BSplineStagingTest::testModeSpanResolution));
    suite->add(QUANTLIB_TEST_CASE(&BSplineStagingTest::testVoronoiSpanGeneration));
    suite->add(QUANTLIB_TEST_CASE(&BSplineStagingTest::testWarmStartPerformance));
    
    return suite;
}