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

#ifndef quantlib_test_bspline_staging_hpp
#define quantlib_test_bspline_staging_hpp

#include <boost/test/unit_test.hpp>
#include <ql/math/interpolations/bsplineinterpolation.hpp>

// Forward declaration for future ModeSpan
struct ModeSpan {
    QuantLib::Real start;
    QuantLib::Real end;
    QuantLib::BSplineInterpolation::Mode mode;
};

class BSplineStagingTest {
public:
    // Test current push/pop behavior and timing
    static void testPushPopBehavior();
    
    // Test that mixed modes currently fail
    static void testMixedModeFailure();
    
    // Test staged problem separation (future)
    static void testStagedProblemSeparation();
    
    // Test parameter mapping for both hard and LS modes
    static void testParameterMapping();
    
    // Test mode span resolution strategies
    static void testModeSpanResolution();
    
    // Test Voronoi span generation from points
    static void testVoronoiSpanGeneration();
    
    // Test warm-start performance benefits
    static void testWarmStartPerformance();
    
    static boost::unit_test_framework::test_suite* suite();
};

#endif