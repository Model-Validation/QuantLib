/*
 * Minimal test for single segment B-spline interpolation
 */

#include <ql/math/interpolations/bsplineinterpolation/bsplineinterpolation.hpp>
#include <ql/math/interpolations/bsplineinterpolation/bsplinestructure.hpp>
#include <ql/math/interpolations/bsplineinterpolation/splinesegment.hpp>
#include <ql/math/interpolations/bsplineinterpolation/splineconstraints.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace QuantLib;
using namespace std;

void testSingleSegment() {
    cout << "========================================" << endl;
    cout << "TEST: Single segment interpolation" << endl;
    cout << "========================================" << endl;

    // Create a single linear segment
    vector<Real> knots = {0.0, 1.0, 2.0};
    Size degree = 1;
    
    cout << "\n1. Creating segment:" << endl;
    auto seg = ext::make_shared<BSplineSegment>(
        knots, degree,
        vector<Integer>(),  // default knot indices
        BSplineSegment::SmoothnessDefault,
        BSplineSegment::TransformDefault,
        BSplineSegment::SideRight,
        1,    // requiredPoints
        true  // isGlobal
    );
    
    Size totalDof = seg->getNumVariables();
    cout << "  Segment DOF = " << totalDof << endl;
    
    cout << "\n2. Creating SplineConstraints (no pre-existing constraints):" << endl;
    
    // Create empty P matrix for regularization
    Eigen::SparseMatrix<Real> P(totalDof, totalDof);
    P.setIdentity();
    P *= 1e-6;  // Small regularization
    
    cout << "  P matrix dimensions: " << P.rows() << "x" << P.cols() << endl;
    
    // No pre-existing constraints - but we need to add a dummy constraint
    // because SplineConstraints expects at least one 
    // OR: Use a different constructor that doesn't require constraints
    vector<Eigen::Triplet<Real>> triplets;
    vector<SplineConstraints::ConstraintType> types;
    vector<Real> b;
    vector<Real> c(totalDof, 0.0);
    
    // Check if we can pass empty constraints
    if (triplets.empty()) {
        cout << "  WARNING: No constraints, adding dummy constraint" << endl;
        // Add a dummy constraint that doesn't affect the solution
        // This constraint says: 0*c[0] + 0*c[1] + 0*c[2] = 0
        // But this might cause issues, so let's try without constraints
    }
    
    cout << "  Total DOF: " << totalDof << endl;
    cout << "  Pre-existing constraints: 0" << endl;
    
    // Create SplineConstraints
    cout << "  Creating SplineConstraints object..." << endl;
    auto constraints = ext::make_shared<SplineConstraints>(
        totalDof,
        P,
        triplets,
        b,
        c,
        types,
        false  // fitData = false for exact interpolation
    );
    cout << "  SplineConstraints created successfully" << endl;
    
    cout << "\n3. Creating BSplineStructure:" << endl;
    vector<ext::shared_ptr<BSplineSegment>> segments = {seg};
    
    cout << "  Number of segments: " << segments.size() << endl;
    
    BSplineStructure structure(
        segments,
        constraints,
        false,  // useSegmentNodes
        false   // rejectZeroNode
    );
    cout << "  BSplineStructure created successfully" << endl;
    
    // Test interpolation with minimal data
    cout << "\n4. Testing interpolation:" << endl;
    vector<Real> x = {0.0, 1.0, 2.0};
    vector<Real> y = {0.0, 0.5, 1.0};
    
    cout << "  Data points:" << endl;
    for (Size i = 0; i < x.size(); ++i) {
        cout << "    (" << x[i] << ", " << y[i] << ")" << endl;
    }
    
    // NOTE: Don't call addInterpolationNodes here - interpolate() does it internally!
    // structure.addInterpolationNodes(x);  // WRONG - causes double addition
    
    cout << "  Calling interpolate..." << endl;
    Eigen::VectorXd solution;
    try {
        solution = structure.interpolate(x, y);
        cout << "  Interpolation completed successfully" << endl;
    } catch (const std::exception& e) {
        cout << "  [ERROR] Interpolation failed: " << e.what() << endl;
        return;
    }
    
    cout << "\n5. Solution coefficients:" << endl;
    for (int i = 0; i < solution.size(); ++i) {
        cout << "  c[" << i << "] = " << solution[i] << endl;
    }
    
    // Evaluate at test points
    cout << "\n6. Evaluating at test points:" << endl;
    for (Size i = 0; i < x.size(); ++i) {
        Real val = seg->value(solution, x[i], 0, BSplineSegment::SideRight);
        cout << "  f(" << x[i] << ") = " << val << " (expected " << y[i] << ")" << endl;
        Real error = std::abs(val - y[i]);
        if (error < 1e-8) {
            cout << "    [PASS] Error = " << error << endl;
        } else {
            cout << "    [FAIL] Error = " << error << endl;
        }
    }
}

int main() {
    cout << "Single Segment B-spline Test" << endl;
    cout << "============================" << endl;
    
    try {
        testSingleSegment();
    } catch (const std::exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        return 1;
    }
    
    cout << "\n============================" << endl;
    cout << "Test completed" << endl;
    
    return 0;
}