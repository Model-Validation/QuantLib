/*
 * Standalone C++ test for debugging multi-segment join constraint enforcement
 * Compile with: cl /I.. /EHsc debug_join_constraints.cpp ../ql/QuantLib-x64-mt-s.lib
 * Or add to CMake as a test target
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

// Helper to print matrix for debugging
void printMatrix(const string& name, const Eigen::SparseMatrix<Real>& mat) {
    cout << name << " (" << mat.rows() << "x" << mat.cols() << "):" << endl;
    for (int i = 0; i < mat.rows(); ++i) {
        cout << "  Row " << i << ": ";
        for (int j = 0; j < mat.cols(); ++j) {
            Real val = mat.coeff(i, j);
            if (std::abs(val) > 1e-10) {
                cout << "[" << j << "]=" << val << " ";
            }
        }
        cout << endl;
    }
}

// Helper to print vector
void printVector(const string& name, const Eigen::VectorXd& vec) {
    cout << name << " (size " << vec.size() << "): ";
    for (int i = 0; i < vec.size(); ++i) {
        if (std::abs(vec[i]) > 1e-10) {
            cout << "[" << i << "]=" << vec[i] << " ";
        }
    }
    cout << endl;
}

void testSimpleJoin() {
    cout << "========================================" << endl;
    cout << "TEST: Simple 2-segment join constraint" << endl;
    cout << "========================================" << endl;

    // Create two linear segments
    vector<Real> knots1 = {0.0, 1.0, 2.0};
    vector<Real> knots2 = {2.0, 3.0, 4.0};
    Size degree = 1;
    
    // Create segments
    cout << "\n1. Creating segments:" << endl;
    auto seg1 = ext::make_shared<BSplineSegment>(
        knots1, degree,
        vector<Integer>(),  // default knot indices
        BSplineSegment::SmoothnessDefault,
        BSplineSegment::TransformDefault,
        BSplineSegment::SideRight,
        1,    // requiredPoints
        true  // isGlobal
    );
    
    auto seg2 = ext::make_shared<BSplineSegment>(
        knots2, degree,
        vector<Integer>(),
        BSplineSegment::SmoothnessDefault,
        BSplineSegment::TransformDefault,
        BSplineSegment::SideRight,
        1,
        true
    );
    
    cout << "  Segment 1: DOF = " << seg1->getNumVariables() << endl;
    cout << "  Segment 2: DOF = " << seg2->getNumVariables() << endl;
    
    // Create join constraint manually
    cout << "\n2. Creating join constraint at x=2.0:" << endl;
    
    // Evaluate basis at join point from each segment
    Eigen::VectorXd left_basis = seg1->evaluateAll(2.0, -1, BSplineSegment::SideRight);
    Eigen::VectorXd right_basis = seg2->evaluateAll(2.0, -1, BSplineSegment::SideLeft);
    
    cout << "  Left basis (from seg1): ";
    for (int i = 0; i < left_basis.size(); ++i) {
        cout << left_basis[i] << " ";
    }
    cout << endl;
    
    cout << "  Right basis (from seg2): ";
    for (int i = 0; i < right_basis.size(); ++i) {
        cout << right_basis[i] << " ";
    }
    cout << endl;
    
    // Build constraint row
    Size totalDof = seg1->getNumVariables() + seg2->getNumVariables();
    vector<Real> constraintRow(totalDof, 0.0);
    
    // Left segment coefficients
    for (Size i = 0; i < seg1->getNumVariables(); ++i) {
        constraintRow[i] = left_basis[i];
    }
    
    // Right segment coefficients (negated for equality)
    Size offset = seg1->getNumVariables();
    for (Size i = 0; i < seg2->getNumVariables(); ++i) {
        constraintRow[offset + i] = -right_basis[i];
    }
    
    cout << "  Constraint row: ";
    for (Real val : constraintRow) {
        cout << val << " ";
    }
    cout << endl;
    
    // Create SplineConstraints with join constraint
    cout << "\n3. Creating SplineConstraints:" << endl;
    
    // Create empty P matrix for now
    Eigen::SparseMatrix<Real> P(totalDof, totalDof);
    P.setIdentity();
    P *= 1e-6;  // Small regularization
    
    cout << "  P matrix dimensions: " << P.rows() << "x" << P.cols() << endl;
    
    // Create A matrix with join constraint
    vector<Eigen::Triplet<Real>> triplets;
    for (Size j = 0; j < totalDof; ++j) {
        if (std::abs(constraintRow[j]) > 1e-10) {
            triplets.push_back(Eigen::Triplet<Real>(0, j, constraintRow[j]));
        }
    }
    
    cout << "  Number of triplets for A matrix: " << triplets.size() << endl;
    
    // Create constraint types (equality)
    vector<SplineConstraints::ConstraintType> types = {
        SplineConstraints::ConstraintType::Equal
    };
    
    // RHS vector (0 for continuity)
    vector<Real> b = {0.0};
    
    // Linear term
    vector<Real> c(totalDof, 0.0);
    
    cout << "  Total DOF: " << totalDof << endl;
    cout << "  Constraints: 1 equality (join)" << endl;
    cout << "  b vector size: " << b.size() << endl;
    cout << "  c vector size: " << c.size() << endl;
    cout << "  types vector size: " << types.size() << endl;
    
    // Create SplineConstraints
    cout << "  Creating SplineConstraints object..." << endl;
    auto constraints = ext::make_shared<SplineConstraints>(
        totalDof,
        P,
        triplets,
        b,
        c,
        types,
        false  // fitData = false for now
    );
    cout << "  SplineConstraints created successfully" << endl;
    
    // DEBUG: Check what SplineConstraints sees
    cout << "\n4. SplineConstraints internal state:" << endl;
    cout << "  Num variables: " << constraints->getNumVariables() << endl;
    cout << "  Num constraints: " << constraints->getNConstraints() << endl;
    
    // Create BSplineStructure
    cout << "\n5. Creating BSplineStructure:" << endl;
    vector<ext::shared_ptr<BSplineSegment>> segments = {seg1, seg2};
    
    cout << "  Number of segments: " << segments.size() << endl;
    cout << "  Constraints pointer valid: " << (constraints != nullptr) << endl;
    
    BSplineStructure structure(
        segments,
        constraints,
        false,  // useSegmentNodes
        false   // rejectZeroNode
    );
    cout << "  BSplineStructure created successfully" << endl;
    
    // Test interpolation with minimal data
    cout << "\n6. Testing interpolation:" << endl;
    vector<Real> x = {0.0, 2.0, 4.0};  // Include join point
    vector<Real> y = {0.0, 1.0, 2.0};
    
    cout << "  Data points:" << endl;
    for (Size i = 0; i < x.size(); ++i) {
        cout << "    (" << x[i] << ", " << y[i] << ")" << endl;
    }
    
    // NOTE: Don't call addInterpolationNodes - interpolate() does it internally!
    // Calling it here would double the parameters and cause dimension mismatch
    
    // Interpolate
    cout << "  Calling interpolate..." << endl;
    Eigen::VectorXd solution;
    try {
        solution = structure.interpolate(x, y);
        cout << "  Interpolation completed successfully" << endl;
    } catch (const std::exception& e) {
        cout << "  [ERROR] Interpolation failed: " << e.what() << endl;
        return;
    }
    
    cout << "\n7. Solution coefficients:" << endl;
    for (int i = 0; i < solution.size(); ++i) {
        cout << "  c[" << i << "] = " << solution[i] << endl;
    }
    
    // Evaluate at join point
    cout << "\n8. Evaluating at join point x=2.0:" << endl;
    
    // From left (seg1)
    Real leftVal = seg1->value(solution.head(seg1->getNumVariables()), 2.0, 0, 
                               BSplineSegment::SideRight);
    
    // From right (seg2)
    Real rightVal = seg2->value(solution.tail(seg2->getNumVariables()), 2.0, 0,
                                BSplineSegment::SideLeft);
    
    cout << "  From left:  " << leftVal << endl;
    cout << "  From right: " << rightVal << endl;
    cout << "  Gap: " << std::abs(rightVal - leftVal) << endl;
    
    if (std::abs(rightVal - leftVal) < 1e-8) {
        cout << "  [PASS] Join is continuous!" << endl;
    } else {
        cout << "  [FAIL] Join is NOT continuous!" << endl;
    }
}

void testLSMode() {
    cout << "\n========================================" << endl;
    cout << "TEST: LS mode with join constraint" << endl;
    cout << "========================================" << endl;
    
    // Similar setup but with more data points for LS
    vector<Real> knots1 = {0.0, 1.0, 2.0};
    vector<Real> knots2 = {2.0, 3.0, 4.0};
    Size degree = 2;  // Quadratic for more DOF
    
    auto seg1 = ext::make_shared<BSplineSegment>(
        knots1, degree,
        vector<Integer>(),
        BSplineSegment::SmoothnessDefault,
        BSplineSegment::TransformDefault,
        BSplineSegment::SideRight,
        1, true
    );
    
    auto seg2 = ext::make_shared<BSplineSegment>(
        knots2, degree,
        vector<Integer>(),
        BSplineSegment::SmoothnessDefault,
        BSplineSegment::TransformDefault,
        BSplineSegment::SideRight,
        1, true
    );
    
    Size totalDof = seg1->getNumVariables() + seg2->getNumVariables();
    cout << "  Total DOF: " << totalDof << " (seg1: " << seg1->getNumVariables() 
         << ", seg2: " << seg2->getNumVariables() << ")" << endl;
    
    // Data for LS fitting
    vector<Real> x = {0.0, 0.5, 1.0, 1.5, 2.5, 3.0, 3.5, 4.0};
    vector<Real> y = {0.0, 0.3, 0.8, 0.9, 1.1, 1.5, 1.7, 2.0};
    
    cout << "  Data points: " << x.size() << " (overdetermined for LS)" << endl;
    
    // Create join constraint
    Eigen::VectorXd left_basis = seg1->evaluateAll(2.0, -1, BSplineSegment::SideRight);
    Eigen::VectorXd right_basis = seg2->evaluateAll(2.0, -1, BSplineSegment::SideLeft);
    
    vector<Eigen::Triplet<Real>> triplets;
    for (Size i = 0; i < seg1->getNumVariables(); ++i) {
        if (std::abs(left_basis[i]) > 1e-10) {
            triplets.push_back(Eigen::Triplet<Real>(0, i, left_basis[i]));
        }
    }
    Size offset = seg1->getNumVariables();
    for (Size i = 0; i < seg2->getNumVariables(); ++i) {
        if (std::abs(right_basis[i]) > 1e-10) {
            triplets.push_back(Eigen::Triplet<Real>(0, offset + i, -right_basis[i]));
        }
    }
    
    // Empty P for now (will be set by LS)
    Eigen::SparseMatrix<Real> P(totalDof, totalDof);
    
    vector<SplineConstraints::ConstraintType> types = {
        SplineConstraints::ConstraintType::Equal
    };
    vector<Real> b = {0.0};
    vector<Real> c(totalDof, 0.0);
    
    // Create with fitData=true for LS mode
    auto constraints = ext::make_shared<SplineConstraints>(
        totalDof, P, triplets, b, c, types,
        true  // fitData = true for LS
    );
    
    cout << "\n  SplineConstraints created with fitData=true" << endl;
    cout << "  Pre-existing constraints: 1 (join)" << endl;
    
    // Create structure and interpolate
    vector<ext::shared_ptr<BSplineSegment>> segments = {seg1, seg2};
    BSplineStructure structure(segments, constraints, false, false);
    
    // This should use LS mode internally
    Eigen::VectorXd solution = structure.interpolate(x, y);
    
    cout << "\n  LS solution obtained" << endl;
    
    // Check join continuity
    Real leftVal = seg1->value(solution.head(seg1->getNumVariables()), 2.0, 0,
                               BSplineSegment::SideRight);
    Real rightVal = seg2->value(solution.tail(seg2->getNumVariables()), 2.0, 0,
                                BSplineSegment::SideLeft);
    
    cout << "\n  Join point evaluation:" << endl;
    cout << "    From left:  " << leftVal << endl;
    cout << "    From right: " << rightVal << endl;
    cout << "    Gap: " << std::abs(rightVal - leftVal) << endl;
    
    if (std::abs(rightVal - leftVal) < 1e-8) {
        cout << "    [PASS] Join is continuous in LS mode!" << endl;
    } else {
        cout << "    [FAIL] Join NOT continuous in LS mode!" << endl;
        cout << "    This is the bug we're trying to fix!" << endl;
    }
}

int main() {
    cout << "Multi-segment Join Constraint Debug Test" << endl;
    cout << "========================================" << endl;
    
    try {
        testSimpleJoin();
        testLSMode();
    } catch (const std::exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        return 1;
    }
    
    cout << "\n========================================" << endl;
    cout << "Debug hints for fixing the issue:" << endl;
    cout << "1. Check if constraints are passed to SCS solver" << endl;
    cout << "2. Check if P matrix combination preserves constraints" << endl;
    cout << "3. Check if SCS respects equality constraints with LS objective" << endl;
    cout << "========================================" << endl;
    
    return 0;
}