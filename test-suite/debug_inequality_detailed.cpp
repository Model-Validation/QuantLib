/*
 * Detailed test to debug inequality constraint violations in hard mode
 * This test adds detailed output to track where constraints are lost
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

void printVector(const string& name, const vector<Real>& vec) {
    cout << name << " [" << vec.size() << "]: ";
    for (Real val : vec) {
        cout << val << " ";
    }
    cout << endl;
}

void testInequalityConstraint(bool useLS) {
    cout << "\n========================================" << endl;
    cout << "TEST: Inequality constraint y(0.5) >= 0.7" << endl;
    cout << "Mode: " << (useLS ? "LS (fitData=true)" : "HARD (fitData=false)") << endl;
    cout << "========================================" << endl;
    
    // Simple linear B-spline from 0 to 1
    vector<Real> knots = {0.0, 0.0, 1.0, 1.0}; // Linear B-spline knots
    Size degree = 1;
    
    // Create segment
    auto segment = ext::make_shared<BSplineSegment>(
        knots, degree,
        vector<Integer>(),
        BSplineSegment::SmoothnessDefault,
        BSplineSegment::TransformDefault,
        BSplineSegment::SideRight,
        1,    // requiredPoints
        true  // isGlobal
    );
    
    Size numVars = segment->getNumVariables();
    cout << "\n1. Segment setup:" << endl;
    cout << "   Knots: ";
    printVector("", knots);
    cout << "   Degree: " << degree << endl;
    cout << "   Variables: " << numVars << endl;
    
    // Evaluate basis at constraint point
    Real constraintX = 0.5;
    Real constraintValue = 0.7;
    Eigen::VectorXd basis = segment->evaluateAll(constraintX, -1, BSplineSegment::SideRight);
    
    cout << "\n2. Basis at x=" << constraintX << ":" << endl;
    cout << "   ";
    for (int i = 0; i < basis.size(); ++i) {
        cout << "B[" << i << "]=" << basis[i] << " ";
    }
    cout << endl;
    
    // Build inequality constraint: basis * coeffs >= constraintValue
    // In SCS format: -basis * coeffs <= -constraintValue (negate for <=)
    vector<Eigen::Triplet<Real>> A_triplets;
    for (Size i = 0; i < numVars; ++i) {
        if (std::abs(basis[i]) > 1e-10) {
            // Note: SCS expects Ax <= b for inequalities, so we negate
            A_triplets.push_back(Eigen::Triplet<Real>(0, i, -basis[i]));
        }
    }
    
    cout << "\n3. Inequality constraint (SCS format: Ax <= b):" << endl;
    cout << "   A row: ";
    for (const auto& t : A_triplets) {
        cout << "[" << t.col() << "]=" << t.value() << " ";
    }
    cout << endl;
    cout << "   b value: " << -constraintValue << " (negated)" << endl;
    
    // Create constraint system
    Eigen::SparseMatrix<Real> P(numVars, numVars);
    if (useLS) {
        // Small regularization for LS mode
        P.setIdentity();
        P *= 1e-6;
    } else {
        // Zero P for hard mode (no objective, just constraints)
        // Actually, we need some regularization even in hard mode
        P.setIdentity();
        P *= 1e-8;
    }
    
    vector<Real> b = {-constraintValue};  // RHS for inequality
    vector<Real> c(numVars, 0.0);  // No linear term
    
    vector<SplineConstraints::ConstraintType> types = {
        SplineConstraints::ConstraintType::LessEqual  // Inequality
    };
    
    cout << "\n4. Creating SplineConstraints:" << endl;
    cout << "   P matrix: " << (useLS ? "Identity * 1e-6" : "Identity * 1e-8") << endl;
    cout << "   Constraints: 1 inequality" << endl;
    cout << "   fitData: " << (useLS ? "true" : "false") << endl;
    
    auto constraints = ext::make_shared<SplineConstraints>(
        numVars, P, A_triplets, b, c, types, useLS
    );
    
    cout << "   Constraints created: " << constraints->getNConstraints() << " constraints" << endl;
    cout << "   Variables: " << constraints->getNumVariables() << endl;
    
    // Create BSplineStructure
    vector<ext::shared_ptr<BSplineSegment>> segments = {segment};
    BSplineStructure structure(segments, constraints, false, false);
    
    // Data points for interpolation
    vector<Real> x, y;
    if (useLS) {
        // Multiple points for LS
        x = {0.0, 0.25, 0.75, 1.0};
        y = {0.0, 0.25, 0.75, 1.0};
    } else {
        // Minimal points for hard mode
        x = {0.0, 1.0};
        y = {0.0, 1.0};
    }
    
    cout << "\n5. Interpolation data:" << endl;
    cout << "   Points: ";
    for (Size i = 0; i < x.size(); ++i) {
        cout << "(" << x[i] << "," << y[i] << ") ";
    }
    cout << endl;
    
    // Interpolate
    cout << "\n6. Calling interpolate..." << endl;
    Eigen::VectorXd solution;
    try {
        solution = structure.interpolate(x, y);
        cout << "   Interpolation successful" << endl;
    } catch (const std::exception& e) {
        cout << "   [ERROR] Interpolation failed: " << e.what() << endl;
        return;
    }
    
    cout << "\n7. Solution coefficients:" << endl;
    for (int i = 0; i < solution.size(); ++i) {
        cout << "   c[" << i << "] = " << solution[i] << endl;
    }
    
    // Evaluate at constraint point
    Real value = segment->value(solution, constraintX, 0, BSplineSegment::SideRight);
    
    cout << "\n8. Evaluation at x=" << constraintX << ":" << endl;
    cout << "   y(" << constraintX << ") = " << value << endl;
    cout << "   Constraint: y >= " << constraintValue << endl;
    cout << "   Satisfied: " << (value >= constraintValue - 1e-6 ? "YES" : "NO") << endl;
    
    if (value < constraintValue - 1e-6) {
        cout << "   [FAIL] Constraint violated by " << (constraintValue - value) << endl;
        
        // Debug: manually check constraint
        Real manual = 0;
        for (int i = 0; i < basis.size(); ++i) {
            manual += basis[i] * solution[i];
        }
        cout << "   Manual evaluation: " << manual << endl;
    } else {
        cout << "   [PASS] Constraint satisfied" << endl;
    }
}

int main() {
    cout << "Inequality Constraint Debug Test (Detailed)" << endl;
    cout << "===========================================" << endl;
    
    try {
        // Test hard mode (should fail)
        testInequalityConstraint(false);
        
        // Test LS mode (should work)
        testInequalityConstraint(true);
        
    } catch (const std::exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        return 1;
    }
    
    cout << "\n===========================================" << endl;
    cout << "Summary: Hard mode ignores inequality constraints" << endl;
    cout << "This appears to be because in hard mode, when there" << endl;
    cout << "are no interpolation constraints, the inequality" << endl;
    cout << "constraints are not being properly enforced." << endl;
    cout << "===========================================" << endl;
    
    return 0;
}