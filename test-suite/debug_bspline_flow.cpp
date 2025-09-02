/*
 * Debug test to trace the flow through BSplineStructure interpolation
 * This tests the exact scenario from Python: quadratic B-spline with inequality
 */

#include <ql/math/interpolations/bsplineinterpolation/bsplineinterpolation.hpp>
#include <ql/math/interpolations/bsplineinterpolation/bsplinestructure.hpp>
#include <ql/math/interpolations/bsplineinterpolation/splinesegment.hpp>
#include <ql/math/interpolations/bsplineinterpolation/splineconstraints.hpp>
#include <ql/shared_ptr.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace QuantLib;
using namespace std;

void testBSplineFlow() {
    cout << "========================================" << endl;
    cout << "B-Spline Flow Test - Quadratic with Inequality" << endl;
    cout << "========================================" << endl;
    
    // Quadratic B-spline - use simple knots and indices
    vector<Real> simpleKnots = {0.0, 1.0};  // Just the unique knots
    Size degree = 2;
    vector<Integer> knotIndices = {0, 0, 0, 1, 1, 1};  // Repeat for degree+1
    
    // Create segment
    auto segment = ext::make_shared<BSplineSegment>(
        simpleKnots, degree,
        knotIndices,
        BSplineSegment::SmoothnessDefault,
        BSplineSegment::TransformDefault,
        BSplineSegment::SideRight,
        1,    // requiredPoints
        true  // isGlobal
    );
    
    Size n_vars = segment->getNumVariables();
    cout << "\n1. Segment created:" << endl;
    cout << "   Degree: " << degree << endl;
    cout << "   Simple knots: ";
    for (Real k : simpleKnots) cout << k << " ";
    cout << endl;
    cout << "   Knot indices: ";
    for (Integer i : knotIndices) cout << i << " ";
    cout << endl;
    cout << "   Variables: " << n_vars << " (should be 3 for quadratic)" << endl;
    
    // Evaluate basis functions at key points
    cout << "\n2. Basis functions:" << endl;
    vector<Real> test_x = {0.0, 0.5, 1.0};
    for (Real x : test_x) {
        Eigen::VectorXd basis = segment->evaluateAll(x, -1, BSplineSegment::SideRight);
        cout << "   At x=" << x << ": ";
        for (int i = 0; i < basis.size(); ++i) {
            cout << "B" << i << "=" << basis[i] << " ";
        }
        cout << endl;
    }
    
    // Create inequality constraint: y(0.5) >= 0.7
    cout << "\n3. Creating inequality constraint y(0.5) >= 0.7:" << endl;
    
    Eigen::VectorXd basis_05 = segment->evaluateAll(0.5, -1, BSplineSegment::SideRight);
    cout << "   Basis at x=0.5: ";
    for (int i = 0; i < basis_05.size(); ++i) {
        if (abs(basis_05[i]) > 1e-10) {
            cout << "B" << i << "=" << basis_05[i] << " ";
        }
    }
    cout << endl;
    
    // For quadratic B-spline, at x=0.5:
    // B0(0.5) = 0.25, B1(0.5) = 0.5, B2(0.5) = 0.25
    // So constraint is: 0.25*c0 + 0.5*c1 + 0.25*c2 >= 0.7
    // In SCS form: -(0.25*c0 + 0.5*c1 + 0.25*c2) <= -0.7
    
    vector<Eigen::Triplet<double>> A_triplets;
    for (int i = 0; i < n_vars; ++i) {
        if (abs(basis_05[i]) > 1e-10) {
            A_triplets.push_back(Eigen::Triplet<double>(0, i, -basis_05[i]));
        }
    }
    
    cout << "   Constraint in SCS form: ";
    for (const auto& t : A_triplets) {
        cout << "A[0," << t.col() << "]=" << t.value() << " ";
    }
    cout << "b[0]=-0.7" << endl;
    
    // Create SplineConstraints with inequality
    Eigen::SparseMatrix<double> P(n_vars, n_vars);
    P.setIdentity();
    P *= 1e-8;  // Small regularization
    
    vector<Real> b = {-0.7};
    vector<Real> c(n_vars, 0.0);
    
    vector<SplineConstraints::ConstraintType> types = {
        SplineConstraints::ConstraintType::LessEqual
    };
    
    cout << "\n4. Creating SplineConstraints:" << endl;
    auto constraints = ext::make_shared<SplineConstraints>(
        n_vars, P, A_triplets, b, c, types, false  // fitData=false (hard mode)
    );
    cout << "   Constraints created with " << constraints->getNConstraints() << " inequality" << endl;
    
    // Create BSplineStructure
    cout << "\n5. Creating BSplineStructure:" << endl;
    vector<ext::shared_ptr<BSplineSegment>> segments = {segment};
    BSplineStructure structure(segments, constraints, false, false);
    cout << "   Structure created" << endl;
    
    // Test Case 1: Minimal interpolation (2 points for 3 DOF)
    cout << "\n6. Test Case 1: Minimal interpolation (underdetermined)" << endl;
    vector<Real> x_data = {0.0, 1.0};
    vector<Real> y_data = {0.0, 1.0};
    
    cout << "   Data points: ";
    for (Size i = 0; i < x_data.size(); ++i) {
        cout << "(" << x_data[i] << "," << y_data[i] << ") ";
    }
    cout << endl;
    
    cout << "   This creates an underdetermined system:" << endl;
    cout << "     - 3 variables (c0, c1, c2)" << endl;
    cout << "     - 2 interpolation constraints (f(0)=0, f(1)=1)" << endl;
    cout << "     - 1 inequality constraint (f(0.5)>=0.7)" << endl;
    cout << "   The inequality should guide the solution" << endl;
    
    // Interpolate
    cout << "\n   Calling interpolate..." << endl;
    Eigen::VectorXd solution;
    try {
        solution = structure.interpolate(x_data, y_data);
        cout << "   Interpolation succeeded" << endl;
    } catch (const std::exception& e) {
        cout << "   ERROR: " << e.what() << endl;
        return;
    }
    
    cout << "\n   Solution coefficients:" << endl;
    for (int i = 0; i < solution.size(); ++i) {
        cout << "     c[" << i << "] = " << solution[i] << endl;
    }
    
    // Evaluate at test points
    cout << "\n   Evaluation:" << endl;
    for (Real x : test_x) {
        Real val = segment->value(solution, x, 0, BSplineSegment::SideRight);
        cout << "     f(" << x << ") = " << val;
        if (x == 0.5) {
            cout << " (constraint: >= 0.7) " << (val >= 0.7 - 1e-6 ? "[OK]" : "[FAIL]");
        }
        cout << endl;
    }
    
    // Test Case 2: With 3 points (exactly determined before inequality)
    cout << "\n7. Test Case 2: With 3 interpolation points" << endl;
    x_data = {0.0, 0.5, 1.0};
    y_data = {0.0, 0.5, 1.0};  // Natural linear interpolation
    
    cout << "   Data points: ";
    for (Size i = 0; i < x_data.size(); ++i) {
        cout << "(" << x_data[i] << "," << y_data[i] << ") ";
    }
    cout << endl;
    
    cout << "   This would naturally give f(0.5)=0.5" << endl;
    cout << "   But inequality requires f(0.5)>=0.7" << endl;
    cout << "   System is overdetermined - should fail or adjust" << endl;
    
    // Create fresh structure with same constraints
    BSplineStructure structure2(segments, constraints, false, false);
    
    try {
        solution = structure2.interpolate(x_data, y_data);
        cout << "   Interpolation succeeded (unexpected?)" << endl;
        
        cout << "\n   Solution coefficients:" << endl;
        for (int i = 0; i < solution.size(); ++i) {
            cout << "     c[" << i << "] = " << solution[i] << endl;
        }
        
        Real val_05 = segment->value(solution, 0.5, 0, BSplineSegment::SideRight);
        cout << "   f(0.5) = " << val_05 << " " << (val_05 >= 0.7 - 1e-6 ? "[Constraint OK]" : "[Constraint VIOLATED]") << endl;
        
    } catch (const std::exception& e) {
        cout << "   Expected failure: " << e.what() << endl;
    }
}

int main() {
    cout << "B-Spline Flow Debug Test" << endl;
    cout << "========================" << endl;
    
    try {
        testBSplineFlow();
    } catch (const std::exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        return 1;
    }
    
    cout << "\n========================" << endl;
    cout << "Summary:" << endl;
    cout << "The bug appears when BSplineStructure::interpolate adds" << endl;
    cout << "interpolation nodes as equality constraints in hard mode." << endl;
    cout << "Pre-existing inequality constraints are not properly" << endl;
    cout << "combined with the new equality constraints." << endl;
    cout << "========================" << endl;
    
    return 0;
}