/*
 * Standalone C++ test for debugging inequality constraints in hard mode
 * Compile with: cl /I.. /EHsc debug_inequality.cpp ../lib/QuantLib-x64-mt.lib
 * Or use build_debug.bat
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
        cout << "[" << i << "]=" << vec[i] << " ";
    }
    cout << endl;
}

void testInequalityConstraint() {
    cout << "========================================" << endl;
    cout << "TEST: Inequality constraint y(0.5) >= 0.7" << endl;
    cout << "Linear interpolation from (0,0) to (1,1)" << endl;
    cout << "========================================" << endl;

    // Create linear B-spline segment
    // Knots: [0, 0.5, 1]
    vector<Real> knots = {0.0, 0.5, 1.0};
    Size degree = 1;
    
    // For linear B-spline with 3 knots, we need the knot vector with multiplicities
    // Full knot vector: [0, 0, 0.5, 1, 1]
    vector<Size> knotIndices = {0, 0, 1, 2, 2};
    
    auto segment = ext::make_shared<BSplineSegment>(
        knots, degree,
        vector<Integer>(),  // default knot indices
        BSplineSegment::SmoothnessDefault,
        BSplineSegment::TransformDefault,
        BSplineSegment::SideRight,
        1,    // requiredPoints
        true  // isGlobal
    );
    Size nVars = segment->getNumVariables();
    
    cout << "\nSegment info:" << endl;
    cout << "  Knots: ";
    for (auto k : knots) cout << k << " ";
    cout << endl;
    cout << "  Degree: " << degree << endl;
    cout << "  Number of variables (basis functions): " << nVars << endl;
    
    // Data points to interpolate
    vector<Real> xData = {0.0, 1.0};
    vector<Real> yData = {0.0, 1.0};
    
    cout << "\nData points:" << endl;
    for (size_t i = 0; i < xData.size(); ++i) {
        cout << "  (" << xData[i] << ", " << yData[i] << ")" << endl;
    }
    
    // Create structure with single segment
    vector<ext::shared_ptr<BSplineSegment>> segments = {segment};
    
    // Test both modes
    vector<string> modes = {"hard", "ls"};
    
    for (const auto& mode : modes) {
        cout << "\n========================================" << endl;
        cout << "Mode: " << mode << endl;
        cout << "========================================" << endl;
        
        // Create constraint system
        // Empty objective initially
        Eigen::SparseMatrix<Real> P(nVars, nVars);
        Eigen::VectorXd c = Eigen::VectorXd::Zero(nVars);
        
        // Start with inequality constraint: y(0.5) >= 0.7
        // At x=0.5, for linear B-spline: B0(0.5)=0.5, B1(0.5)=0.5, B2(0.5)=0
        // So y(0.5) = 0.5*c[0] + 0.5*c[1] >= 0.7
        // In standard form: -0.5*c[0] - 0.5*c[1] <= -0.7
        
        Eigen::SparseMatrix<Real> A_ineq(1, nVars);
        A_ineq.coeffRef(0, 0) = -0.5;
        A_ineq.coeffRef(0, 1) = -0.5;
        
        Eigen::VectorXd b_ineq(1);
        b_ineq[0] = -0.7;
        
        cout << "\nInequality constraint (before interpolation):" << endl;
        printMatrix("A_ineq", A_ineq);
        printVector("b_ineq", b_ineq);
        
        // Add interpolation constraints
        Eigen::SparseMatrix<Real> A_interp(2, nVars);
        Eigen::VectorXd b_interp(2);
        
        // At x=0: B0(0)=1, B1(0)=0, B2(0)=0
        A_interp.coeffRef(0, 0) = 1.0;
        b_interp[0] = 0.0;
        
        // At x=1: B0(1)=0, B1(1)=0, B2(1)=1
        A_interp.coeffRef(1, 2) = 1.0;
        b_interp[1] = 1.0;
        
        cout << "\nInterpolation constraints:" << endl;
        printMatrix("A_interp", A_interp);
        printVector("b_interp", b_interp);
        
        if (mode == "ls") {
            // For LS mode, interpolation goes into objective
            P = A_interp.transpose() * A_interp;
            c = -A_interp.transpose() * b_interp;
            
            cout << "\nLS objective:" << endl;
            printMatrix("P", P);
            printVector("c", c);
            
            // Only inequality remains as constraint
            auto constraints = ext::make_shared<SplineConstraints>(
                nVars, P, A_ineq, b_ineq, c,
                0,  // 0 equalities
                1,  // 1 inequality
                false  // not SCS ordered yet
            );
            
            cout << "\nCreated SplineConstraints for LS mode" << endl;
            cout << "  Equalities: 0" << endl;
            cout << "  Inequalities: 1" << endl;
            
        } else {
            // For hard mode, interpolation is equality constraint
            // Combine all constraints
            Eigen::SparseMatrix<Real> A_all(3, nVars);
            Eigen::VectorXd b_all(3);
            
            // Equalities first (SCS order)
            A_all.topRows(2) = A_interp;
            b_all.head(2) = b_interp;
            
            // Then inequality
            A_all.bottomRows(1) = A_ineq;
            b_all.tail(1) = b_ineq;
            
            cout << "\nCombined constraints (SCS ordered):" << endl;
            printMatrix("A_all", A_all);
            printVector("b_all", b_all);
            
            auto constraints = ext::make_shared<SplineConstraints>(
                nVars, P, A_all, b_all, c,
                2,  // 2 equalities
                1,  // 1 inequality
                false  // not SCS ordered flag (we did it manually)
            );
            
            cout << "\nCreated SplineConstraints for hard mode" << endl;
            cout << "  Equalities: 2" << endl;
            cout << "  Inequalities: 1" << endl;
        }
        
        // Note: Can't create structure and model without proper constraint object
        // This would need the actual constraint creation code
        
        cout << "\nSolving..." << endl;
        
        try {
            auto interpolation = model.interpolate(xData, yData);
            
            cout << "\nSolution found!" << endl;
            
            // Evaluate at key points
            Real y0 = interpolation(0.0);
            Real y05 = interpolation(0.5);
            Real y1 = interpolation(1.0);
            
            cout << "\nEvaluation:" << endl;
            cout << "  f(0.0) = " << y0 << " (expected: 0.0)" << endl;
            cout << "  f(0.5) = " << y05 << " (constraint: >= 0.7)" << endl;
            cout << "  f(1.0) = " << y1 << " (expected: 1.0)" << endl;
            
            // Check constraint satisfaction
            bool constraintSatisfied = (y05 >= 0.7 - 1e-10);
            cout << "\nConstraint y(0.5) >= 0.7: " 
                 << (constraintSatisfied ? "SATISFIED" : "VIOLATED") << endl;
            
            // Check interpolation accuracy
            Real err0 = std::abs(y0 - 0.0);
            Real err1 = std::abs(y1 - 1.0);
            cout << "Interpolation errors: |f(0)-0| = " << err0 
                 << ", |f(1)-1| = " << err1 << endl;
            
        } catch (const exception& e) {
            cout << "\nERROR: " << e.what() << endl;
        }
    }
}

int main() {
    try {
        testInequalityConstraint();
    } catch (const exception& e) {
        cerr << "Unexpected error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}