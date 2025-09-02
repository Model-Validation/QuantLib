/*
 * Direct test of SplineConstraints to isolate inequality constraint bug
 * This bypasses BSplineStructure to test constraint handling directly
 */

#include <ql/math/interpolations/bsplineinterpolation/splineconstraints.hpp>
#include <ql/shared_ptr.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace QuantLib;
using namespace std;

void testDirectConstraints() {
    cout << "========================================" << endl;
    cout << "Direct SplineConstraints Test" << endl;
    cout << "========================================" << endl;
    
    // Simple 3-variable problem: c[0], c[1], c[2]
    Size n = 3;
    
    // Test 1: Only inequality constraint
    cout << "\nTest 1: Only inequality constraint -c[1] <= -0.7" << endl;
    cout << "========================================" << endl;
    
    // P matrix (small regularization)
    Eigen::SparseMatrix<Real> P(n, n);
    P.setIdentity();
    P *= 1e-8;
    
    // Inequality: -c[1] <= -0.7 (c[1] >= 0.7)
    vector<Eigen::Triplet<double>> A_triplets;
    A_triplets.push_back(Eigen::Triplet<double>(0, 1, -1.0));
    
    vector<Real> b = {-0.7};
    vector<Real> c(n, 0.0);
    
    vector<SplineConstraints::ConstraintType> types = {
        SplineConstraints::ConstraintType::LessEqual
    };
    
    cout << "Creating constraints with:" << endl;
    cout << "  Variables: " << n << endl;
    cout << "  Inequalities: 1" << endl;
    cout << "  Constraint: -c[1] <= -0.7" << endl;
    
    auto constraints = ext::make_shared<SplineConstraints>(
        n, P, A_triplets, b, c, types, false  // fitData=false
    );
    
    cout << "Constraints created successfully" << endl;
    cout << "  Total constraints: " << constraints->getNConstraints() << endl;
    
    // Solve without parameters
    cout << "\nSolving..." << endl;
    int status = constraints->solve();
    
    if (status == 1) {
        cout << "Solution found!" << endl;
        Eigen::VectorXd solution = constraints->getSolution();
        for (int i = 0; i < solution.size(); ++i) {
            cout << "  c[" << i << "] = " << solution[i] << endl;
        }
        
        // Check constraint
        Real constraint_value = -solution[1];
        cout << "\nConstraint check: -c[1] = " << constraint_value 
             << " <= -0.7? " << (constraint_value <= -0.7 + 1e-6 ? "YES" : "NO") << endl;
    } else {
        cout << "Solver failed with status: " << status << endl;
    }
    
    // Test 2: Inequality + equality constraints
    cout << "\n\nTest 2: Mixed constraints" << endl;
    cout << "========================================" << endl;
    cout << "Constraints in SCS order:" << endl;
    cout << "  EQ0: c[0] = 0" << endl;
    cout << "  EQ1: c[2] = 1" << endl;
    cout << "  INEQ0: -c[1] <= -0.7" << endl;
    
    // Build constraint matrix in SCS order (equalities first)
    vector<Eigen::Triplet<double>> A_mixed_triplets;
    // Equality: c[0] = 0
    A_mixed_triplets.push_back(Eigen::Triplet<double>(0, 0, 1.0));
    // Equality: c[2] = 1
    A_mixed_triplets.push_back(Eigen::Triplet<double>(1, 2, 1.0));
    // Inequality: -c[1] <= -0.7
    A_mixed_triplets.push_back(Eigen::Triplet<double>(2, 1, -1.0));
    
    vector<Real> b_mixed = {0.0, 1.0, -0.7};
    
    // Create with 2 equalities, 1 inequality
    Size numEq = 2;
    Size numIneq = 1;
    
    cout << "\nCreating mixed constraints:" << endl;
    cout << "  Variables: " << n << endl;
    cout << "  Equalities: " << numEq << endl;
    cout << "  Inequalities: " << numIneq << endl;
    
    // For the constructor with numEq/numIneq, we need to convert triplets to dense matrix
    vector<vector<double>> A_dense(numEq + numIneq, vector<double>(n, 0.0));
    for (const auto& t : A_mixed_triplets) {
        A_dense[t.row()][t.col()] = t.value();
    }
    
    // Convert P to dense as well
    vector<vector<double>> P_dense(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        P_dense[i][i] = 1e-8;
    }
    
    auto constraints2 = ext::make_shared<SplineConstraints>(
        n, P_dense, A_dense, b_mixed, c,
        numEq, numIneq, false,  // fitData=false
        1e-9, 1e-9, 1e-6  // tolerances
    );
    
    cout << "Constraints created successfully" << endl;
    cout << "  Total constraints: " << constraints2->getNConstraints() << endl;
    
    // Solve
    cout << "\nSolving mixed system..." << endl;
    status = constraints2->solve();
    
    if (status == 1) {
        cout << "Solution found!" << endl;
        Eigen::VectorXd solution = constraints2->getSolution();
        for (int i = 0; i < solution.size(); ++i) {
            cout << "  c[" << i << "] = " << solution[i] << endl;
        }
        
        // Check all constraints
        cout << "\nConstraint checks:" << endl;
        cout << "  c[0] = " << solution[0] << " (should be 0): " 
             << (abs(solution[0]) < 1e-6 ? "OK" : "FAIL") << endl;
        cout << "  c[2] = " << solution[2] << " (should be 1): " 
             << (abs(solution[2] - 1.0) < 1e-6 ? "OK" : "FAIL") << endl;
        cout << "  c[1] = " << solution[1] << " >= 0.7? " 
             << (solution[1] >= 0.7 - 1e-6 ? "OK" : "FAIL") << endl;
    } else {
        cout << "Solver failed with status: " << status << endl;
    }
}

int main() {
    cout << "Direct SplineConstraints Test" << endl;
    cout << "=============================" << endl;
    
    try {
        testDirectConstraints();
    } catch (const std::exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        return 1;
    }
    
    cout << "\n=============================" << endl;
    cout << "Test completed" << endl;
    
    return 0;
}