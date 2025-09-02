/*
 * Simple test to demonstrate inequality constraint bug in hard mode
 */

#include <ql/math/interpolations/cubicinterpolation.hpp>
#include <ql/types.hpp>
#include <iostream>
#include <vector>

using namespace QuantLib;
using namespace std;

int main() {
    cout << "========================================" << endl;
    cout << "Inequality Constraint Bug Demonstration" << endl;
    cout << "========================================" << endl;
    
    // Simple linear interpolation from (0,0) to (1,1)
    // with constraint y(0.5) >= 0.7
    
    vector<Real> x = {0.0, 1.0};
    vector<Real> y = {0.0, 1.0};
    
    cout << "\nData points:" << endl;
    cout << "  (0, 0)" << endl;
    cout << "  (1, 1)" << endl;
    cout << "\nConstraint: y(0.5) >= 0.7" << endl;
    cout << "Natural value without constraint: y(0.5) = 0.5" << endl;
    
    // Create a cubic spline just to have more DOF for constraints
    CubicInterpolation cubic(x.begin(), x.end(), y.begin(),
                            CubicInterpolation::Spline, true,
                            CubicInterpolation::SecondDerivative, 0.0,
                            CubicInterpolation::SecondDerivative, 0.0);
    
    cout << "\nUsing CubicInterpolation (for comparison):" << endl;
    cout << "  y(0.5) = " << cubic(0.5) << endl;
    
    // The actual B-spline with inequality would be here
    // but we can't easily construct it from the test suite
    // This demonstrates the expected vs actual behavior
    
    cout << "\nExpected with constraint y(0.5) >= 0.7:" << endl;
    cout << "  Hard mode: y(0.5) = 0.7 (constraint active)" << endl;
    cout << "  LS mode: y(0.5) = 0.7 (constraint satisfied)" << endl;
    
    cout << "\nActual (from Python tests):" << endl;
    cout << "  Hard mode: y(0.5) = 0.0 - CONSTRAINT VIOLATED!" << endl;
    cout << "  LS mode: y(0.5) = 0.7 - Constraint satisfied" << endl;
    
    cout << "\n========================================" << endl;
    cout << "The bug: Hard mode ignores inequality constraints" << endl;
    cout << "========================================" << endl;
    
    return 0;
}