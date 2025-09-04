/* Simplified debug inspection methods for B-spline structures */

#ifndef quantlib_bspline_debug_inspection_hpp
#define quantlib_bspline_debug_inspection_hpp

#include <ql/math/interpolations/bsplineinterpolation/bsplineinterpolation.hpp>
#include <ql/math/interpolations/bsplineinterpolation/bsplinestructure.hpp>
#include <ql/math/interpolations/bsplineinterpolation/splineconstraints.hpp>
#include <vector>
#include <iostream>

namespace QuantLib {

    // Simplified inspection that only uses public methods
    class BSplineDebugInspector {
    public:
        // For BSplineStructure - get basic info
        static Size getNumberOfSegments(const BSplineStructure& structure) {
            // Use public method if available
            return structure.getSplineSegmentsSwig().size();
        }
        
        // For BSplineInterpolation - already has public get_coefficients
        static std::vector<Real> getSolutionCoefficients(
            const BSplineInterpolation& interp) {
            return interp.get_coefficients();
        }
        
        // Print basic info
        static void printBasicInfo(const BSplineStructure& structure) {
            std::cout << "\n=== BSplineStructure Info ===" << std::endl;
            std::cout << "Number of segments: " << getNumberOfSegments(structure) << std::endl;
            std::cout << "Total variables: " << structure.getNumVariablesSwig() << std::endl;
        }
    };

}

#endif