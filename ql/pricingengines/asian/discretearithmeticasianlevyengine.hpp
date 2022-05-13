/*
 Copyright (C) 2020 Skandinaviska Enskilda Banken AB (publ)

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file ql/pricingengines/asian/discretearithmeticasianlevyengine.hpp
    \brief Levy moment-matching discrete Asian option Engine
    \ingroup asianengines
*/

#ifndef quantlib_discrete_arithmetic_asian_levy_engine_hpp
#define quantlib_discrete_arithmetic_asian_levy_engine_hpp

#include <ql/instruments/asianoption.hpp>
#include <ql/processes/blackscholesprocess.hpp>

namespace QuantLib {

    /*! Levy two moment-matching discrete Asian option engine
        Analytical pricing based on the two-moment Levy approximation.
        References: "Option Pricing Formulas, Second Edition", E.G. Haug, 2006, pp. 192-195.

        \test
        - the correctness of the returned value is tested by reproducing
          results in literature with flat volatility term structures.
        - the pricing of trades with guaranteed exercise/OTM is also tested.
    */
    class DiscreteArithmeticAsianLevyEngine : public DiscreteAveragingAsianOption::engine {
      public:
        explicit DiscreteArithmeticAsianLevyEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess> process)
        : process_(std::move(process)) {
            registerWith(process_);
        }

        void calculate() const override;

      private:
        ext::shared_ptr<GeneralizedBlackScholesProcess> process_;
    };

}

#endif
