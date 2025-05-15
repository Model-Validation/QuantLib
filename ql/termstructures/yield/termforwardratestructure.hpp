/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004 StatPro Italia srl
 Copyright (C) 2009 Ferdinando Ametrano
 Copyright (C) 2024 SEB AB Sverrir Thorvaldsson

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

/*! \file forwardstructure.hpp
    \brief Forward-based yield term structure
*/

#ifndef quantlib_term_forward_rate_structure_hpp
#define quantlib_term_forward_rate_structure_hpp

#include "ql/indexes/iborindex.hpp"
#include <ql/termstructures/yieldtermstructure.hpp>

namespace QuantLib {

    //! %Termed forward-rate term structure
    /*! This abstract class acts as an adapter to YieldTermStructure allowing
        the programmer to implement only the <tt>forwardImpl(Time)</tt> method
        in derived classes.

        Zero yields and discounts cannot be calculated from termed forwards.

        Termed forward rates are assumed to be simple compounding.

        \ingroup yieldtermstructures
    */
    class TermForwardRateStructure : public YieldTermStructure {
      public:
        /*! \name Constructors
            See the TermStructure documentation for issues regarding
            constructors.
        */
        //@{
        explicit TermForwardRateStructure(const ext::shared_ptr<IborIndex>& index);
        explicit TermForwardRateStructure(const Date& referenceDate,
                                          const ext::shared_ptr<IborIndex>& index);

        Rate zeroYieldImpl(Time t) const override;
        DiscountFactor discountImpl(Time t) const override;
        //@}
    };

    inline Rate TermForwardRateStructure::zeroYieldImpl(Time t) const {
        QL_FAIL("Cannot get zero yield rate for a termed forward curve");
    }

    inline DiscountFactor TermForwardRateStructure::discountImpl(Time t) const {
        QL_FAIL("Cannot get discount factor for a termed forward curve");
    }
}

#endif
