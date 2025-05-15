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
#include <ql/termstructures/yield/termforwardratestructure.hpp>

namespace QuantLib {

    TermForwardRateStructure::TermForwardRateStructure(const ext::shared_ptr<IborIndex>& index)
    : YieldTermStructure(index) {
        isTermForward_ = true;
        supportsDiscount_ = false;
        supportsZero_ = false;
    }

    TermForwardRateStructure::TermForwardRateStructure(const Date& referenceDate,
        const ext::shared_ptr<IborIndex>& index)
    : YieldTermStructure(referenceDate, index) {
        isTermForward_ = true;
        supportsDiscount_ = false;
        supportsZero_ = false;
    }
}
