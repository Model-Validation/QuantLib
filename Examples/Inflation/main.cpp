#include <ql/currencies/europe.hpp>
#include <ql/indexes/inflationindex.hpp>
#include <ql/indexes/region.hpp>
#include <ql/indexes/inflation/ukrpi.hpp>
#include <ql/termstructures/bootstraphelper.hpp>
#include <ql/termstructures/inflation/inflationhelpers.hpp>
#include <ql/termstructures/inflation/piecewisezeroinflationcurve.hpp>
#include <ql/termstructures/inflation/piecewiseyoyinflationcurve.hpp>
#include <ql/termstructures/inflationtermstructure.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/calendars/unitedkingdom.hpp>
#include <ql/time/date.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/daycounters/thirty360.hpp>
#include <ql/time/frequency.hpp>
#include <ql/time/period.hpp>
#include <ql/time/schedule.hpp>
#include <ql/types.hpp>

#include <matplot/matplot.h>

#include <vector>
#include <iostream>


#define LENGTH(a) (sizeof(a) / sizeof(a[0]))


using namespace QuantLib;

// Coupon SWAP data
struct Coupon {
    Date maturity;
    Rate rate;
};

ext::shared_ptr<YieldTermStructure> nominalTermStructure(Date studyDate, Rate baseRate) {
    return ext::shared_ptr<YieldTermStructure>(
            new FlatForward(studyDate, baseRate, Actual360()));
}

template <class T>
std::vector<ext::shared_ptr<BootstrapHelper<T>>>
makeHelpers(const std::vector<Coupon>& coupons,
            std::function<ext::shared_ptr<BootstrapHelper<T>>(const Handle<Quote>&, const Date&)>
                makeHelper) {

    std::vector<ext::shared_ptr<BootstrapHelper<T>>> instruments;
    for (Coupon coupon : coupons) {
        Handle<Quote> quote(ext::shared_ptr<Quote>(new SimpleQuote(coupon.rate / 100.0)));
        ext::shared_ptr<BootstrapHelper<T>> anInstrument = makeHelper(quote, coupon.maturity);
        instruments.push_back(anInstrument);
    }

    return instruments;
}

int main() {
    // Evaluation date, AsOfDate
    Date studyDate(13, February, 2024);
    std::cout << "studyDate:" << studyDate << '\n';
    Calendar calendar = UnitedKingdom();
    BusinessDayConvention bdc = ModifiedFollowing;
    DayCounter dc = Thirty360(Thirty360::European);
    Frequency frequency = Monthly;
    // Common to YoY and ZC
    Actual360 dc360;
    // Use the same for both YoY and ZC inflation swap curves
    Period observationLag = Period(2, Months);
    
    // KPI Index
    RelinkableHandle<ZeroInflationTermStructure> relinkableZeroTermStructureHandle;
    Period availabilityLag(25, Days);
    ZeroInflationIndex ukrpi("UKRPI", UKRegion(), false, Monthly, availabilityLag, GBPCurrency(), relinkableZeroTermStructureHandle);
    //UKRPI ukrpi(relinkableZeroTermStructureHandle);
    
    
    boost::shared_ptr<ZeroInflationIndex> pUkrpiIndex =
        boost::make_shared<ZeroInflationIndex>(ukrpi);

    Real rpi[] = {294.6, 296.0, 296.9, 301.1, 301.9, 304.0, 305.5, 307.4, 308.6, 312.0,
                  314.3, 317.7, 317.7, 320.2, 323.5, 334.6, 337.1, 340.0, 343.2, 345.2,
                  347.6, 356.2, 358.3, 360.4, 360.3, 364.5, 367.2, 372.8, 375.3, 376.4,
                  374.2, 376.6, 378.4, 377.8, 377.3, 379.0, 378.0};

    Date from(1, January, 2021);
    Date to(1, January, 2024);
    Schedule rpiSchedule = MakeSchedule().from(from).to(to).withFrequency(Monthly);

    for (Size i = 0; i < LENGTH(rpi); ++i) {
        ukrpi.addFixing(rpiSchedule[i], rpi[i]);
    }

    Date baseDate = pUkrpiIndex->lastFixingDate();
    std::cout << "baseDate: " << baseDate << '\n';
    TimeSeries<Real> ts = ukrpi.timeSeries();
    std::cout << "ts.size(): " << ts.size() << '\n';


    ///////////////
    // Seasonality
        std::vector<Rate> seasonalityFactors = {
        -0.820367,
        0.458629,
        0.003204,
        0.324084,
        0.060618,
        0.011060,
        -0.229372,
        0.320663,
        -0.126170,
        -0.166436,
        -0.119117,
        0.305323
    };

    Date seasonalityBaseDate  = to; // Where RPI stops
    ext::shared_ptr<MultiplicativePriceSeasonality> nonUnitSeasonality =
        ext::make_shared<MultiplicativePriceSeasonality>(seasonalityBaseDate, Monthly, seasonalityFactors);

    // YoY Inflation Swaps
    bool interpolatedYoy = false;
    YoYInflationIndex yyIndex(pUkrpiIndex, interpolatedYoy);
    boost::shared_ptr<YoYInflationIndex> pYoyIndex =
        boost::make_shared<YoYInflationIndex>(yyIndex);
        Rate yyRate = yyIndex.fixing(Date(1, January, 2024));
    std::cout << "1 Jan 2024: " << yyRate * 100 << '%' << std::endl;

    // YoY Swaps (first 2 years to get monthly granularity)

    std::vector<Coupon> yoySwap = {
        {Date(15, May, 2024), 4.25500},
        {Date(13, June, 2024), 2.98000},
        {Date(15, July, 2024), 2.67250},
        {Date(14, August, 2024), 2.66500},
        {Date(13, September, 2024), 3.38750},
        {Date(14, October, 2024), 3.29500},
        {Date(13, November, 2024), 2.90500},
        {Date(13, December, 2024), 3.32500},
        {Date(15, January, 2025), 3.61000},
        {Date(13, February, 2025), 3.64500},
        {Date(13, March, 2025), 3.44628},
        {Date(15, April, 2025), 3.29883},
        {Date(13, May, 2025), 3.10817},
        {Date(13, June, 2025), 3.72133},
        {Date(15, July, 2025), 3.75031},
        {Date(13, August, 2025), 3.65739},
        {Date(15, September, 2025), 3.78288},
        {Date(14, October, 2025), 3.81064},
        {Date(13, November, 2025), 3.84212},
        {Date(15, December, 2025), 3.77549},
        {Date(13, January, 2026), 3.7825},
    };


    //Rate baseYoYRate = 0.03645;//yyIndex.fixing(Date(1, January, 2024));
    Rate baseYoYRate = 0.056014;//yyIndex.fixing(Date(1, January, 2024));
    std::cout << "baseYoYRate: " << baseYoYRate << '\n';
    ext::shared_ptr<YieldTermStructure> yoyNominalTS = nominalTermStructure(studyDate, baseYoYRate);

    auto makeHelperYoY = [&](const Handle<Quote>& quote, const Date& maturity) {
        return ext::make_shared<YearOnYearInflationSwapHelper>(
                    quote,
                    observationLag,
                    maturity,
                    calendar,
                    bdc,
                    dc,
                    pYoyIndex,
                    Handle<YieldTermStructure>(yoyNominalTS));
    };

    auto helpersYoY = makeHelpers<YoYInflationTermStructure>(yoySwap, makeHelperYoY);

    auto pYoYIS = ext::make_shared<PiecewiseYoYInflationCurve<Linear>>(
        studyDate, baseDate, baseYoYRate, frequency, true, dc, helpersYoY);

    pYoYIS->update();
    pYoYIS->setSeasonality(nonUnitSeasonality);
    
    ///////////////
    // ZCIS //
    // Zero-Coupon Inflation Swaps

    std::vector<Coupon> zcSwap = {
        {Date(13, February, 2026), 3.75625},
        {Date(13, February, 2027), 3.78750},
        {Date(13, February, 2028), 3.79250},
        {Date(13, February, 2029), 3.79250},
        {Date(13, February, 2030), 3.80625},
        {Date(13, February, 2031), 3.75250},
        {Date(13, February, 2032), 3.69000},
        {Date(14, February, 2033), 3.64000},
        {Date(14, February, 2034), 3.60000},
        {Date(13, February, 2036), 3.53875},
        {Date(14, February, 2039), 3.47125},
        {Date(15, February, 2044), 3.39125},
        {Date(15, February, 2049), 3.29375},
        {Date(13, February, 2054), 3.23000},
        {Date(13, February, 2064), 3.12500},
        {Date(13, February, 2074), 3.10250},
    };


    // boost::shared_ptr<YieldTermStructure> yieldTermStructure(
    //        new FlatForward(studyDate, 0.05, Actual360()));
    //    Handle<YieldTermStructure> nts(yieldTermStructure);
    Rate zcBaseRate = 0.03756250;
    ext::shared_ptr<YieldTermStructure> nts = nominalTermStructure(studyDate, zcBaseRate);

    auto makeHelper = [&](const Handle<Quote>& quote, const Date& maturity) {
        return ext::make_shared<ZeroCouponInflationSwapHelper>(quote,
                                                               observationLag,
                                                               maturity,
                                                               calendar,
                                                               bdc,
                                                               dc,
                                                               pUkrpiIndex,
                                                               CPI::Linear,
                                                               Handle<YieldTermStructure>(nts));
    };

    std::vector<ext::shared_ptr<BootstrapHelper<ZeroInflationTermStructure>>> helpers =
        makeHelpers<ZeroInflationTermStructure>(zcSwap, makeHelper);

    ext::shared_ptr<PiecewiseZeroInflationCurve<Linear> > pZCIS =
        ext::make_shared<PiecewiseZeroInflationCurve<Linear>>(
            studyDate,
            baseDate,
            frequency,
            dc,
            helpers
            );
    
    relinkableZeroTermStructureHandle.linkTo(pZCIS);


    std::vector<Time> times = pZCIS->times();  // year fractions
    std::vector<Real> rates = pZCIS->data();
    std::vector<Date> c = pZCIS->dates();
    std::vector<std::pair<Date, Rate> > nodes = pZCIS->nodes();

    std::cout << "ZCIS baseDate(): " << pZCIS->baseDate() << '\n';
    std::cout << "ZCIS maxDate(): " << pZCIS->maxDate() << '\n';

    //////////////////////
    // Plot Inflation curve

    std::vector<Date> dates = {
        Date(1, January, 2024),  Date(1, February, 2024), Date(1, March, 2024),
        Date(1, April, 2024),    Date(1, May, 2024),      Date(1, June, 2024),
        Date(1, July, 2024),     Date(1, August, 2024),   Date(1, September, 2024),
        Date(1, October, 2024),  Date(1, November, 2024), Date(1, December, 2024),
        Date(1, January, 2025),  Date(1, February, 2025), Date(1, March, 2025),
        Date(1, April, 2025),    Date(1, May, 2025),      Date(1, June, 2025),
        Date(1, July, 2025),     Date(1, August, 2025),   Date(1, September, 2025),
        Date(1, October, 2025),  
        Date(1, November, 2025), 
        Date(1, December, 2025),
        Date(1, December, 2026), 
        Date(1, December, 2027), 
        Date(1, December, 2028), 
        Date(1, December, 2029),
        Date(1, December, 2030),
        Date(1, December, 2031),
        Date(1, December, 2032), 
        Date(1, December, 2033),
        Date(1, December, 2035), 
        Date(1, December, 2038), 
        Date(1, December, 2043),
        Date(1, December, 2048), 
        Date(1, December, 2053), 
        Date(1, December, 2063),
        Date(1, December, 2073),
    };

    //Date startDay(1, January, 2024); // Must be after or equal base date
    //Date toDay(1, September, 2033);
    //Schedule printSchedule = MakeSchedule().from(startDay).to(toDay).withFrequency(Monthly);

    std::vector<Time> x(dates.size());
    std::vector<Rate> y(dates.size());

    for (auto day : dates) {
        Time t = dc360.yearFraction(dates[0], day);
        Real spotRate;
        if (day < Date(1, November, 2025)) {
            std::cout << "YoY benchmarch:" << '\n';
            spotRate = pYoYIS->yoyRate(day);
        }
        else {
            std::cout << "ZC benchmarch:" << '\n';
            spotRate = pZCIS->zeroRate(day);
        }

        std::cout << "Date: " << day << '\n';
        std::cout << "Spot rate: " << spotRate << '\n';
        
        x.emplace_back(t);
        y.emplace_back(spotRate);
    }
    
   
    // inflation index
    Date iiStart(1, January, 2024); // Must be after or equal base date
    Date iiTo(1, September, 2073);
    Schedule inflationSchedule = MakeSchedule().from(iiStart).to(iiTo).withFrequency(Monthly);
    std::vector<Time> xt(inflationSchedule.size());
    std::vector<Rate> ii(inflationSchedule.size());
    std::cout << "printschedule.size: " << inflationSchedule.size() << '\n' ;
    for (auto day : inflationSchedule) {
        Time t = dc360.yearFraction(dates[0], day);
        Rate index = pUkrpiIndex->fixing(day);
        //std::cout << "Time fraction: " << day << '\n';
        //std::cout << "Inflation index: " << index << '\n';
        xt.emplace_back(t);
        ii.emplace_back(index);
    }

    //////// plot ////////
    using namespace matplot;

    figure();
    title("Spot rate, annual comp");
    xlabel("Years");
    plot(x, y, "-o");
    
    figure();
    title("Inflation index");
    xlabel("Years");
    plot(xt, ii, "-o");
    
    show();
    
    //std::cout << pZCIS->zeroRate(studyDate) << '\n';
    /*
    std::cout << "---times---" << '\n';
    for (auto time : times) {
        std::cout << time << '\n';
    }

    std::cout << "---dates---" << '\n';
    for (auto date : dates) {
        std::cout << date << '\n';
    }

    std::cout << "---rates---" << '\n';
    for (auto rate : rates) {
        std::cout << rate << '\n';
    }

    std::cout << "---nodes---" << '\n';
    for (auto node : nodes) {
        std::cout << node.first << '\n';
        std::cout << node.second * 100 << '%' << '\n';
    }*/

    return 0;
}
