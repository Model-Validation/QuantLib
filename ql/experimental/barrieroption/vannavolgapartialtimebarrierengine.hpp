
#ifndef quantlib_vanna_volga_partial_time_barrier_engine_hpp
#define quantlib_vanna_volga_partial_time_barrier_engine_hpp

#include <ql/experimental/fx/deltavolquote.hpp>
#include <ql/handle.hpp>
#include <ql/instruments/partialtimebarrieroption.hpp>
#include <ql/quote.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/types.hpp>

namespace QuantLib {
    class VannaVolgaPartialTimeBarrierEngine : public PartialTimeBarrierOption::engine {
      public:
        VannaVolgaPartialTimeBarrierEngine(Handle<DeltaVolQuote> atmVol,
                                           Handle<DeltaVolQuote> vol25Put,
                                           Handle<DeltaVolQuote> vol25Call,
                                           Handle<Quote> spotFX,
                                           Handle<YieldTermStructure> domesticTS,
                                           Handle<YieldTermStructure> foreignTS,
                                           bool adaptVanDelta = false,
                                           Real bsPriceWithSmile = 0.0);

        void calculate() const override;

      private:
        const Handle<DeltaVolQuote> atmVol_;
        const Handle<DeltaVolQuote> vol25Put_;
        const Handle<DeltaVolQuote> vol25Call_;
        const Time T_;
        const Handle<Quote> spotFX_;
        const Handle<YieldTermStructure> domesticTS_;
        const Handle<YieldTermStructure> foreignTS_;
        const bool adaptVanDelta_;
        const Real bsPriceWithSmile_;
    };
}

#endif
