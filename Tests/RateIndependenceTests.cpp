// Rate independence of the analog-memory stages.
//
// The emulation stages store physical time constants (inductor tail, dielectric
// absorption, PSU sag, noise shaping). Those must be functions of *seconds*,
// not of the internal processing rate — otherwise switching the oversampling
// mode changes the sound of the emulation, which is the defect this suite
// guards against.
//
// Every case below drives the same input at two or three processing rates and
// compares the response on a common time grid.

#include "Tests/TestHarness.h"

#include "Dsp/Rate.h"
#include "Dsp/Stages.h"
#include "Dsp/Transformer.h"

#include <cmath>
#include <vector>

namespace {

constexpr double kRates[] = {48000.0 * 8.0, 48000.0 * 16.0, 48000.0 * 32.0};

// Sample a stage's step response at fixed wall-clock instants so responses at
// different rates can be compared directly.
//
// The sample at i == 0 is deliberately skipped: a discrete one-pole answers a
// step with (1 - pole) * x, which approximates the continuous response at
// t = dt and therefore scales with the rate. That is sampling-grid error, not
// a rate-dependent time constant, and it is gone by the next grid point.
template <typename Stage, typename Setup>
std::vector<double> stepResponse (double rate, float step, double durationSec,
                                  double sampleEverySec, Setup&& setup)
{
    Stage stage;
    setup (stage, rate);
    stage.reset();

    const long long total = static_cast<long long> (durationSec * rate);
    const long long stride = static_cast<long long> (sampleEverySec * rate);

    std::vector<double> samples;
    for (long long i = 0; i < total; ++i)
    {
        const double y = static_cast<double> (stage.process (step));
        if (stride > 0 && i > 0 && (i % stride) == 0)
            samples.push_back (y);
    }
    return samples;
}

double relativeSpread (const std::vector<std::vector<double>>& runs)
{
    double spread = 0.0;
    const size_t count = runs.front().size();
    for (size_t k = 0; k < count; ++k)
    {
        double lo = runs.front()[k], hi = runs.front()[k];
        for (const auto& run : runs)
        {
            lo = std::min (lo, run[k]);
            hi = std::max (hi, run[k]);
        }
        const double scale = std::max (1e-9, std::abs (hi) + std::abs (lo));
        spread = std::max (spread, 2.0 * (hi - lo) / scale);
    }
    return spread;
}

} // namespace

EON_TEST_CASE("rescalePole keeps the pole frequency across rates")
{
    constexpr double authored = 0.9;                 // 48 kHz constant
    const double referenceHz = -std::log (authored) * eon::detail::kAuthoringRate
                             / (2.0 * M_PI);

    EON_CHECK_NEAR (eon::detail::rescalePole (authored, eon::detail::kAuthoringRate),
                    authored, 0.0);

    for (const double rate : kRates)
    {
        const double pole = eon::detail::rescalePole (authored, rate);
        EON_CHECK (pole > authored && pole < 1.0);
        const double recoveredHz = -std::log (pole) * rate / (2.0 * M_PI);
        EON_CHECK_NEAR (recoveredHz, referenceHz, 1e-9 * referenceHz);
    }
}

EON_TEST_CASE("rescalePole rejects degenerate input instead of producing NaN")
{
    EON_CHECK_NEAR (eon::detail::rescalePole (1.0, 48000.0), 1.0, 0.0);
    EON_CHECK_NEAR (eon::detail::rescalePole (0.0, 48000.0), 0.0, 0.0);
    EON_CHECK_NEAR (eon::detail::rescalePole (-0.5, 48000.0), -0.5, 0.0);
    EON_CHECK_NEAR (eon::detail::rescalePole (0.9, 0.0), 0.9, 0.0);
    EON_CHECK (std::isfinite (eon::detail::rescalePole (0.9, -1.0)));
}

EON_TEST_CASE("unprepared stages keep their authoring-rate constants")
{
    const eon::ClassAStage classA;
    EON_CHECK_NEAR (classA.pole1, 0.9993, 0.0);
    EON_CHECK_NEAR (classA.pole2, 0.99996, 0.0);
    EON_CHECK_NEAR (classA.pole3, 0.999998, 0.0);
    EON_CHECK_NEAR (classA.sagTh.poleSag, 0.9970, 0.0);

    const eon::AnalogAir air;
    EON_CHECK_NEAR (air.poleLp, 0.90, 0.0);
    EON_CHECK_NEAR (air.polePink, 0.985, 0.0);
    EON_CHECK_NEAR (air.amplitude, 5e-5, 0.0);

    const eon::InductorResonator inductor;
    EON_CHECK_NEAR (inductor.pole1, 0.80, 0.0);
    EON_CHECK_NEAR (inductor.pole2, 0.86, 0.0);

    const eon::SagThermal sag;
    EON_CHECK_NEAR (sag.poleThermal, 0.99990, 0.0);
    EON_CHECK_NEAR (sag.poleDrift, 0.99998, 0.0);

    const eon::JilesAtherton transformer;
    EON_CHECK_NEAR (transformer.noiseScale, 1.0, 0.0);
}

EON_TEST_CASE("ClassA dielectric-absorption step response is rate independent")
{
    const auto setup = [] (eon::ClassAStage& stage, double rate) { stage.prepare (rate); };

    std::vector<std::vector<double>> runs;
    for (const double rate : kRates)
        runs.push_back (stepResponse<eon::ClassAStage> (rate, 0.5f, 0.5, 0.005, setup));

    EON_CHECK (relativeSpread (runs) < 0.02);
}

EON_TEST_CASE("Inductor tail step response is rate independent")
{
    const auto setup = [] (eon::InductorResonator& stage, double rate) { stage.prepare (rate); };
    const auto drive = [] (eon::InductorResonator& stage, float x)
    {
        return stage.process (x, 1600.f, 1.6f, 1.f);
    };

    // The stage takes extra per-sample arguments, so drive it directly rather
    // than through the generic helper.
    std::vector<std::vector<double>> runs;
    for (const double rate : kRates)
    {
        eon::InductorResonator stage;
        setup (stage, rate);
        stage.reset();

        std::vector<double> samples;
        const long long total = static_cast<long long> (0.2 * rate);
        const long long stride = static_cast<long long> (0.002 * rate);
        for (long long i = 0; i < total; ++i)
        {
            drive (stage, 0.25f);
            if (stride > 0 && i > 0 && (i % stride) == 0)
                samples.push_back (stage.lp1);
        }
        runs.push_back (samples);
    }

    EON_CHECK (relativeSpread (runs) < 0.02);
}

EON_TEST_CASE("PSU sag and thermal step response are rate independent")
{
    const auto setup = [] (eon::SagThermal& stage, double rate) { stage.prepare (rate); };

    std::vector<std::vector<double>> sagRuns, thermalRuns;
    for (const double rate : kRates)
    {
        eon::SagThermal stage;
        setup (stage, rate);
        stage.reset();

        std::vector<double> sagSamples, thermalSamples;
        const long long total = static_cast<long long> (1.0 * rate);
        const long long stride = static_cast<long long> (0.01 * rate);
        for (long long i = 0; i < total; ++i)
        {
            stage.tick (0.5);
            if (stride > 0 && i > 0 && (i % stride) == 0)
            {
                sagSamples.push_back (stage.sag);
                thermalSamples.push_back (stage.thermal);
            }
        }
        sagRuns.push_back (sagSamples);
        thermalRuns.push_back (thermalSamples);
    }

    EON_CHECK (relativeSpread (sagRuns) < 0.02);
    EON_CHECK (relativeSpread (thermalRuns) < 0.02);
}

EON_TEST_CASE("analog air keeps its in-band noise power across rates")
{
    // Total broadband RMS legitimately rises with rate (the noise is shaped and
    // there is more bandwidth), so compare the low-frequency power: a 1 kHz
    // one-pole smoother is applied to both runs before measuring.
    const auto inBandRms = [] (double rate)
    {
        eon::AnalogAir air;
        air.prepare (rate);
        air.reset();

        const double pole = std::exp (-2.0 * M_PI * 1000.0 / rate);
        double lp = 0.0, sumSq = 0.0;
        const long long total = static_cast<long long> (0.5 * rate);
        for (long long i = 0; i < total; ++i)
        {
            lp = lp * pole + static_cast<double> (air.process()) * (1.0 - pole);
            sumSq += lp * lp;
        }
        return std::sqrt (sumSq / static_cast<double> (total));
    };

    const double eco = inBandRms (kRates[0]);
    const double god = inBandRms (kRates[2]);

    EON_CHECK (eco > 0.0 && god > 0.0);
    const double ratio = god / eco;
    EON_CHECK (ratio > 0.8 && ratio < 1.25);
}

EON_TEST_CASE("analog air amplitude scales with sqrt of the rate")
{
    eon::AnalogAir eco; eco.prepare (kRates[0]);
    eon::AnalogAir god; god.prepare (kRates[2]);

    EON_CHECK_NEAR (god.amplitude / eco.amplitude, 2.0, 1e-9);
}

EON_TEST_CASE("DC blocker keeps its 20 Hz gain across rates")
{
    const auto gainAt20Hz = [] (double rate)
    {
        eon::DCBlocker blocker;
        blocker.prepare (rate, 18.0);
        blocker.reset();

        const double w = 2.0 * M_PI * 20.0 / rate;
        double peak = 0.0;
        const long long total = static_cast<long long> (2.0 * rate);
        const long long settle = total / 2;
        for (long long i = 0; i < total; ++i)
        {
            const float x = static_cast<float> (std::sin (w * static_cast<double> (i)));
            const double y = static_cast<double> (blocker.process (x));
            if (i >= settle)
                peak = std::max (peak, std::abs (y));
        }
        return peak;
    };

    const double eco = gainAt20Hz (kRates[0]);
    const double god = gainAt20Hz (kRates[2]);
    EON_CHECK (eco > 0.0);
    EON_CHECK_NEAR (god / eco, 1.0, 0.02);
}

EON_TEST_CASE("prepared stages stay finite on extreme input")
{
    for (const double rate : kRates)
    {
        eon::ClassAStage classA; classA.prepare (rate); classA.reset();
        eon::InductorResonator inductor; inductor.prepare (rate); inductor.reset();
        eon::AnalogAir air; air.prepare (rate); air.reset();
        eon::JilesAtherton transformer; transformer.prepare (rate); transformer.reset();

        for (int i = 0; i < 4096; ++i)
        {
            const float x = (i % 2 == 0) ? 40.f : -40.f;
            EON_CHECK (std::isfinite (classA.process (x)));
            EON_CHECK (std::isfinite (inductor.process (x, 3200.f, 5.f, 10.f)));
            EON_CHECK (std::isfinite (transformer.process (x)));
        }
        EON_CHECK (std::isfinite (air.process()));
    }
}
