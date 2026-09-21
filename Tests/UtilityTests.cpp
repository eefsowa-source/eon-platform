#include "Dsp/Rng.h"
#include "Dsp/Solvers.h"
#include "Dsp/Wdf.h"
#include "Tests/TestHarness.h"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>

EON_TEST_CASE("utility suite smoke test")
{
    eon::Rng first(1234);
    eon::Rng second(1234);
    EON_CHECK(first.nextU64() == second.nextU64());
    EON_CHECK_NEAR(first.next(), second.next(), 0.0);
    EON_CHECK_NEAR(0.1 + 0.2, 0.3, 1e-15);
}

EON_TEST_CASE("Lambert W log solver stays in range from tiny to huge inputs")
{
    const double underflow = eon::lambertW0log(-1000.0);
    EON_CHECK(std::isfinite(underflow));
    EON_CHECK(underflow >= 0.0);

    EON_CHECK_NEAR(eon::lambertW0(1.0), 0.56714329040978387, 1e-14);
    EON_CHECK_NEAR(eon::lambertW0(1e300) + std::log(eon::lambertW0(1e300)),
                   std::log(1e300), 1e-12);
    for (const double logX : {-10.0, 0.0, 1.0, 100.0, 1000.0})
    {
        const double w = eon::lambertW0log(logX);
        EON_CHECK(std::isfinite(w));
        EON_CHECK(w > 0.0);
        EON_CHECK_NEAR(w + std::log(w), logX, 1e-12 * std::max(1.0, std::abs(logX)));
    }

    EON_CHECK(eon::lambertW0log(-std::numeric_limits<double>::infinity()) == 0.0);
    EON_CHECK(eon::lambertW0log(std::numeric_limits<double>::infinity())
              == std::numeric_limits<double>::infinity());
    EON_CHECK(std::isnan(eon::lambertW0log(std::numeric_limits<double>::quiet_NaN())));
}

EON_TEST_CASE("resistive current source reflects its incident wave")
{
    eon::WdfISourceRes source(1000.0);
    source.Is = 0.001;
    source.incident(0.3);
    EON_CHECK_NEAR(source.emitted(), -1.7, 1e-15);

    source.reset();
    EON_CHECK_NEAR(source.emitted(), -2.0, 1e-15);
}

EON_TEST_CASE("capacitor voltage uses latest incident and prior reflected waves")
{
    eon::WdfCapacitor capacitor(1e-6);
    capacitor.setSampleRate(48000.0);

    EON_CHECK_NEAR(capacitor.emitted(), 0.0, 0.0);
    EON_CHECK_NEAR(capacitor.voltage(), 0.0, 0.0);
    capacitor.incident(1.0);
    EON_CHECK_NEAR(capacitor.voltage(), 0.5, 0.0);
    EON_CHECK_NEAR(capacitor.emitted(), 1.0, 0.0);

    capacitor.incident(0.4);
    EON_CHECK_NEAR(capacitor.voltage(), 0.7, 1e-15);
    EON_CHECK_NEAR(capacitor.emitted(), 0.4, 0.0);

    capacitor.reset();
    EON_CHECK_NEAR(capacitor.voltage(), 0.0, 0.0);
    EON_CHECK_NEAR(capacitor.emitted(), 0.0, 0.0);
}

EON_TEST_CASE("diode pair remains finite under large and changing drives")
{
    eon::WdfDiodePair diodePair;
    diodePair.R = 4700.0;
    for (const double input : {1.0, 2.0, 20.0, 0.0})
    {
        diodePair.incident(input);
        const double reflected = diodePair.emitted();
        EON_CHECK(std::isfinite(reflected));
        EON_CHECK(std::isfinite(diodePair.vPrev));
    }

    diodePair.reset();
    diodePair.incident(0.0);
    EON_CHECK_NEAR(diodePair.emitted(), 0.0, 0.0);
    EON_CHECK_NEAR(diodePair.vPrev, 0.0, 0.0);

    diodePair.incident(1.0);
    const double positive = diodePair.emitted();
    diodePair.reset();
    diodePair.incident(-1.0);
    EON_CHECK_NEAR(diodePair.emitted(), -positive, 1e-10);
}
