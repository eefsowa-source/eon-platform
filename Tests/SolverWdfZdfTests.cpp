#include "TestHarness.h"

#include "Dsp/Solvers.h"
#include "Dsp/Wdf.h"
#include "Dsp/Zdf.h"

#include <cmath>

void runSolverWdfZdfTests()
{
    {
        const double w = eon::lambertW0log(-1000.0);
        eon::test::check(std::isfinite(w) && w >= 0.0,
                         "Lambert W log-domain solver handles underflow-scale input");
    }

    {
        eon::WdfISourceRes source(1000.0);
        source.Is = 0.001;
        source.incident(0.3);
        eon::test::check(eon::test::near(source.emitted(), -1.7, 1.0e-12),
                         "WDF current source emits its Norton wave");
    }

    {
        eon::WdfCapacitor capacitor(1.0e-6);
        capacitor.setSampleRate(48000.0);
        capacitor.incident(1.0);
        eon::test::check(eon::test::near(capacitor.voltage(), 0.5, 1.0e-12),
                         "WDF capacitor reports half the incident plus reflected wave");
    }

    {
        eon::WdfDiodePair diode;
        diode.R = 4700.0;
        for (double input : {1.0, 2.0, 20.0, 0.0})
        {
            diode.incident(input);
            eon::test::check(std::isfinite(diode.emitted()),
                             "WDF diode pair remains finite under large input and recovery");
        }
        diode.reset();
        diode.incident(0.0);
        eon::test::check(eon::test::near(diode.emitted(), 0.0, 1.0e-12),
                         "WDF diode reset clears the warm start");
    }

    {
        eon::Ladder4 ladder;
        ladder.setCutoff(12000.0, 48000.0);
        ladder.setResonance(4.0);
        ladder.reset();
        const double input = 2.0;
        const double output = ladder.process(input);
        const double G = ladder.g / (1.0 + ladder.g);
        const double G4 = G * G * G * G;
        const double u = output / G4;
        const double residual = u - std::tanh(input * ladder.drive - ladder.k * output);
        eon::test::check(std::abs(residual) < 1.0e-9,
                         "saturated ladder satisfies its implicit feedback equation");

        for (double cutoff : {20.0, 1000.0, 12000.0, 23000.0})
        {
            for (double resonance : {0.0, 2.0, 4.0})
            {
                eon::Ladder4 sweep;
                sweep.setCutoff(cutoff, 48000.0);
                sweep.setResonance(resonance);
                sweep.reset();
                for (int i = -400; i <= 400; ++i)
                {
                    const double value = sweep.process(static_cast<double>(i) / 100.0);
                    eon::test::check(std::isfinite(value) && std::abs(value) < 20.0,
                                     "saturated ladder remains finite and bounded");
                }
            }
        }
    }
}
