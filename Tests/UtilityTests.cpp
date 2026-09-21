#include "Dsp/Rng.h"
#include "Dsp/Solvers.h"
#include "Dsp/Wdf.h"
#include "Dsp/Zdf.h"
#include "Tests/TestHarness.h"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>

template <typename Solver>
static bool reportsSolveFailure(const Solver& solver)
{
    if constexpr (requires { solver.solveSucceeded; })
        return ! solver.solveSucceeded;
    return false;
}

template <typename Solver>
static bool reportsSolveSuccess(const Solver& solver)
{
    if constexpr (requires { solver.solveSucceeded; })
        return solver.solveSucceeded;
    return false;
}

static double ladderFeedbackStateTerm(double G, const double states[4])
{
    const double stateSum = G * G * G * states[0]
                          + G * G * states[1]
                          + G * states[2]
                          + states[3];
    return (1.0 - G) * stateSum;
}

EON_TEST_CASE("utility suite smoke test")
{
    eon::Rng first(1234);
    eon::Rng second(1234);
    EON_CHECK(first.nextU64() == second.nextU64());
    EON_CHECK_NEAR(first.next(), second.next(), 0.0);
    EON_CHECK_NEAR(0.1 + 0.2, 0.3, 1e-15);
}

EON_TEST_CASE("four-pole ladder solves tanh saturation inside feedback loop")
{
    eon::Ladder4 ladder;
    ladder.setCutoff(12000.0, 48000.0);
    ladder.setResonance(4.0);
    ladder.reset();

    constexpr double input = 2.0;
    const double G = ladder.g / (1.0 + ladder.g);
    const double G4 = G * G * G * G;
    const double stateTerm = ladderFeedbackStateTerm(G, ladder.s);
    const double output = ladder.process(input);
    const double u = (output - stateTerm) / G4;
    const double residual = u - std::tanh(input * ladder.drive
                                        - ladder.k * (G4 * u + stateTerm));

    EON_CHECK_NEAR(residual, 0.0, 1e-9);
}

EON_TEST_CASE("four-pole ladder feedback residual holds across stateful settings")
{
    constexpr double sampleRate = 48000.0;
    const double cutoffs[] = {20.0, 1000.0, 12000.0, 23000.0};
    const double resonances[] = {0.0, 2.0, 4.0};
    const double drives[] = {-4.0, 4.0};
    const double inputs[] = {-0.75, 0.2, 0.9, -0.4, 0.1};

    for (const double cutoff : cutoffs)
    {
        for (const double resonance : resonances)
        {
            for (const double drive : drives)
            {
                eon::Ladder4 ladder;
                ladder.setCutoff(cutoff, sampleRate);
                ladder.setResonance(resonance);
                ladder.drive = drive;
                ladder.reset();

                const double G = ladder.g / (1.0 + ladder.g);
                const double G4 = G * G * G * G;
                for (const double input : inputs)
                {
                    double preStates[4];
                    for (int i = 0; i < 4; ++i)
                        preStates[i] = ladder.s[i];
                    const double stateTerm = ladderFeedbackStateTerm(G, preStates);
                    const double output = ladder.process(input);
                    const double u = (output - stateTerm) / G4;
                    const double residual = u - std::tanh(input * ladder.drive
                                                        - ladder.k * (G4 * u + stateTerm));
                    const double conditioningTolerance = 64.0
                        * std::numeric_limits<double>::epsilon()
                        * (std::abs(output) + std::abs(stateTerm)) / G4;

                    EON_CHECK(std::isfinite(output));
                    EON_CHECK(std::abs(output) <= 2.0);
                    EON_CHECK(std::abs(residual) <= std::max(1e-9, conditioningTolerance));
                    for (const double state : ladder.s)
                    {
                        EON_CHECK(std::isfinite(state));
                        EON_CHECK(std::abs(state) <= 4.0);
                    }
                }
            }
        }
    }
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

EON_TEST_CASE("diode pair keeps extreme finite drives inside the diode voltage range")
{
    eon::WdfDiodePair diodePair;
    diodePair.R = 4700.0;
    const double maxDiodeVoltage = diodePair.Vt
                                 * (std::log(std::numeric_limits<double>::max())
                                    - std::log(diodePair.Is));

    for (const double input : {1e100, -1e100, 1e300, -1e300})
    {
        diodePair.incident(input);
        const double reflected = diodePair.emitted();
        EON_CHECK(std::isfinite(reflected));
        EON_CHECK(input > 0.0 ? reflected < 0.0 : reflected > 0.0);
        EON_CHECK(std::abs(diodePair.vPrev) <= maxDiodeVoltage);
        EON_CHECK_NEAR(reflected / input, -1.0, 1e-12);
    }

    diodePair.incident(1.0);
    EON_CHECK(std::isfinite(diodePair.emitted()));
    EON_CHECK(diodePair.emitted() < 0.0);
    diodePair.reset();
    EON_CHECK_NEAR(diodePair.vPrev, 0.0, 0.0);
    diodePair.incident(0.0);
    EON_CHECK_NEAR(diodePair.emitted(), 0.0, 0.0);
}

EON_TEST_CASE("diode pair default budget accepts representable high-drive roots")
{
    eon::WdfDiodePair diodePair;
    diodePair.R = 4700.0;
    for (const double input : {1e48, 1e90, -1e90, -1e48})
    {
        diodePair.incident(input);
        const double reflected = diodePair.emitted();
        const double q = diodePair.vPrev / diodePair.Vt;
        const double current = std::copysign(
            std::exp(std::log(diodePair.Is) + std::abs(q)), q);
        const double voltageResidual = (diodePair.vPrev - input)
                                     + diodePair.R * current;

        EON_CHECK(reportsSolveSuccess(diodePair));
        EON_CHECK(input > 0.0 ? reflected < 0.0 : reflected > 0.0);
        EON_CHECK(std::abs(voltageResidual) <= 1e-12 * std::abs(input));
    }
}

EON_TEST_CASE("diode pair reports invalid input and recovers solve status")
{
    eon::WdfDiodePair diodePair;
    diodePair.R = 0.0;
    diodePair.incident(1e100);
    (void) diodePair.emitted();
    EON_CHECK(reportsSolveFailure(diodePair));

    diodePair.R = 4700.0;
    diodePair.Is = std::numeric_limits<double>::quiet_NaN();
    diodePair.incident(-1e100);
    (void) diodePair.emitted();
    EON_CHECK(reportsSolveFailure(diodePair));

    diodePair.Is = 2.52e-9;
    diodePair.incident(std::numeric_limits<double>::infinity());
    (void) diodePair.emitted();
    EON_CHECK(reportsSolveFailure(diodePair));

    diodePair.incident(1.0);
    (void) diodePair.emitted();
    EON_CHECK(reportsSolveSuccess(diodePair));
}

EON_TEST_CASE("diode pair voltage residual remains valid at tiny positive resistances")
{
    eon::WdfDiodePair diodePair;
    for (const double resistance : {1e-300, std::numeric_limits<double>::denorm_min()})
    {
        diodePair.R = resistance;
        diodePair.reset();
        diodePair.incident(1.0);
        const double reflected = diodePair.emitted();
        const double voltage = diodePair.vPrev;
        const double scaledResidual = voltage - 1.0
                                    + resistance * (2.0 * diodePair.Is
                                                  * std::sinh(voltage / diodePair.Vt));
        EON_CHECK_NEAR(voltage, 1.0, 1e-12);
        EON_CHECK_NEAR(reflected, 1.0, 1e-14);
        EON_CHECK_NEAR(scaledResidual, 0.0, 1e-14);
        EON_CHECK(reportsSolveSuccess(diodePair));
    }
}

EON_TEST_CASE("diode pair brackets high drive using scaled current at tiny resistance")
{
    eon::WdfDiodePair diodePair;
    for (const double resistance : {1e-300, std::numeric_limits<double>::denorm_min()})
    {
        diodePair.R = resistance;
        diodePair.reset();
        diodePair.incident(20.0);
        const double reflected = diodePair.emitted();
        const double voltage = diodePair.vPrev;
        const double q = voltage / diodePair.Vt;
        const double scaledCurrent = std::abs(q) < 20.0
            ? resistance * 2.0 * diodePair.Is * std::sinh(q)
            : std::copysign(std::exp(std::log(resistance) + std::log(diodePair.Is)
                                     + std::abs(q)), q);
        const double scaledResidual = (voltage - 20.0) + scaledCurrent;

        EON_CHECK(reportsSolveSuccess(diodePair));
        EON_CHECK(reflected > 0.0);
        EON_CHECK(voltage > 0.0 && voltage <= 20.0);
        EON_CHECK_NEAR(scaledResidual, 0.0, 1e-12);
    }
}

EON_TEST_CASE("diode pair tiny signals retain their small-signal response")
{
    eon::WdfDiodePair diodePair;
    diodePair.R = 4700.0;
    constexpr double input = 1e-21;
    diodePair.incident(input);
    const double reflected = diodePair.emitted();
    const double expectedVoltage = input / (1.0 + diodePair.R * 2.0 * diodePair.Is / diodePair.Vt);
    const double expectedReflected = 2.0 * expectedVoltage - input;

    EON_CHECK_NEAR(diodePair.vPrev, expectedVoltage, 1e-35);
    EON_CHECK_NEAR(reflected, expectedReflected, 1e-35);
    EON_CHECK(reportsSolveSuccess(diodePair));
}

EON_TEST_CASE("diode pair iterations bound the total solve work")
{
    eon::WdfDiodePair diodePair;
    diodePair.R = 4700.0;
    diodePair.incident(1.0);
    (void) diodePair.emitted();
    const double warmStart = diodePair.vPrev;

    diodePair.iterations = 0;
    diodePair.incident(2.0);
    const double zeroBudgetOutput = diodePair.emitted();
    EON_CHECK_NEAR(diodePair.vPrev, warmStart, 0.0);
    EON_CHECK_NEAR(zeroBudgetOutput, 2.0 * warmStart - 2.0, 0.0);
    EON_CHECK(reportsSolveFailure(diodePair));

    diodePair.iterations = 1;
    diodePair.incident(2.0);
    const double oneStepOutput = diodePair.emitted();
    const double oneStepVoltage = diodePair.vPrev;
    EON_CHECK(oneStepVoltage != warmStart);

    diodePair.iterations = 8;
    diodePair.incident(2.0);
    const double defaultBudgetOutput = diodePair.emitted();
    EON_CHECK(std::abs(defaultBudgetOutput - oneStepOutput) > 1e-10);
    EON_CHECK(reportsSolveSuccess(diodePair));
}
