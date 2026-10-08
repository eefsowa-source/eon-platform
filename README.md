# eon-platform

Shared pure-DSP core for EON audio plugins — header-only C++20, framework-free
(no JUCE dependency). Every public header under `Dsp/` compiles standalone and
is verified by the CTest suite in `Tests/`.

## Features

- **Antialiased nonlinearities** — ADAA1/ADAA2 in double precision with
  ready-made saturators (`TanhSat`, `SoftClipSat`)
- **Oversampling** — stateful 2x/4x/8x half-band polyphase oversampler with
  block-partition invariance (output is identical no matter how the host
  slices its buffers)
- **Rate-independent analog stages** — DC blocker, inductor resonator,
  analog air noise, PSU sag/thermal drift, Class-A stage, and a Jiles–Atherton
  hysteresis transformer, all holding physical time constants under any
  processing rate
- **Circuit solvers** — Koren triode with analytic plate-current derivative
  and optional cathode self-bias, wave-digital ports/adaptors/elements
  including diode and antiparallel diode-pair models, Lambert W helper
- **Zero-delay-feedback filters** — TPT one-pole, state-variable filter, and
  a bracketed four-pole ladder with tanh saturation

## Requirements

- C++20 compiler (tested with GCC, Clang/AppleClang, MSVC)
- CMake ≥ 3.16 — only needed to build the test suite; the library itself is
  just headers

## Getting started

The repository root is the include root:

```cpp
#include "Dsp/Stages.h"      // eon::DCBlocker
#include "Dsp/Triode.h"      // eon::TriodeStage
```

With CMake, either add this repo as a subdirectory/FetchContent and link the
interface target, or add it to your include path directly:

```cmake
add_subdirectory(eon-platform)      # or FetchContent_MakeAvailable(eon-platform)
target_link_libraries(your_target PRIVATE eon_dsp)
```

Consumers never get the platform's own test suite — tests only build when this
repo is the top-level project.

## Quick examples

### DC blocker

```cpp
#include "Dsp/Stages.h"

eon::DCBlocker dc;
dc.prepare (48000.0, 20.0);   // sample rate, cutoff [Hz]
dc.reset();

float y = dc.process (x);     // per sample
```

### Ladder filter

```cpp
#include "Dsp/Zdf.h"

eon::Ladder4 ladder;
ladder.setCutoff (1200.0, 48000.0);  // fc, sample rate
ladder.setResonance (2.0);           // 0..4
ladder.drive = 1.5;                  // pre-gain into the tanh feedback
ladder.reset();

double y = ladder.process (x);
```

### Triode stage under oversampling

```cpp
#include "Dsp/Oversampling.h"
#include "Dsp/Triode.h"

eon::Oversampler os;
os.setStages (2);              // 4x
os.prepare (maxBlockSize);     // largest host-rate block you will pass
os.reset();

eon::TriodeStage triode;
triode.setCathode (2700.f, 25e-6f, 48000.0 * os.factor());  // *processing* rate
triode.reset();                // solves the quiescent operating point

// per block:
os.up (input, n, upsampled);   // upsampled has n * os.factor() slots
for (int i = 0; i < n * os.factor(); ++i)
    upsampled[i] = triode.process (upsampled[i]);
os.down (upsampled, output, n);
```

## Module overview

| Header | Public types | Contract highlights |
|---|---|---|
| `Dsp/Adaa.h` | `ADAA1`, `ADAA2`, `TanhSat`, `SoftClip`, `SoftClipSat` | One instance per channel per nonlinear stage; `reset()` clears antiderivative state |
| `Dsp/Measure.h` | measurement helpers | analysis utilities |
| `Dsp/Oversampling.h` | `Oversampler`, `HalfBand2x`, `HbTaps*` | `prepare(maxBlock)` before use; `setStages(1..3)`; separate up/down histories; output is block-partition invariant |
| `Dsp/Rate.h` | `rescalePole`, `rescaleNoise` | transforms one-pole constants and noise amplitudes between rates |
| `Dsp/Rng.h` | `Rng` | deterministic per-instance PRNG |
| `Dsp/RtGuard.h` | `ScopedDenormalsOff` | scoped FTZ/DAZ for audio callbacks |
| `Dsp/Solvers.h` | solver helpers incl. Lambert W | used by Wdf/Triode internals |
| `Dsp/Stages.h` | `DCBlocker`, `InductorResonator`, `AnalogAir`, `SagThermal`, `ClassAStage` | call `prepare(rateHz)` with the *processing* rate (oversampled rate inside an oversampling loop) — otherwise the stage replays its 48 kHz authoring constants |
| `Dsp/Transformer.h` | `JilesAtherton` | hysteresis transformer; `prepare(rateHz)` rate-compensates Barkhausen noise |
| `Dsp/Triode.h` | `TriodeStage` | Koren model, warm-started Newton; `setCathode()` then `reset()` for Rk‖Ck self-bias; `reset()` again after changing anything that moves the DC operating point |
| `Dsp/Wdf.h` | `WdfPort`, `WdfResistor/Capacitor/Inductor`, sources, `WdfSeries2/3`, `WdfParallel2`, `WdfDiode`, `WdfDiodePair` | `WdfDiodePair`: `R`, `Is`, `Vt` must be finite and > 0; check `solveSucceeded` after `emitted()`; `iterations <= 0` falls back to the warm-start voltage |
| `Dsp/Zdf.h` | `OnePoleTPT`, `SvfTPT`, `Ladder4` | TPT topologies; `Ladder4` solves tanh in the feedback loop with bracketing |

Full runtime contracts (call ordering, state ownership, fallback behavior) live
in [NOTES.md](NOTES.md).

## Building and testing

From the repository root:

```sh
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build-debug --parallel
ctest --test-dir build-debug --output-on-failure
```

Release and sanitizer configurations:

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release && cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure

cmake -S . -B build-sanitize -DEON_DSP_SANITIZERS=ON && cmake --build build-sanitize --parallel
ctest --test-dir build-sanitize --output-on-failure
```

The suite also builds `eon_dsp_header_checks`, which compiles every public
header in isolation.

**Never use `-ffast-math` with this library.** The solvers rely on
`std::isfinite()` and exact NaN/Inf semantics to stay finite; fast-math
breaks that contract. Tests are always compiled with strict floating point
(`-ffp-contract=off -fno-fast-math`, or `/fp:strict` on MSVC).

## Troubleshooting

- **Stage sounds wrong under oversampling** — you probably passed the host
  rate to `prepare()`. Pass the *processing* rate (host rate × oversampling
  factor).
- **Triode sounds different after parameter change** — call `reset()`; the
  quiescent operating point is only solved there.
- **Shared state glitches across channels** — every channel needs its own
  instance of stateful objects (ADAA states, filters, stages).
- **Assert in `Oversampler::up/down`** — the block is larger than
  `prepare(maxBlock)` capacity.
- **Diode-pair output looks flat** — check `solveSucceeded`; invalid
  `R`/`Is`/`Vt` or a non-converged solve falls back deterministically.

## License

MIT — see [LICENSE](LICENSE). SPDX headers mark every public header.

Vendored dependencies under `third_party/` (sst-basic-blocks, sst-cpputils,
simde) keep their own licenses.
