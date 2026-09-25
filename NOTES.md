# eon_dsp — shared pure-DSP core

Framework-free C++20 headers (no JUCE dependency) shared across EON plugin
projects. Include the repository root and include headers as `Dsp/<name>.h`.

| Header | Contents |
|---|---|
| `Dsp/Adaa.h` | Double-precision ADAA1/ADAA2 state and antiderivatives; `TanhSat`, `SoftClip`, and `SoftClipSat` |
| `Dsp/Measure.h` | Measurement and analysis helpers |
| `Dsp/Oversampling.h` | Stateful 2x/4x/8x half-band polyphase oversampler |
| `Dsp/Rate.h` | Rate compensation for per-sample one-pole constants (`rescalePole`, `rescaleNoise`) |
| `Dsp/Rng.h` | Deterministic per-instance PRNG |
| `Dsp/RtGuard.h` | Scoped denormal handling helper |
| `Dsp/Solvers.h` | Numerical solver helpers, including Lambert W |
| `Dsp/Stages.h` | DC blocker, inductor resonator, analog air, bounded OU sag/thermal, Class-A stage |
| `Dsp/Transformer.h` | Jiles-Atherton hysteresis transformer |
| `Dsp/Triode.h` | Koren triode with analytic plate-current derivative and warm-started Newton solve; optional Rk‖Ck cathode self-bias |
| `Dsp/Wdf.h` | Wave-digital ports, adaptors, reactive elements, diode and antiparallel diode-pair models |
| `Dsp/Zdf.h` | Zero-delay-feedback one-pole, state-variable, and ladder filters |

## Runtime contracts

- **Rate independence:** `InductorResonator`, `AnalogAir`, `SagThermal`,
  `ClassAStage`, and `JilesAtherton` store physical time constants and must be
  given `prepare(rateHz)` with the *processing* rate — the oversampled rate
  when the stage runs inside an oversampling loop. Without it they replay the
  48 kHz authoring constants per sample, so their poles (and the noise floor)
  scale with the oversampling factor. `Dsp/Rate.h` holds the transform;
  `Tests/RateIndependenceTests.cpp` is the gate.
- **DCBlocker:** Call `prepare(sampleRate, cutoffHz)` before processing. Both
  values must be finite and strictly positive; the default cutoff is 18 Hz.
  Call `reset()` when starting a new signal/state epoch.
- **Oversampler:** Call `prepare(maxBlock)` with the largest host-rate block
  that will be passed to `up()` or `down()`. Keep every call at or below that
  capacity; scratch storage is sized during prepare. Set the stage count (1–3,
  default 2/4x) before processing. Upsampled output needs `n * factor()` slots
  and must not alias the input. Upsampling and downsampling have separate
  histories; use `up()` for the upward path and `down()` for the return path.
  Histories persist across calls, so processing the same contiguous samples in
  different block partitions produces the same stream as one call. Call
  `reset()` when beginning a fresh stream.
- **ADAA:** Each `ADAA1`, `ADAA2`, `TanhSat`, or `SoftClipSat` stateful instance
  belongs to one channel and one nonlinear stage. Give each channel/stage its
  own object; sharing an instance interleaves histories.
- **TriodeStage:** Call `reset()` after parameter changes that affect its DC
  operating point, including circuit and cathode settings. Reset solves the
  quiescent point with a temporary higher iteration count, then restores the
  configured per-sample `iterations` budget. With cathode self-bias disabled,
  the grounded-cathode behavior is retained; do not describe reset as
  bit-identical to earlier implementations because the corrected operating
  point and solver behavior can change results.
- **WdfDiodePair:** `R`, `Is`, and `Vt` must each be finite and greater than
  zero for a valid solve. `iterations` is the total Newton/refinement budget;
  zero or less preserves the existing warm-start voltage (`vPrev`) and returns
  its reflected wave, `2 * vPrev - aIn`, without solving. Check
  `solveSucceeded` after `emitted()` to tell whether the solve converged or
  returned a finite fallback result. Invalid parameters and other fallback
  paths leave it false.
- **WdfCapacitor:** `voltage()` reports the voltage saved at the last
  scattering operation. It is a held state value and is not recomputed by the
  getter.

## Build and test

Run from the repository root. These commands build and execute the CTest suite
in separate Debug, Release, and AddressSanitizer/UndefinedBehaviorSanitizer
build directories:

```sh
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build-debug --parallel
ctest --test-dir build-debug --output-on-failure

cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure

cmake -S . -B build-sanitize -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"
cmake --build build-sanitize --parallel
ctest --test-dir build-sanitize --output-on-failure
```

## Verification record

2026-09-22, source revision `cb2cc2c`: 35 core tests passed; all 11 public
headers were included in header compile checks; Debug, Release, and
ASan/UBSan CTest runs each passed 2/2 tests. This records core-library build
and test evidence only. It does not establish plugin build/validation, host
loading, visual rendering, or DAW listening results.

## Vendored third-party (added 2026-09-18)

| Path | Contents | Include |
|---|---|---|
| `third_party/sst-basic-blocks` | Surge Synth Team DSP primitives (C++20, header-only): `dsp/` (filters, envelopes, clippers, quads), `modulators/`, `mod-matrix/`, `params/`, `simd/`, `tables/` | `-I third_party/sst-basic-blocks/include` then `#include "sst/basic-blocks/..."` |
| `third_party/sst-cpputils` | header-only utils required by sst-basic-blocks (`ring_buffer.h`, `rtsan_support.h`, allocators) | `-I third_party/sst-cpputils/include` |
| `third_party/simde` | SIMDe — SSE intrinsics shim required by sst-basic-blocks on Apple Silicon | `-I third_party/simde` |

Update with `git -C third_party/<repo> pull` (shallow clones).
All three include paths are needed together on arm64.
`mod-matrix` needs sst-cpputils; without its include path define
`SST_CPPUTILS_UNAVAILABLE=1` to drop those features.

## Duplicate implementations elsewhere (surveyed 2026-09-18)

| Project | File | vs eon_dsp | Verdict |
|---|---|---|---|
| Tube comp | `DSP/AADAWaveshaper.h` | ADAA1-only, **float32 state** — the exact failure mode `Dsp/Adaa.h` documents (F-difference quantization noise) | eon_dsp wins; Tube comp should migrate to `Dsp/Adaa.h` |
| Tube comp | `DSP/KorenTriode.h` + `TriodeStage.h` | Same Koren equations but **numeric (central-diff) derivative**, JUCE `ProcessContext`-coupled | Cathode-sag mechanism was ported 2026-09-18 into `Dsp/Triode.h` as opt-in `setCathode()`; compare current behavior and reset semantics before migration. Not yet ported: grid/output RC coupling networks and bias/harmonic-ratio params |
| EEVE1073 | `Shared/EeFourTimesOversampler.h` | JUCE oversampler wrapper, not pure DSP | Stays per-project (or becomes the JUCE-binding layer later) |

Rule for adopters: migrate to these headers only when the project's own
test gate (pluginval / smoke tests / sound review) can verify the swap.
