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
- **Oversampler eco stage 0:** `ecoStage0` selects a separate filter with its
  own delay history and latency. Change it only while the stream is stopped,
  then call `reset()`; flipping it mid-stream clicks (measured: ~35x the normal
  sample-to-sample step on a 200 Hz sine).
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
- **WdfDiode / WdfVSourceIdeal / WdfISourceRes:** `reset()` clears the stored
  incident wave. `WdfDiode` returns `-a` (total reflection) when `R`, `Is`,
  or `Vt` is invalid and 0 for non-finite input. `WdfISourceRes` evaluates
  `b = a - 2R*Is` with the last incident wave it saw — exact at the root,
  one sample late as a leaf.
- **Rng:** instances sharing the default seed emit identical streams; call
  `reseed()` per channel/stage that needs decorrelated noise.

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

2026-10-23, working tree on top of `d851186`: same 51 core tests pass in
Release + ASan/UBSan and in Debug, all 12 headers compile standalone. This
revision halves the oversampler FIR work: `HalfBand2x` evaluates both polyphase
sums in symmetric-folded form (every tap set is `static constexpr` and now
gated by `hbTapsSymmetric`), replaces the `%` delay-line indexing with a
compare-based wrap and the even-branch line with masked indexing. Measured
up+down roundtrip on one machine: 2x 46.1 -> 26.4 ns, 4x 110.1 -> 61.4 ns,
8x 206.3 -> 118.6 ns per base sample; output is bit-identical to the previous
header across 1,200,000 samples (all stage counts and both `ecoStage0`
settings), so block-partition invariance and `latencySamples()` are unchanged.
`.github/workflows/bump-products.yml` now rejects non-release tags, fails on
`gh` auth/API errors instead of reporting them as "no platform pin", and only
prints "bump PR opened" after a successful `gh pr create`. The 2026-10-22
record below is kept for history.

2026-10-22, working tree on top of `d851186`: 51 core tests pass in Release
+ ASan/UBSan and in Debug, all 12 headers compile standalone. Fixes in this
revision: `InductorResonator` output tap no longer scales with the processing
rate (~1.1 dB gain drift at 8x removed); `ADAA2` uses the three-node centroid
when the node span is < 1e-4 (slow-signal error 6e-3 -> ~5e-8); `DCBlocker::prepare`
rejects +-Inf. Each has a regression test that fails on the previous code.
The 2026-09-22 record below is kept for history.

2026-09-22, source revision `cb2cc2c`: 35 core tests passed; all 11 public
headers were included in header compile checks; Debug, Release, and
ASan/UBSan CTest runs each passed 2/2 tests. This records core-library build
and test evidence only. It does not establish plugin build/validation, host
loading, visual rendering, or DAW listening results.

## Vendored third-party (added 2026-09-18)

**Not part of the library, and not linked by anything.** No header in `Dsp/`,
no test, no CMake target and no CI workflow references `third_party/`, so
nothing here is compiled, linked or shipped by `eon_dsp`.

| Path | Contents | License | State |
|---|---|---|---|
| `third_party/sst-basic-blocks` | Surge Synth Team DSP primitives (header-only: filters, envelopes, clippers, quads, modulators, mod-matrix, params, simd, tables) | **GPL-3.0** | **Removed 2026-10-24** — never used here and incompatible with an MIT build |
| `third_party/sst-cpputils` | header-only utils that `sst-basic-blocks` needed (`ring_buffer.h`, `rtsan_support.h`, allocators) | MIT | Present; unused |
| `third_party/simde` | SIMDe — SSE intrinsics shim that `sst-basic-blocks` needed on Apple Silicon | MIT | Present; unused |

Update with `git -C third_party/<repo> pull` (shallow clones).

### Why `sst-basic-blocks` was removed

It is GPL-3.0 and `eon_dsp` is MIT. Adding its include path to a product
build, or shipping that code inside a plugin binary, makes the combined work
GPL-3.0. Nothing in this repository ever included it, so it was a liability with
no benefit — the gitlink and its `.gitmodules` stanza were dropped in
`third_party/sst-basic-blocks` removal rather than left for someone to
discover later.

`sst-cpputils` and `simde` are MIT and are kept as reference, but they only
existed to serve `sst-basic-blocks`, so they are equally unused here and are
candidates for the same treatment. To adopt a primitive from any of this code,
reimplement it from the published DSP math and gate the result on this
library's own tests.

The history is recoverable if ever needed:
`git show cac8aaf:.gitmodules` and `git ls-tree cac8aaf third_party/`.

## Duplicate implementations elsewhere (surveyed 2026-09-18)

| Project | File | vs eon_dsp | Verdict |
|---|---|---|---|
| Tube comp | `DSP/AADAWaveshaper.h` | ADAA1-only, **float32 state** — the exact failure mode `Dsp/Adaa.h` documents (F-difference quantization noise) | eon_dsp wins; Tube comp should migrate to `Dsp/Adaa.h` |
| Tube comp | `DSP/KorenTriode.h` + `TriodeStage.h` | Same Koren equations but **numeric (central-diff) derivative**, JUCE `ProcessContext`-coupled | Cathode-sag mechanism was ported 2026-09-18 into `Dsp/Triode.h` as opt-in `setCathode()`; compare current behavior and reset semantics before migration. Not yet ported: grid/output RC coupling networks and bias/harmonic-ratio params |
| EEVE1073 | `Shared/EeFourTimesOversampler.h` | JUCE oversampler wrapper, not pure DSP | Stays per-project (or becomes the JUCE-binding layer later) |

Rule for adopters: migrate to these headers only when the project's own
test gate (pluginval / smoke tests / sound review) can verify the swap.
