# eon_dsp — shared pure-DSP core

Framework-free headers (no JUCE dependency) shared across EON plugin
projects. Include as `-I <repo>/eon_dsp`, headers under `Dsp/`.

| Header | Contents |
|---|---|
| `Dsp/Adaa.h` | ADAA1/ADAA2 (double-precision state + antiderivatives), `TanhSat`, `SoftClipSat` |
| `Dsp/Triode.h` | Koren triode + Newton-Raphson plate solve (analytic dIp/dVp, warm start) + optional Rk‖Ck cathode self-bias sag (off = bit-identical; ported from Tube comp's mechanism, analytic derivative replaces its numeric one) |
| `Dsp/Transformer.h` | Jiles-Atherton hysteresis transformer |
| `Dsp/Stages.h` | DC blocker, inductor resonator, analog air, bounded OU sag/thermal, Class-A stage |
| `Dsp/Rng.h` | deterministic per-instance PRNG |

## Runtime contracts

- Call `DCBlocker::prepare(processingRateHz, cutoffHz)` before audio processing.
- Call `Oversampler::prepare(maxBlockSize)` before `up()`/`down()`; output is
  invariant to host block partition.
- Each ADAA object belongs to exactly one logical processing stage and channel.
- `reset()` establishes a finite, zero-input operating state for triodes, WDF
  roots, filters, and oversamplers.

## Test

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

JUCE-dependent pieces (oversampling manager, LookAndFeel) deliberately stay
in each project's own `Source/Shared/`.

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
| Tube comp | `DSP/KorenTriode.h` + `TriodeStage.h` | Same Koren equations but **numeric (central-diff) derivative**, JUCE `ProcessContext`-coupled | **Cathode-sag ported 2026-09-18** into `Dsp/Triode.h` as opt-in `setCathode()` (off = bit-identical). Not yet ported: grid/output RC coupling networks, bias/harmonic-ratio params — port when a comp product lands |
| EEVE1073 | `Shared/EeFourTimesOversampler.h` | JUCE oversampler wrapper, not pure DSP | Stays per-project (or becomes the JUCE-binding layer later) |

Rule for adopters: migrate to these headers only when the project's own
test gate (pluginval / smoke tests / sound review) can verify the swap.
