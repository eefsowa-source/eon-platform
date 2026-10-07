# eon-platform GitHub Issue Checklist

## 우선순위 1: 필수 프로젝트 기초

### Issue #1: License 추가
**Priority**: CRITICAL  
**Type**: Setup  
**Effort**: 5 min

**문제**: 공개 repo이지만 라이선스가 없어서 재사용성이 제약됨

**할 일**:
- [ ] 라이선스 선택 (MIT / Apache 2.0 / GPL 추천)
- [ ] LICENSE 파일 생성
- [ ] 모든 헤더 파일에 라이선스 header 추가
- [ ] README에 라이선스 명시

**참고**: 공개 repo이면서 라이선스가 없으면 법적으로 "사용 권한 불명"이 됨

---

### Issue #2: README.md 작성
**Priority**: CRITICAL  
**Type**: Documentation  
**Effort**: 2-3 hours

**문제**: README가 없어서 프로젝트 목적과 사용법을 알 수 없음

**할 일**:
- [ ] 프로젝트 개요 작성
  - eon-platform이 무엇인가 (shared DSP core for EON plugins)
  - framework-free, header-only 특징
  - C++20 requirement
- [ ] 빠른 시작 가이드
  - CMake 빌드 방법
  - `prepare()` / `reset()` / `process()` 패턴
  - rate compensation 개념 설명
- [ ] API 문서
  - Dsp/의 각 모듈별 간단한 설명
  - public struct/function 목록
  - 각 모듈의 contract (prepare가 필요한지, per-channel stateful 등)
- [ ] 예제 코드
  - "simple lowpass filter"
  - "triode stage with oversampling"
  - "complete signal chain" example
- [ ] Troubleshooting section
  - 흔한 오류 (prepare() 호출 안 함, rate 잘못 넘김 등)

**참고**: NOTES.md / AGENTS.md는 이미 있으므로, README는 초보자 친화적으로

---

### Issue #3: .gitignore 및 프로젝트 파일 정리
**Priority**: MEDIUM  
**Type**: Setup  
**Effort**: 30 min

**문제**: .gitignore이 너무 작고, CMake 빌드 산물들이 commit되지 않음 확인 필요

**할 일**:
- [ ] .gitignore 확장
  ```
  build*/
  *.o
  *.a
  *.so
  *.dylib
  *.dll
  .DS_Store
  *.swp
  *.swo
  *~
  CMakeCache.txt
  CMakeFiles/
  Makefile
  ```
- [ ] .editorconfig 추가 (코드 스타일 일관성)
- [ ] .github/workflows 구조 확인

---

## 우선순위 2: 테스트 강화

### Issue #4: WdfDiodePair extreme-value tests
**Priority**: HIGH  
**Type**: Testing  
**Effort**: 2-3 hours

**문제**: WdfDiodePair는 가장 복잡한 solver인데, edge case 테스트 부족

**할 일**:
- [ ] Tests/WdfEdgeCaseTests.cpp 생성
- [ ] 다음 케이스 추가:
  - `aIn = 0` (should return 0)
  - `aIn = very large` (e.g., 1e6)
  - `aIn = very small` (e.g., 1e-12)
  - `R = very small` (e.g., 1e-6)
  - `R = very large` (e.g., 1e6)
  - `Is = boundary values`
  - `Vt = boundary values`
  - `iterations = 0` fallback (should use vPrev)
  - `vPrev` warm-start convergence verification
  - NaN/Inf input should not crash
- [ ] CMakeLists.txt에 테스트 추가
- [ ] 모든 테스트 통과 확인

**참고**: WdfDiodePair::emitted()는 270줄 이상의 복잡한 코드이므로, 테스트가 특히 중요함

---

### Issue #5: TriodeStage state management tests
**Priority**: HIGH  
**Type**: Testing  
**Effort**: 2 hours

**문제**: cathode self-bias와 plate solve의 coupling이 복잡해서, state bug 가능성 높음

**할 일**:
- [ ] Tests/TriodeStateTests.cpp 생성
- [ ] 다음 시나리오 테스트:
  - `reset()` 전후 `Vq` 변화 검증
  - `setCathode(0, 0, rate)` then `reset()` (grounded cathode)
  - `setCathode(1000, 25e-6, 48000)` then `reset()` (self-biased)
  - process() 호출 후 `VpPrev` warm-start 효과 확인
  - 다양한 `Bplus`, `Rload`, `biasVg` 조합에서 plate solve convergence
  - cathode voltage가 bounded range를 넘지 않는지 확인
  - quiescent operating point가 반복 가능한지 확인 (deterministic)
- [ ] CMakeLists.txt에 테스트 추가

**참고**: Triode.h의 reset()은 bisection + solveVp 조합으로 복잡하므로, 단계별 검증 필수

---

### Issue #6: Oversampler roundtrip & alias tests
**Priority**: HIGH  
**Type**: Testing  
**Effort**: 2-3 hours

**문제**: Oversampler가 제대로 alias를 억제하고, latency를 맞추는지 검증 필요

**할 일**:
- [ ] Tests/OversamplerAdvancedTests.cpp 생성
- [ ] 다음 검사 추가:
  - **Roundtrip exactness**: up() → down() 후 결과가 원본과 비슷한지
    - impulse response roundtrip
    - sine wave roundtrip (fs/4, fs/10, fs/20 등 다양한 주파수)
  - **Alias suppression**: Nyquist 근처 주파수 신호가 제대로 감쇠되는지
    - 0.48*fs 신호 input → alias bin에서 에너지 측정
  - **Latency**: `latencySamples()` 계산이 정확한지
    - delay buffer fill 후 latency 측정
  - **Stage combinations**: 1x / 2x / 4x / 8x 모두 테스트
    - `setStages(1)`, `setStages(2)`, `setStages(3)`
  - **Block size variations**: `prepare(64)` then feed different-sized blocks
    - should work up to capacity
    - should assert beyond capacity
- [ ] CMakeLists.txt에 테스트 추가

**참고**: Oversampler는 "소리가 좋아 보이지만 사실 alias를 방지 못했다"는 버그가 가장 흔하므로, 특히 중요

---

### Issue #7: Rate independence verification
**Priority**: HIGH  
**Type**: Testing  
**Effort**: 2-3 hours

**문제**: Rate.h의 rescalePole / rescaleNoise가 정말 rate-independent인지 검증 필요

**할 일**:
- [ ] Tests/RateIndependenceAdvancedTests.cpp 생성 (existing RateIndependenceTests.cpp 확장)
- [ ] 다음 검사:
  - **Time constant invariance**:
    - 같은 pole @ 48kHz, 96kHz, 192kHz에서 step response time constant 동일한가?
    - tolerance: 0.5% 이내
  - **Noise floor invariance**:
    - pink noise 생성 후, 48kHz vs 96kHz에서 integrated noise power 동일한가?
  - **Multi-stage**:
    - InductorResonator + oversampling @ 48kHz vs 8x 384kHz에서 결과 동일한가?
  - **ClassAStage**:
    - dielectric absorption time constants가 rate-independent인가?
  - **Cascaded**: 여러 rate-compensated stage가 조합될 때도 동일한가?
- [ ] CMakeLists.txt에 테스트 추가

**참고**: rate independence는 eon-platform의 "계약"이므로, 이게 깨지면 오디오 품질이 급격히 떨어짐

---

### Issue #8: NaN/Inf safety tests
**Priority**: MEDIUM  
**Type**: Testing  
**Effort**: 1-2 hours

**문제**: invalid input (NaN, Inf, 0, negative)가 들어올 때 crash 하지 않는지 검증

**할 일**:
- [ ] Tests/NaNInfSafetyTests.cpp 생성
- [ ] 다음 케이스 테스트:
  - Wdf 모든 element: `incident(NaN)`, `incident(Inf)`, `incident(-1e308)`
  - Triode: `process(NaN)`, `process(Inf)`
  - Oversampler: 입력에 NaN/Inf 섞여있을 때
  - all stages with invalid parameters
- [ ] assert 아님 graceful fallback 확인
- [ ] output이 유한한 값인지 확인

**참고**: audio plugin은 host로부터 invalid input을 받을 수 있으므로, 필수

---

## 우선순위 3: API 문서화

### Issue #9: Doxygen/code comments 강화
**Priority**: MEDIUM  
**Type**: Documentation  
**Effort**: 3-4 hours

**문제**: 코드 수준의 주석이 이론은 있지만, 사용 방법이 명확하지 않음

**할 일**:
- [ ] Dsp/의 각 struct에 doxygen comment 추가
  ```cpp
  /// @class WdfPort
  /// @brief Base port for wave digital filter elements.
  /// @details
  ///   Manages incident and reflected waves. Subclasses must implement
  ///   emitted() for reflected wave calculation.
  ///   
  /// @note Must be used within a WdfTree; standalone ports are unsafe.
  struct WdfPort { ... };
  ```
- [ ] 각 함수에 @param, @return, @warning 추가
- [ ] Rate.h, Zdf.h, Wdf.h의 이론적 배경 수식 추가 (주석)
- [ ] Doxygen config 파일 생성 (CMakeLists.txt에서 `doxygen_add_docs` 사용)
- [ ] docs/ 폴더에 generated HTML 배포

**참고**: 고급 DSP 코드는 "수학적 설명"이 없으면 유지보수 불가능

---

### Issue #10: API contract documentation
**Priority**: MEDIUM  
**Type**: Documentation  
**Effort**: 2 hours

**문제**: prepare() / reset() / process() 호출 순서와 조건이 명시 부족

**할 일**:
- [ ] Dsp/API_CONTRACTS.md 생성
- [ ] 각 모듈별로 작성:
  ```markdown
  ## DCBlocker
  - **Initialization**: None required, unprepared default R=0.9997
  - **Preparation**: Call prepare(sampleRate, cutoffHz) before processing
  - **Reset**: Call reset() when starting new signal epoch
  - **Threading**: One instance per channel (not thread-safe)
  - **Statefullness**: Carries x1, y1 between samples
  
  ## WdfDiodePair
  - **Required**: R, Is, Vt must be finite and > 0
  - **Iterations**: 8 default, can be 0 for fallback to vPrev
  - **Warm-start**: Uses vPrev to accelerate convergence
  - **solveSucceeded**: Check after emitted() to verify convergence
  - **Invalid input**: aIn=NaN returns 0.0, sets solveSucceeded=false
  ```
- [ ] README에서 링크

**참고**: API misuse는 고급 코드의 가장 흔한 오류 원인

---

## 우선순위 4: 예제 & 통합

### Issue #11: Add minimal working examples
**Priority**: MEDIUM  
**Type**: Examples  
**Effort**: 2-3 hours

**문제**: repo에 예제가 없어서 어떻게 쓰는지 불명

**할 일**:
- [ ] examples/ 폴더 생성
- [ ] examples/01_simple_dcblocker.cpp
  ```cpp
  #include "Dsp/Stages.h"
  
  int main() {
      eon::DCBlocker blocker;
      blocker.prepare(48000.0, 20.0);  // 48kHz, 20Hz cutoff
      blocker.reset();
      
      float x = 0.5f;
      float y = blocker.process(x);  // process single sample
      return 0;
  }
  ```
- [ ] examples/02_zdf_ladder.cpp
  - simple lowpass with cutoff sweep
- [ ] examples/03_triode_stage.cpp
  - with and without cathode self-bias
- [ ] examples/04_oversampler_usage.cpp
  - up/down with rate compensation
- [ ] examples/05_complete_signal_chain.cpp
  - combines multiple stages with oversampling
- [ ] CMakeLists.txt에서 examples 빌드 옵션 추가
- [ ] README에서 각 예제 링크

**참고**: "동작하는 최소 예제"는 사용자가 시작하는 데 가장 중요

---

### Issue #12: Add header compile isolation tests
**Priority**: MEDIUM  
**Type**: Testing  
**Effort**: 1 hour

**문제**: 현재 CMakeLists.txt에 header compile checks가 있지만, 더 엄격하게 검증 가능

**할 일**:
- [ ] Tests/HeaderCompile/ 확장
  - 각 헤더를 **완전 독립적으로** compile해보기
  - include order 다양하게 섞어보기
  - forward declaration만 사용할 때도 테스트
- [ ] `#pragma once` 검증
- [ ] include guard 중복 확인
- [ ] CMakeLists.txt 강화

**참고**: header-only library는 include order 버그가 나중에 터질 수 있음

---

## 우선순위 5: CI/CD & 배포

### Issue #13: GitHub Actions CI workflow 강화
**Priority**: MEDIUM  
**Type**: CI/CD  
**Effort**: 2 hours

**문제**: .github/workflows가 있는지 확인 후 강화 필요

**할 일**:
- [ ] .github/workflows/ci.yml 확인 or 생성
- [ ] matrix build 추가:
  ```yaml
  strategy:
    matrix:
      os: [ubuntu-latest, macos-latest, windows-latest]
      compiler: [gcc, clang, msvc]
      build_type: [Debug, Release]
      cxx_standard: [20]
  ```
- [ ] 다음 단계 실행:
  - `cmake -S . -B build-${{ matrix.build_type }}`
  - `cmake --build build-${{ matrix.build_type }} --parallel`
  - `ctest --test-dir build-${{ matrix.build_type }} --output-on-failure`
- [ ] AddressSanitizer/UndefinedBehaviorSanitizer job 추가
- [ ] Doxygen docs generation & deploy
- [ ] code coverage 보고 (codecov)

**참고**: multi-platform build는 portable audio code에 필수

---

### Issue #14: Changelog & versioning
**Priority**: MEDIUM  
**Type**: Release  
**Effort**: 1 hour

**문제**: 버전 관리 없이 계속 개발되고 있음

**할 일**:
- [ ] CHANGELOG.md 생성
  ```markdown
  ## [Unreleased]
  
  ## [0.1.0] - 2026-10-07
  ### Added
  - Initial eon-platform release
  - Dsp/Rate.h rate compensation
  - Wdf.h wave digital filter
  - Triode.h Koren triode model
  - Oversampling.h half-band FIR
  
  ### Known Issues
  - WdfDiodePair edge case coverage incomplete
  ```
- [ ] version.h or CMakeLists.txt VERSION 정의
- [ ] semantic versioning 정책 문서화
- [ ] GitHub releases 연동

**참고**: "이 코드가 언제 변했는가"는 사용자 입장에서 중요함

---

## 우선순위 6: 고급 검증

### Issue #15: Fuzzing/random property tests
**Priority**: LOW  
**Type**: Testing  
**Effort**: 3-4 hours

**문제**: random input에서 crash나 incorrect behavior 가능성

**할 일**:
- [ ] Tests/FuzzTests.cpp 생성 (libFuzzer or random)
- [ ] 수천 개의 random input sweep:
  - random `aIn`, `R`, `Is`, `Vt` for WdfDiodePair
  - random `freqHz`, `q`, `drive` for filters
  - random block sizes for oversampler
- [ ] invariants 검증:
  - output is finite
  - state doesn't grow unbounded
  - reset() always stabilizes
- [ ] longevity test (장시간 런)

**참고**: high-end audio는 "한 번 깨면 끝"이므로, robustness 중요

---

### Issue #16: Benchmark & performance profiling
**Priority**: LOW  
**Type**: Performance  
**Effort**: 2-3 hours

**문제**: DSP code의 성능 baseline 없음

**할 일**:
- [ ] benchmarks/ 폴더 생성
- [ ] Google Benchmark 또는 simple timer로:
  - WdfDiodePair::emitted() per-call time
  - TriodeStage::process() throughput
  - Oversampler::up/down throughput
  - 전체 signal chain benchmark
- [ ] CMakeLists.txt에 benchmark target 추가
- [ ] CI에서 benchmark 실행 & 결과 저장
- [ ] regression detection

**참고**: "얼마나 빠른가"는 플러그인 호스트 부하 결정

---

## 우선순위 7: 커뮤니티 & 지원

### Issue #17: Contributing guide 작성
**Priority**: LOW  
**Type**: Community  
**Effort**: 1 hour

**할 일**:
- [ ] CONTRIBUTING.md 생성
  ```markdown
  # Contributing to eon-platform
  
  ## Development Setup
  - Requires C++20 compiler
  - CMake 3.16+
  - Run: cmake -S . -B build && cmake --build build
  
  ## Code Style
  - 4-space indentation
  - noexcept for pure DSP functions
  - comprehensive comments for complex algorithms
  
  ## Before Submitting PR
  - All tests pass: ctest --test-dir build
  - No new warnings (with -Wall -Wextra -Wpedantic)
  - Update CHANGELOG.md
  - Add tests for new code
  ```
- [ ] DEVELOPMENT.md (deeper guide)
- [ ] Code of Conduct

---

### Issue #18: Discussion/feedback channels
**Priority**: LOW  
**Type**: Community  
**Effort**: minimal

**할 일**:
- [ ] GitHub Discussions 활성화
- [ ] README에 discussion link
- [ ] "버그 보고" vs "기능 제안" issue template 생성

---

## 요약 테이블

| # | 제목 | 우선순위 | 종류 | 노력 | 상태 |
|---|------|---------|------|------|------|
| 1 | License 추가 | CRITICAL | Setup | 5m | TODO |
| 2 | README.md 작성 | CRITICAL | Docs | 2-3h | TODO |
| 3 | .gitignore 정리 | MEDIUM | Setup | 30m | TODO |
| 4 | WdfDiodePair tests | HIGH | Test | 2-3h | TODO |
| 5 | TriodeStage tests | HIGH | Test | 2h | TODO |
| 6 | Oversampler tests | HIGH | Test | 2-3h | TODO |
| 7 | Rate independence tests | HIGH | Test | 2-3h | TODO |
| 8 | NaN/Inf safety tests | MEDIUM | Test | 1-2h | TODO |
| 9 | Doxygen comments | MEDIUM | Docs | 3-4h | TODO |
| 10 | API contracts | MEDIUM | Docs | 2h | TODO |
| 11 | Minimal examples | MEDIUM | Examples | 2-3h | TODO |
| 12 | Header isolation tests | MEDIUM | Test | 1h | TODO |
| 13 | GitHub Actions CI | MEDIUM | CI/CD | 2h | TODO |
| 14 | Changelog & versioning | MEDIUM | Release | 1h | TODO |
| 15 | Fuzzing tests | LOW | Test | 3-4h | TODO |
| 16 | Benchmarking | LOW | Perf | 2-3h | TODO |
| 17 | Contributing guide | LOW | Docs | 1h | TODO |
| 18 | Discussion channels | LOW | Community | minimal | TODO |

---

## 추천 로드맵

**Phase 1 (1-2주)**: 기본 기초
- Issue #1: License
- Issue #2: README
- Issue #3: .gitignore

**Phase 2 (2-3주)**: 테스트 강화
- Issue #4, #5, #6, #7: 핵심 모듈 테스트
- Issue #8: Safety tests

**Phase 3 (1-2주)**: 문서화
- Issue #9, #10: Code docs
- Issue #11, #12: Examples & headers

**Phase 4 (1주)**: CI/Release
- Issue #13, #14: CI & versioning

**Phase 5 (optional)**: 고급
- Issue #15, #16, #17, #18: 커뮤니티 & 성능
