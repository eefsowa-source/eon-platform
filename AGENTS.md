# eon_dsp

- `Dsp/`의 공유 헤더가 라이브러리 계약이며 `Tests/`가 검증 위치다. `CMakeLists.txt`는 인터페이스 라이브러리와 CTest 대상을 정의한다.
- DSP 동작을 바꾸면 해당 테스트와 헤더 검사 대상을 빌드하고 `ctest --test-dir <빌드 디렉터리> --output-on-failure`를 실행한다.
- 플러그인 저장소에 복사된 코어와의 관계는 `NOTES.md` 및 `/Users/sungha/Documents/Codex/audio-knowledge/validation/shared-dsp-core-vendoring.md`를 확인한다. 복사본을 조용히 동기화하지 않는다.
