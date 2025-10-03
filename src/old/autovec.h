#pragma once
#include <vector>
#include <functional>
#include <algorithm>
#include <cstddef>

namespace Makebr {

// 접근 시 자동 확장되는 벡터 + "증가 단위(step)" 설정
template <typename T>
class AutoVec {
public:
  using Factory = std::function<T(std::size_t)>;

  // 기본 채움값 (step 지정 가능)
  explicit AutoVec(const T& def = T{}, std::size_t step = 1)
    : def_(def), growth_step_(step ? step : 1) {}

  // 팩토리(새 원소 생성기) (step 지정 가능)
  explicit AutoVec(Factory make, std::size_t step = 1)
    : make_(std::move(make)), growth_step_(step ? step : 1) {}

  // 초기 길이 + 기본값 (step 지정 가능)
  AutoVec(std::size_t n, const T& def, std::size_t step = 1)
    : v_(n, def), def_(def), growth_step_(step ? step : 1) {}

  // 읽기/쓰기: 없으면 자동 확장(증가 단위 적용)
  T& operator[](std::size_t i) {
    if (i >= v_.size()) grow_to(i + 1);
    return v_[i];
  }
  // 읽기 전용: 확장하지 않음
  const T& operator[](std::size_t i) const { return v_[i]; }

  std::size_t size() const { return v_.size(); }
  void clear() { v_.clear(); }
  void reserve(std::size_t n) { v_.reserve(n); }

  // 명시적 리사이즈: 증가 단위 적용
  void resize(std::size_t n) { grow_to(n); }

  // 증가 단위 설정/조회
  void set_growth_step(std::size_t step) { growth_step_ = step ? step : 1; }
  std::size_t growth_step() const { return growth_step_; }

  // 내부 벡터 접근(필요시)
  const std::vector<T>& vec() const { return v_; }
        std::vector<T>& vec()       { return v_; }

private:
  static std::size_t round_up(std::size_t x, std::size_t step) {
    return (x + step - 1) / step * step;
  }

  void grow_to(std::size_t new_len) {
    // 요청 길이를 step의 배수로 올림
    const std::size_t target_len = round_up(new_len, growth_step_);
    if (target_len <= v_.size()) return;

    // capacity는 1.5배 정도로 기하성장 → 재할당 횟수 감소
    std::size_t cap = v_.capacity();
    if (cap < target_len) {
      std::size_t add = std::max<std::size_t>(8, cap/2);
      std::size_t new_cap = std::max(target_len, cap + add);
      v_.reserve(new_cap);
    }

    std::size_t old = v_.size();
    v_.resize(target_len);

    if (make_) {
      for (std::size_t k = old; k < target_len; ++k) v_[k] = make_(k);
    } else {
      std::fill(v_.begin() + old, v_.end(), def_);
    }
  }

  std::vector<T> v_;
  T        def_{};                  // factory 미사용 시 채움값
  Factory  make_{nullptr};          // 새 원소 생성기(필요 시)
  std::size_t growth_step_ = 1;     // 증가 단위(최소 1)
};

} // namespace Makebr
