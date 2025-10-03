#pragma once
#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>

namespace slog {

// Extract Filename
inline const char* basename_cstr(const char* path) {
  const char* base = path;
  for (const char* p = path; *p; ++p) {
    if (*p == '/' || *p == '\\') base = p + 1;
  }
  return base;
}

// Time Stamp
inline std::string timestamp_yy_mm_dd_ss() {
  using namespace std::chrono;
  const auto now = system_clock::now();
  const auto tt  = system_clock::to_time_t(now);
  std::tm tm{};
#if defined(_WIN32)
  localtime_s(&tm, &tt);
#else
  localtime_r(&tt, &tm);
#endif
  std::ostringstream oss;
  oss << '[' << std::put_time(&tm, "%y:%m:%d:%S") << ']'; // YY:MM:DD:SS
  return oss.str();
}

struct Config {
#ifdef LOGGER_THREAD_SAFE
  std::mutex mu;
#endif
};

inline Config& cfg() { static Config c; return c; }

// Log by line
class Line {
public:
  Line(const char* code, int line) : code_(code), line_(line) {}
  ~Line() {
#ifdef LOGGER_THREAD_SAFE
    std::lock_guard<std::mutex> lock(cfg().mu);
#endif
    std::ostream& out = std::cout;
    out << timestamp_yy_mm_dd_ss()
        << " [" << code_ << "] "
        << "[" << line_ << "] "
        << "' " << ss_.str() << "'\n";
  }
  template<typename T>
  Line& operator<<(const T& v) { ss_ << v; return *this; }

private:
  const char* code_;
  int line_;
  std::ostringstream ss_;
};

} // namespace slog


// [YY:MM:DD:SS] [main.cpp] [123] ' message'
#define LOG()   slog::Line(slog::basename_cstr(__FILE__), __LINE__)
// [YY:MM:DD:SS] [MYCODE] [123] ' message'
#define LOGC(x) slog::Line((x), __LINE__)
