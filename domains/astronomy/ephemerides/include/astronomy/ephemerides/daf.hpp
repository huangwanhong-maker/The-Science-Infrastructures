// NAIF's Double precision Array Files: the container of the JPL planetary ephemerides (SPK) and
// of binary orientation data (PCK).
//
// A DAF is a sequence of 1024-byte records: a file record (identification, the numbers ND of
// doubles and NI of integers in each summary, the first and last summary records, the binary
// format), then doubly linked summary records, each with up to 125 / (ND + ⌈NI/2⌉) summaries and
// a name record after it, and the arrays themselves, addressed in doubles from 1 (NAIF, "DAF
// Required Reading", N0067). A summary holds ND doubles (for SPK and PCK, the first and last
// instant covered, TDB seconds past J2000) and NI integers (for SPK: target, centre, frame, data
// type, first and last address; for PCK: body frame, base frame, data type, first and last).
//
// The reader takes the file's bytes, or a path the caller gives; it accepts the little-endian
// IEEE format ("LTL-IEEE") on a little-endian machine, which is every machine the infrastructure
// runs on, and refuses anything else.

#pragma once

#include <bit>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace astronomy::ephemerides {

struct DafSummary {
  std::vector<double> doubles;  // ND of them
  std::vector<std::int32_t> ints;  // NI of them
  std::string name;
};

class Daf {
 public:
  static constexpr std::size_t RECORD = 1024;

  explicit Daf(std::vector<char> bytes) : bytes_(std::move(bytes)) {
    if constexpr (std::endian::native != std::endian::little) throw std::runtime_error("DAF: a little-endian machine is required");
    if (bytes_.size() < RECORD) throw std::invalid_argument("DAF: shorter than one record");
    id_ = std::string(bytes_.data(), 8);
    if (id_.rfind("DAF/", 0) != 0 && id_.rfind("NAIF/DAF", 0) != 0) throw std::invalid_argument("DAF: not a DAF file (\"" + id_ + "\")");
    const std::string format(bytes_.data() + 88, 8);
    if (format != "LTL-IEEE") throw std::invalid_argument("DAF: binary format " + format + " is not supported");
    nd_ = int32_at(8);
    ni_ = int32_at(12);
    if (nd_ < 0 || ni_ < 2 || nd_ > 124) throw std::invalid_argument("DAF: implausible ND and NI");
    const std::size_t words = std::size_t(nd_) + std::size_t(ni_ + 1) / 2;  // a summary, in doubles
    std::int32_t record = int32_at(76);  // FWARD
    int guard = 0;
    while (record > 0) {
      if (++guard > 100000) throw std::invalid_argument("DAF: the summary records loop");
      const std::size_t base = std::size_t(record - 1) * RECORD;
      if (base + RECORD > bytes_.size()) throw std::invalid_argument("DAF: a summary record lies beyond the file");
      const double next = double_at(base), count = double_at(base + 16);
      const std::size_t name_base = base + RECORD;
      const std::size_t name_length = 8 * words;
      for (int i = 0; i < int(count); i++) {
        const std::size_t off = base + 24 + std::size_t(i) * words * 8;
        DafSummary s;
        for (int k = 0; k < nd_; k++) s.doubles.push_back(double_at(off + 8 * std::size_t(k)));
        for (int k = 0; k < ni_; k++) s.ints.push_back(int32_at(off + 8 * std::size_t(nd_) + 4 * std::size_t(k)));
        if (name_base + (i + 1) * name_length <= bytes_.size()) {
          std::string n(bytes_.data() + name_base + i * name_length, name_length);
          n.erase(n.find_last_not_of(' ') + 1);
          s.name = n;
        }
        summaries_.push_back(std::move(s));
      }
      record = std::int32_t(next);
    }
  }
  static Daf read(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::invalid_argument("DAF: cannot open " + path);
    std::vector<char> bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    return Daf(std::move(bytes));
  }

  const std::string& id() const { return id_; }
  int nd() const { return nd_; }
  int ni() const { return ni_; }
  const std::vector<DafSummary>& summaries() const { return summaries_; }
  // The double at a 1-based address.
  double at(std::int64_t address) const { return double_at(std::size_t(address - 1) * 8); }
  // A pointer to the doubles from a 1-based address (the file's bytes are aligned to 8 by records).
  const double* doubles(std::int64_t address) const {
    const std::size_t off = std::size_t(address - 1) * 8;
    if (off + 8 > bytes_.size()) throw std::out_of_range("DAF: address beyond the file");
    return reinterpret_cast<const double*>(bytes_.data() + off);
  }

 private:
  std::int32_t int32_at(std::size_t off) const {
    std::int32_t v;
    std::memcpy(&v, bytes_.data() + off, 4);
    return v;
  }
  double double_at(std::size_t off) const {
    if (off + 8 > bytes_.size()) throw std::out_of_range("DAF: read beyond the file");
    double v;
    std::memcpy(&v, bytes_.data() + off, 8);
    return v;
  }
  std::vector<char> bytes_;
  std::string id_;
  int nd_ = 0, ni_ = 0;
  std::vector<DafSummary> summaries_;
};

// A Chebyshev series Σ c_k T_k(s) and its derivative with respect to s, at s ∈ [−1, 1], by the
// recurrences T_{k+1} = 2s T_k − T_{k−1} and T′_{k+1} = 2T_k + 2s T′_k − T′_{k−1}.
inline std::pair<double, double> chebyshev(const double* c, int n, double s) {
  double t0 = 1, t1 = s, d0 = 0, d1 = 1;
  double value = c[0] + (n > 1 ? c[1] * s : 0.0), derivative = n > 1 ? c[1] : 0.0;
  for (int k = 2; k < n; k++) {
    const double t2 = 2 * s * t1 - t0, d2 = 2 * t1 + 2 * s * d1 - d0;
    value += c[k] * t2;
    derivative += c[k] * d2;
    t0 = t1, t1 = t2, d0 = d1, d1 = d2;
  }
  return {value, derivative};
}

}  // namespace astronomy::ephemerides
