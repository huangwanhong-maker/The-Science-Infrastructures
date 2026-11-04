// DICOM images read and placed in the patient's frame: the slices of a series sorted along their normal, their
// spacing taken from their positions, the grid's geometry in LPS, and NIfTI's RAS carried to it (ADR 0070).
//
// The patient frame of DICOM (PS3.3 C.7.6.2.1.1) is LPS: x toward the patient's left, y posterior, z superior, in
// millimetres. A slice's Image Position (Patient) (0020,0032) is the centre of its first pixel; its Image
// Orientation (Patient) (0020,0037) gives the row direction (along which the column index increases) and then
// the column direction (along which the row index increases). Pixel Spacing (0028,0030) is given row spacing
// first: the distance between rows, then between columns (PS3.3 10.7.1.3). So with i the column index and j the
// row index, as the grid stores them, spacing_i is Pixel Spacing's second value and spacing_j its first: reading
// them in order transposes an anisotropic pixel, the trap the checks guard.
//
// A series' slices are sorted by the projection of their positions on the normal n = row × column, and the
// spacing between slices is taken from successive positions, not from Slice Thickness (0018,0050), which states
// the reconstructed thickness and may differ from the spacing (overlapping or gapped reconstructions). Uneven
// spacing (a missing slice) and positions that move within the plane (gantry tilt, a sheared grid) are reported;
// a Geometry is made only of an even, unsheared series.
//
// Files (PS3.10): the 128-byte preamble and "DICM", the file meta group in explicit VR little endian, then the
// data set in implicit or explicit VR little endian (transfer syntaxes 1.2.840.10008.1.2 and .1.2.1). Compressed
// transfer syntaxes are refused with their UID: the data library decompresses them. Sequences, of defined or
// undefined length, are skipped. Values are scaled by Rescale Slope and Intercept (CT's Hounsfield units).
//
// NIfTI's RAS (eikon's read_nifti states it) differs from LPS by negating x and y; `to_lps` carries a geometry
// across through the anatomy domain's frames. Checks: tests/check_dicom.cpp.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <anatomy/frames/frames.hpp>
#include <eikon/grid/grid.hpp>

namespace medical_imaging::dicom {

using eikon::grid::Geometry;
using eikon::grid::Image;
using eikon::grid::Vec3;

// What a slice says of itself.
struct Slice {
  std::string modality, series;
  int instance = 0;
  Vec3 position{0, 0, 0};                         // Image Position (Patient), mm, LPS
  Vec3 row{1, 0, 0}, column{0, 1, 0};             // Image Orientation (Patient)
  double row_spacing = 1, column_spacing = 1;     // Pixel Spacing: between rows, between columns
  double thickness = 0;
  std::size_t rows = 0, columns = 0;
  double slope = 1, intercept = 0;
  std::vector<float> values;                      // rows × columns, the column index fastest, rescaled
};

namespace detail {
inline std::uint16_t u16(const std::vector<char>& b, std::size_t at) { std::uint16_t v; std::memcpy(&v, b.data() + at, 2); return v; }
inline std::uint32_t u32(const std::vector<char>& b, std::size_t at) { std::uint32_t v; std::memcpy(&v, b.data() + at, 4); return v; }
inline bool long_vr(const char* vr) {
  static const char* L[] = {"OB", "OD", "OF", "OL", "OV", "OW", "SQ", "SV", "UC", "UN", "UR", "UT", "UV"};
  for (auto l : L)
    if (vr[0] == l[0] && vr[1] == l[1]) return true;
  return false;
}
inline std::vector<double> decimals(const std::string& s) {
  std::vector<double> out;
  std::size_t at = 0;
  while (at <= s.size()) {
    const std::size_t end = std::min(s.find('\\', at), s.size());
    const std::string part = s.substr(at, end - at);
    if (part.find_first_not_of(" \t\0", 0, 3) != std::string::npos) out.push_back(std::stod(part));
    at = end + 1;
  }
  return out;
}
inline std::string trimmed(std::string s) {
  while (!s.empty() && (s.back() == ' ' || s.back() == '\0')) s.pop_back();
  std::size_t a = 0;
  while (a < s.size() && s[a] == ' ') a++;
  return s.substr(a);
}
}  // namespace detail

// One slice read from a DICOM file (uncompressed, little endian).
inline Slice read_slice(const std::filesystem::path& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("dicom: cannot read " + path.string());
  std::vector<char> b((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  if (b.size() < 132 || std::memcmp(b.data() + 128, "DICM", 4) != 0) throw std::runtime_error("dicom: not a DICOM Part 10 file: " + path.string());
  std::size_t at = 132;
  bool explicit_vr = true;
  std::map<std::uint32_t, std::pair<std::size_t, std::size_t>> found;  // tag → (offset, length)
  std::string syntax;
  bool meta = true;
  const auto skip_items = [&](auto&& self, std::size_t& p, bool expl) -> void {
    // within an undefined-length sequence: items until the sequence delimiter
    while (p + 8 <= b.size()) {
      const std::uint32_t tag = (std::uint32_t(detail::u16(b, p)) << 16) | detail::u16(b, p + 2);
      const std::uint32_t len = detail::u32(b, p + 4);
      p += 8;
      if (tag == 0xFFFEE0DD) return;  // sequence delimiter
      if (tag != 0xFFFEE000) throw std::runtime_error("dicom: a malformed sequence in " + path.string());
      if (len != 0xFFFFFFFF) { p += len; continue; }
      // an item of undefined length: elements until the item delimiter
      while (p + 8 <= b.size()) {
        const std::uint32_t t = (std::uint32_t(detail::u16(b, p)) << 16) | detail::u16(b, p + 2);
        if (t == 0xFFFEE00D) { p += 8; break; }
        std::uint32_t l;
        bool sq = false;
        if (expl && t >> 16 != 0xFFFE) {
          const char* vr = b.data() + p + 4;
          sq = vr[0] == 'S' && vr[1] == 'Q';
          if (detail::long_vr(vr)) l = detail::u32(b, p + 8), p += 12;
          else l = detail::u16(b, p + 6), p += 8;
        } else {
          l = detail::u32(b, p + 4), p += 8;
        }
        if (l == 0xFFFFFFFF) { (void)sq; self(self, p, expl); }
        else p += l;
      }
    }
  };
  while (at + 8 <= b.size()) {
    if (meta && detail::u16(b, at) != 0x0002) {
      meta = false;
      if (syntax == "1.2.840.10008.1.2") explicit_vr = false;
      else if (syntax == "1.2.840.10008.1.2.1" || syntax.empty()) explicit_vr = true;
      else throw std::runtime_error("dicom: transfer syntax " + syntax + " (compressed or big endian) is not read here: " + path.string());
    }
    const std::uint32_t tag = (std::uint32_t(detail::u16(b, at)) << 16) | detail::u16(b, at + 2);
    std::uint32_t len;
    bool sq = false;
    if (meta || explicit_vr) {
      const char* vr = b.data() + at + 4;
      sq = vr[0] == 'S' && vr[1] == 'Q';
      if (detail::long_vr(vr)) len = detail::u32(b, at + 8), at += 12;
      else len = detail::u16(b, at + 6), at += 8;
    } else {
      len = detail::u32(b, at + 4), at += 8;
    }
    if (len == 0xFFFFFFFF) {
      if (tag == 0x7FE00010) throw std::runtime_error("dicom: encapsulated (compressed) pixel data in " + path.string());
      skip_items(skip_items, at, meta || explicit_vr);
      continue;
    }
    (void)sq;
    if (at + len > b.size()) throw std::runtime_error("dicom: an element runs past the end of " + path.string());
    found[tag] = {at, len};
    if (tag == 0x00020010) syntax = detail::trimmed(std::string(b.data() + at, len));
    at += len;
  }
  const auto text = [&](std::uint32_t tag) -> std::string {
    auto it = found.find(tag);
    return it == found.end() ? "" : detail::trimmed(std::string(b.data() + it->second.first, it->second.second));
  };
  const auto us = [&](std::uint32_t tag) -> std::optional<std::uint16_t> {
    auto it = found.find(tag);
    if (it == found.end() || it->second.second < 2) return std::nullopt;
    return detail::u16(b, it->second.first);
  };
  Slice s;
  s.modality = text(0x00080060);
  s.series = text(0x0020000E);
  if (auto t = text(0x00200013); !t.empty()) s.instance = std::stoi(t);
  const auto ipp = detail::decimals(text(0x00200032)), iop = detail::decimals(text(0x00200037)), ps = detail::decimals(text(0x00280030));
  if (ipp.size() != 3 || iop.size() != 6 || ps.size() != 2) throw std::runtime_error("dicom: position, orientation or pixel spacing missing in " + path.string());
  s.position = {ipp[0], ipp[1], ipp[2]};
  s.row = {iop[0], iop[1], iop[2]};
  s.column = {iop[3], iop[4], iop[5]};
  s.row_spacing = ps[0], s.column_spacing = ps[1];
  if (auto t = detail::decimals(text(0x00180050)); !t.empty()) s.thickness = t[0];
  if (auto t = detail::decimals(text(0x00281053)); !t.empty()) s.slope = t[0];
  if (auto t = detail::decimals(text(0x00281052)); !t.empty()) s.intercept = t[0];
  s.rows = us(0x00280010).value_or(0), s.columns = us(0x00280011).value_or(0);
  const int bits = us(0x00280100).value_or(16), signed_ = us(0x00280103).value_or(0), samples = us(0x00280002).value_or(1);
  if (samples != 1) throw std::runtime_error("dicom: only one sample per pixel is read: " + path.string());
  auto px = found.find(0x7FE00010);
  if (px == found.end()) throw std::runtime_error("dicom: no pixel data in " + path.string());
  const std::size_t n = s.rows * s.columns;
  if (px->second.second < n * std::size_t(bits / 8)) throw std::runtime_error("dicom: short pixel data in " + path.string());
  s.values.resize(n);
  for (std::size_t q = 0; q < n; q++) {
    double v;
    const std::size_t o = px->second.first + q * std::size_t(bits / 8);
    if (bits == 8) v = signed_ ? double(std::int8_t(b[o])) : double(std::uint8_t(b[o]));
    else if (bits == 16) v = signed_ ? double(std::int16_t(detail::u16(b, o))) : double(detail::u16(b, o));
    else if (bits == 32) v = signed_ ? double(std::int32_t(detail::u32(b, o))) : double(detail::u32(b, o));
    else throw std::runtime_error("dicom: " + std::to_string(bits) + " bits allocated is not read here");
    s.values[q] = float(s.slope * v + s.intercept);
  }
  return s;
}

// ---- series ------------------------------------------------------------------------------------------------

struct SeriesReport {
  std::size_t slices = 0;
  double spacing = 0;            // the median distance between successive slices along the normal
  double spacing_spread = 0;     // the largest departure of a step from the median
  double shear = 0;              // the largest in-plane movement of a step (gantry tilt)
  double thickness = 0;          // Slice Thickness as stated (may differ from the spacing)
  bool even = true, unsheared = true;
  std::vector<std::size_t> gaps; // steps (k → k + 1) longer than 1.5 spacings
};

// The slices sorted along their normal and the report of their spacing.
inline SeriesReport sort_series(std::vector<Slice>& slices, double tolerance = 0.01) {
  SeriesReport r;
  r.slices = slices.size();
  if (slices.empty()) return r;
  const Vec3 n = {slices[0].row[1] * slices[0].column[2] - slices[0].row[2] * slices[0].column[1],
                  slices[0].row[2] * slices[0].column[0] - slices[0].row[0] * slices[0].column[2],
                  slices[0].row[0] * slices[0].column[1] - slices[0].row[1] * slices[0].column[0]};
  const auto along = [&](const Slice& s) { return s.position[0] * n[0] + s.position[1] * n[1] + s.position[2] * n[2]; };
  std::sort(slices.begin(), slices.end(), [&](const Slice& a, const Slice& b) { return along(a) < along(b); });
  r.thickness = slices[0].thickness;
  if (slices.size() < 2) { r.spacing = slices[0].thickness > 0 ? slices[0].thickness : 1; return r; }
  std::vector<double> steps;
  for (std::size_t k = 1; k < slices.size(); k++) {
    const Vec3 d{slices[k].position[0] - slices[k - 1].position[0], slices[k].position[1] - slices[k - 1].position[1], slices[k].position[2] - slices[k - 1].position[2]};
    const double a = d[0] * n[0] + d[1] * n[1] + d[2] * n[2];
    steps.push_back(a);
    const Vec3 inplane{d[0] - a * n[0], d[1] - a * n[1], d[2] - a * n[2]};
    r.shear = std::max(r.shear, std::sqrt(inplane[0] * inplane[0] + inplane[1] * inplane[1] + inplane[2] * inplane[2]));
  }
  auto sorted = steps;
  std::nth_element(sorted.begin(), sorted.begin() + std::ptrdiff_t(sorted.size() / 2), sorted.end());
  r.spacing = sorted[sorted.size() / 2];
  for (std::size_t k = 0; k < steps.size(); k++) {
    r.spacing_spread = std::max(r.spacing_spread, std::fabs(steps[k] - r.spacing));
    if (steps[k] > 1.5 * r.spacing) r.gaps.push_back(k);
  }
  r.even = r.spacing_spread <= tolerance * r.spacing + 1e-6;
  r.unsheared = r.shear <= tolerance * r.spacing + 1e-6;
  return r;
}

// The grid of a sorted, even and unsheared series, in LPS; throws otherwise.
inline Geometry series_geometry(const std::vector<Slice>& slices, const SeriesReport& r) {
  if (slices.empty()) throw std::invalid_argument("dicom: an empty series");
  if (!r.even) throw std::invalid_argument("dicom: uneven slice spacing (a missing slice?); resample before making a grid");
  if (!r.unsheared) throw std::invalid_argument("dicom: slice positions move within the plane (gantry tilt); resample before making a grid");
  const auto& s = slices[0];
  Geometry g;
  g.size = {s.columns, s.rows, slices.size()};
  g.spacing = {s.column_spacing, s.row_spacing, r.spacing};  // PixelSpacing is row spacing first
  const Vec3 n = {s.row[1] * s.column[2] - s.row[2] * s.column[1], s.row[2] * s.column[0] - s.row[0] * s.column[2], s.row[0] * s.column[1] - s.row[1] * s.column[0]};
  for (int a = 0; a < 3; a++) g.direction[a][0] = s.row[a], g.direction[a][1] = s.column[a], g.direction[a][2] = n[a];
  g.origin = s.position;
  return g;
}

// The image of a series: its slices read, sorted, checked and stacked.
inline Image<float> read_series(const std::vector<std::filesystem::path>& files, SeriesReport* report = nullptr) {
  std::vector<Slice> slices;
  for (const auto& f : files) slices.push_back(read_slice(f));
  const auto r = sort_series(slices);
  if (report) *report = r;
  Image<float> im(series_geometry(slices, r));
  const std::size_t plane = slices[0].rows * slices[0].columns;
  for (std::size_t k = 0; k < slices.size(); k++) {
    if (slices[k].values.size() != plane) throw std::invalid_argument("dicom: slices of different sizes in one series");
    std::copy(slices[k].values.begin(), slices[k].values.end(), im.data.begin() + std::ptrdiff_t(k * plane));
  }
  return im;
}

// A geometry in NIfTI's RAS carried to LPS: x and y of the origin and of every direction negated.
inline Geometry to_lps(const Geometry& ras) {
  Geometry g = ras;
  g.origin = anatomy::frames::convert(ras.origin, anatomy::frames::Frame::RAS, anatomy::frames::Frame::LPS);
  for (int c = 0; c < 3; c++) {
    const auto d = anatomy::frames::convert({ras.direction[0][c], ras.direction[1][c], ras.direction[2][c]}, anatomy::frames::Frame::RAS, anatomy::frames::Frame::LPS);
    for (int a = 0; a < 3; a++) g.direction[a][c] = d[a];
  }
  return g;
}

}  // namespace medical_imaging::dicom
