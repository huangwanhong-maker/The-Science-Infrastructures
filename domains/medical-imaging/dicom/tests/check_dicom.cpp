// Checks of DICOM series and frames (dicom.hpp; ADR 0070).
//
// Files are written here, in explicit and implicit VR little endian, for a series whose geometry is known: an
// oblique orientation, a non-square pixel of 0.6 mm between columns and 0.8 mm between rows (so that reading
// Pixel Spacing in the wrong order is caught), 2.5 mm between slices, the files shuffled, and a sequence of
// undefined length to skip.
// - The grid read is the one intended: every voxel's world position, to 10⁻⁶ mm (the precision of the decimal
//   strings written, ten significant digits; DICOM's DS holds at most sixteen characters); values rescaled to Hounsfield
//   units; the same from implicit VR.
// - A missing slice is reported as a gap and an uneven series refused; a tilted gantry as shear.
// - A compressed transfer syntax is refused with its UID.
// - NIfTI's RAS carried to LPS: every voxel's position negated in x and y.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <calx/verification/check.hpp>
#include <medical_imaging/dicom/dicom.hpp>

using calx::verification::Suite;
namespace md = medical_imaging::dicom;
namespace eg = eikon::grid;

namespace {

struct Writer {
  std::vector<char> b;
  bool explicit_vr = true;
  void tag(std::uint16_t g, std::uint16_t e) { b.push_back(char(g & 0xFF)), b.push_back(char(g >> 8)), b.push_back(char(e & 0xFF)), b.push_back(char(e >> 8)); }
  void u16(std::uint16_t v) { b.push_back(char(v & 0xFF)), b.push_back(char(v >> 8)); }
  void u32(std::uint32_t v) { for (int k = 0; k < 4; k++) b.push_back(char(v >> (8 * k) & 0xFF)); }
  void element(std::uint16_t g, std::uint16_t e, const char* vr, std::string value, bool force_explicit = false) {
    if (value.size() % 2) value += (std::string(vr) == "UI" ? '\0' : ' ');
    tag(g, e);
    const bool ex = explicit_vr || force_explicit;
    const bool lng = std::string(vr) == "OW" || std::string(vr) == "OB" || std::string(vr) == "SQ" || std::string(vr) == "UN";
    if (ex) {
      b.push_back(vr[0]), b.push_back(vr[1]);
      if (lng) u16(0), u32(std::uint32_t(value.size()));
      else u16(std::uint16_t(value.size()));
    } else {
      u32(std::uint32_t(value.size()));
    }
    b.insert(b.end(), value.begin(), value.end());
  }
  void us(std::uint16_t g, std::uint16_t e, std::uint16_t v) { element(g, e, "US", std::string{char(v & 0xFF), char(v >> 8)}); }
  // a sequence of undefined length holding one item of undefined length with one element
  void sequence(std::uint16_t g, std::uint16_t e) {
    tag(g, e);
    if (explicit_vr) b.push_back('S'), b.push_back('Q'), u16(0);
    u32(0xFFFFFFFF);
    tag(0xFFFE, 0xE000), u32(0xFFFFFFFF);
    element(0x0008, 0x1150, "UI", "1.2.3.4");
    tag(0xFFFE, 0xE00D), u32(0);
    tag(0xFFFE, 0xE0DD), u32(0);
  }
};

std::string ds(double v) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%.10g", v);
  return buf;
}

void write_slice(const std::filesystem::path& p, bool explicit_vr, const std::string& syntax, const eg::Vec3& ipp, const eg::Vec3& row, const eg::Vec3& col,
                 int instance, const std::vector<std::int16_t>& stored) {
  Writer w;
  w.b.assign(128, 0);
  for (char c : std::string("DICM")) w.b.push_back(c);
  w.explicit_vr = true;
  w.element(0x0002, 0x0010, "UI", syntax);
  w.explicit_vr = explicit_vr;
  w.element(0x0008, 0x0060, "CS", "CT");
  w.sequence(0x0008, 0x1140);
  w.element(0x0018, 0x0050, "DS", "3");
  w.element(0x0020, 0x000E, "UI", "1.2.826.0.1.3680043.2.1125.1");
  w.element(0x0020, 0x0013, "IS", std::to_string(instance));
  w.element(0x0020, 0x0032, "DS", ds(ipp[0]) + "\\" + ds(ipp[1]) + "\\" + ds(ipp[2]));
  w.element(0x0020, 0x0037, "DS", ds(row[0]) + "\\" + ds(row[1]) + "\\" + ds(row[2]) + "\\" + ds(col[0]) + "\\" + ds(col[1]) + "\\" + ds(col[2]));
  w.us(0x0028, 0x0002, 1);
  w.us(0x0028, 0x0010, 6);  // rows
  w.us(0x0028, 0x0011, 5);  // columns
  w.element(0x0028, 0x0030, "DS", "0.8\\0.6");  // between rows, between columns
  w.us(0x0028, 0x0100, 16);
  w.us(0x0028, 0x0103, 1);
  w.element(0x0028, 0x1052, "DS", "-1024");
  w.element(0x0028, 0x1053, "DS", "1");
  std::string px(stored.size() * 2, '\0');
  std::memcpy(px.data(), stored.data(), px.size());
  w.element(0x7FE0, 0x0010, "OW", px);
  std::ofstream(p, std::ios::binary).write(w.b.data(), std::streamsize(w.b.size()));
}

}  // namespace

int main() {
  Suite s;
  const auto dir = std::filesystem::temp_directory_path() / "medical-imaging-check-dicom";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  // an oblique frame: row and column directions turned about (1, 2, 2)/3 by 0.4 rad
  const double c = std::cos(0.4), sn = std::sin(0.4), t = 1 - c, x = 1.0 / 3, y = 2.0 / 3, z = 2.0 / 3;
  const double R[3][3] = {{t * x * x + c, t * x * y - sn * z, t * x * z + sn * y}, {t * x * y + sn * z, t * y * y + c, t * y * z - sn * x}, {t * x * z - sn * y, t * y * z + sn * x, t * z * z + c}};
  const eg::Vec3 row{R[0][0], R[1][0], R[2][0]}, col{R[0][1], R[1][1], R[2][1]}, n{R[0][2], R[1][2], R[2][2]};
  const eg::Vec3 o{-120.5, 33.25, 410.0};
  const double dz = 2.5;
  const auto stored = [](int k) {
    std::vector<std::int16_t> v(30);
    for (int q = 0; q < 30; q++) v[std::size_t(q)] = std::int16_t(1000 + 100 * k + q);
    return v;
  };
  const auto series = [&](const std::string& name, bool explicit_vr, std::vector<int> ks, double tilt) {
    std::vector<std::filesystem::path> files;
    std::vector<int> order = ks;
    std::reverse(order.begin(), order.end());
    std::rotate(order.begin(), order.begin() + 2, order.end());  // shuffled
    for (int k : order) {
      const eg::Vec3 p{o[0] + k * dz * n[0] + k * tilt * row[0], o[1] + k * dz * n[1] + k * tilt * row[1], o[2] + k * dz * n[2] + k * tilt * row[2]};
      const auto f = dir / (name + "-" + std::to_string(k) + ".dcm");
      write_slice(f, explicit_vr, explicit_vr ? "1.2.840.10008.1.2.1" : "1.2.840.10008.1.2", p, row, col, k + 1, stored(k));
      files.push_back(f);
    }
    return files;
  };

  for (int ex = 1; ex >= 0; ex--) {
    const auto files = series(ex ? "explicit" : "implicit", ex != 0, {0, 1, 2, 3, 4, 5, 6}, 0);
    md::SeriesReport r;
    const auto im = md::read_series(files, &r);
    double worst = 0;
    bool values = true;
    for (std::size_t k = 0; k < 7; k++)
      for (std::size_t j = 0; j < 6; j++)
        for (std::size_t i = 0; i < 5; i++) {
          const auto w = im.geometry.world(double(i), double(j), double(k));
          for (int a = 0; a < 3; a++) worst = std::max(worst, std::fabs(w[a] - (o[a] + double(i) * 0.6 * row[a] + double(j) * 0.8 * col[a] + double(k) * dz * n[a])));
          values = values && im.at(i, j, k) == float(1000 + 100 * int(k) + int(i + 5 * j) - 1024);
        }
    const std::string vr = ex ? "explicit VR" : "implicit VR";
    s.check(vr + ": every voxel where the series puts it, Pixel Spacing read row spacing first", worst, 1e-6, "mm");
    s.require(vr + ": values rescaled to Hounsfield units, slices in order", values && r.even && r.unsheared && r.slices == 7);
    s.check(vr + ": the spacing from the positions (2.5 mm), not Slice Thickness (3 mm)", r.spacing - 2.5, 1e-6, "mm");
  }
  {
    std::vector<md::Slice> sl;
    for (const auto& f : series("gap", true, {0, 1, 2, 4, 5}, 0)) sl.push_back(md::read_slice(f));
    const auto r = md::sort_series(sl);
    bool refused = false;
    try { md::series_geometry(sl, r); } catch (const std::invalid_argument&) { refused = true; }
    s.require("a missing slice: a gap after the third, the series refused as uneven", !r.even && r.gaps == std::vector<std::size_t>{2} && refused);
    std::vector<md::Slice> tl;
    for (const auto& f : series("tilt", true, {0, 1, 2, 3}, 0.4)) tl.push_back(md::read_slice(f));
    const auto rt = md::sort_series(tl);
    s.check("a tilted gantry: the in-plane step reported as shear", rt.shear - 0.4, 1e-6, "mm");
    s.require("and the series refused as sheared", !rt.unsheared);
  }
  {
    const auto f = dir / "compressed.dcm";
    write_slice(f, true, "1.2.840.10008.1.2.4.90", o, row, col, 1, stored(0));
    std::string what;
    try { md::read_slice(f); } catch (const std::runtime_error& e) { what = e.what(); }
    s.require("a compressed transfer syntax refused with its UID", what.find("1.2.840.10008.1.2.4.90") != std::string::npos, what);
  }
  {
    eg::Geometry ras;
    ras.size = {4, 3, 2};
    ras.spacing = {0.7, 0.9, 2.0};
    ras.origin = {10, -20, 30};
    for (int a = 0; a < 3; a++)
      for (int b = 0; b < 3; b++) ras.direction[a][b] = R[a][b];
    const auto lps = md::to_lps(ras);
    double worst = 0;
    for (std::size_t v = 0; v < ras.count(); v++) {
      const auto ijk = ras.unlinear(v);
      const auto p = ras.world(double(ijk[0]), double(ijk[1]), double(ijk[2])), q = lps.world(double(ijk[0]), double(ijk[1]), double(ijk[2]));
      worst = std::max({worst, std::fabs(q[0] + p[0]), std::fabs(q[1] + p[1]), std::fabs(q[2] - p[2])});
    }
    s.check("NIfTI's RAS carried to LPS: x and y negated at every voxel", worst, 1e-12, "mm");
  }
  std::filesystem::remove_all(dir);
  return s.done();
}
