// Checks of the loss-gain-fusion model against CytoGPS's own (Abrams et al. 2019, Bioinformatics
// 35, 5365-5366): the six karyotypes of RCytoGPS 1.2.14's examples, each clone's bands of loss, gain
// and fusion as the CytoGPS web service computed them on 2020-04-10, and the detailed forms of their
// derivatives; CytoGPS's bands (868, a few named at another resolution than NCBI's 850 level, and
// p10/q10 for the centromere) are compared with NCBI's by the prefix order of band names.

#include <genomics/karyotype/interpret.hpp>

#include <calx/verification/check.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace genomics::karyotype;

namespace {

std::string read_file(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  std::stringstream s;
  s << f.rdbuf();
  return s.str();
}

std::vector<std::string> split(const std::string& s, char sep) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string x;
  while (std::getline(ss, x, sep))
    if (!x.empty()) out.push_back(x);
  return out;
}

bool centromeric(const std::string& name) {
  return name.size() > 3 && (name.substr(name.size() - 3) == "p10" || name.substr(name.size() - 3) == "q10");
}

// CytoGPS's band names as sets of NCBI's 850-level bands: a band whose name begins the other's, or
// whose name the other begins.
std::set<std::size_t> to_ncbi(const std::vector<std::string>& names, const Ideogram& g) {
  std::set<std::size_t> out;
  for (const auto& n : names) {
    if (centromeric(n)) continue;
    const std::size_t a = n.find_first_of("pq");
    const int c = chromosome_index(n.substr(0, a));
    const auto d = band_digits(n.substr(a + 1));
    if (c < 0 || !d) continue;
    for (std::size_t b = g.begin_of(c); b < g.end_of(c); b++) {
      const Band& x = g.bands[b];
      if (x.arm != n[a]) continue;
      if (x.digits.compare(0, d->size(), *d) == 0 || d->compare(0, x.digits.size(), x.digits) == 0) out.insert(b);
    }
  }
  return out;
}

std::string names(const std::set<std::size_t>& s, const Ideogram& g) {
  std::string out;
  for (auto b : s) out += (out.empty() ? "" : ",") + g.bands[b].full_name();
  return out;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const char* env = std::getenv("SCIENCE_CACHE");
  const std::filesystem::path path = (env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE)) / "ncbi-ideogram" / "ideogram_9606_GCF_000001305.16_850_V1.tsv";
  if (!std::filesystem::exists(path)) {
    s.skip("the ideogram is not in the data cache (python data/acquisition/fetch.py ncbi-ideogram): CytoGPS's model is not compared");
    return s.done();
  }
  const Ideogram g = read_ideogram(read_file(path), 850);

  std::ifstream f(std::filesystem::path(GENOMICS_FIXTURES) / "karyotype-cytogps.tsv");
  std::string line;
  std::size_t clones = 0;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::vector<std::string> cols;
    std::stringstream ss(line);
    std::string c;
    while (std::getline(ss, c, '\t')) cols.push_back(c);
    while (cols.size() < 7) cols.push_back("");
    const std::string& karyotype = cols[0];
    const Parsed p = parse(karyotype, Edition::iscn2016);
    if (cols[1] != "Success") {
      s.require(karyotype + ": refused by CytoGPS (" + cols[1] + ") and here: " + p.errors(), !p.ok());
      continue;
    }
    const Interpretation in = interpret(p.description, g);
    const std::size_t k = static_cast<std::size_t>(std::stoi(cols[2]));
    if (!p.ok() || k >= in.lines.size()) {
      s.require(karyotype + ": accepted by CytoGPS and here", false, p.errors());
      continue;
    }
    clones++;
    const LineInterpretation& L = in.lines[k];
    const auto loss = to_ncbi(split(cols[3], ','), g), gain = to_ncbi(split(cols[4], ','), g), fusion = to_ncbi(split(cols[5], ','), g);
    const std::set<std::size_t> my_loss(L.lgf.loss.begin(), L.lgf.loss.end()), my_gain(L.lgf.gain.begin(), L.lgf.gain.end()),
        my_fusion(L.lgf.fusion.begin(), L.lgf.fusion.end());
    std::set<std::string> cen, my_cen;
    for (const auto& n : split(cols[5], ','))
      if (centromeric(n)) cen.insert(n);
    for (const auto& n : L.lgf.fusion_labels)
      if (centromeric(n)) my_cen.insert(n);
    const std::string tag = karyotype + (in.lines.size() > 1 ? " (clone " + std::to_string(k + 1) + ")" : "");
    s.require(tag + ": the bands lost are CytoGPS's (" + std::to_string(loss.size()) + ")", my_loss == loss, names(my_loss, g) + " vs " + names(loss, g));
    s.require(tag + ": the bands gained are CytoGPS's (" + std::to_string(gain.size()) + ")", my_gain == gain, names(my_gain, g) + " vs " + names(gain, g));
    s.require(tag + ": the bands of fusion are CytoGPS's (" + std::to_string(fusion.size() + cen.size()) + ")", my_fusion == fusion && my_cen == cen,
              names(my_fusion, g) + " vs " + names(fusion, g));
    // The detailed forms CytoGPS writes for its derivatives, as der(2;13)(2qter->2q10::13q10->13qter).
    const auto want = split(cols[6], '|');
    if (!want.empty()) {
      std::vector<std::string> mine;
      LongFormOptions o;
      o.unicode = false;
      for (const auto& ch : L.chromosomes)
        if (ch.abnormal && ch.source >= 0 && L.abnormalities[static_cast<std::size_t>(ch.source)].type == Type::der) {
          const Abnormality& a = L.abnormalities[static_cast<std::size_t>(ch.source)];
          std::string head = "der(";
          for (std::size_t i = 0; i < a.chromosomes.size(); i++) head += (i ? ";" : "") + a.chromosomes[i].name();
          mine.push_back(head + ")(" + long_form(ch, o) + ")");
        }
      std::string got;
      for (const auto& m : mine) got += m + " ";
      s.require(tag + ": the derivatives written as CytoGPS writes them, " + cols[6], mine == want, got);
    }
  }
  s.require("every clone CytoGPS computed was compared (" + std::to_string(clones) + ")", clones == 10);
  return s.done();
}
