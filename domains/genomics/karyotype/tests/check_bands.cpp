// Checks of the ideogram (bands.hpp), held to its source files in the data library's cache
// (ncbi-ideogram): NCBI's GRCh38 levels of 400, 550 and 850 bands read band for band, each
// chromosome covered from its first base to its last without gaps, the bands in order outward from
// the centromere on each arm, the coarser levels the prefix classes of the finer with the same
// coordinates, NCBI's 850 level identical to UCSC's hg38 cytoBand, GRCh37 read from NCBI's files
// whatever their convention of starts and identical to UCSC's hg19; the bands holding published
// coordinates named as the publications name them; faults of a file reported with their lines.

#include <genomics/karyotype/bands.hpp>

#include <calx/verification/check.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

using namespace genomics::karyotype;

namespace {

std::string read_file(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  std::stringstream s;
  s << f.rdbuf();
  return s.str();
}

std::size_t data_lines(const std::string& text) {
  std::size_t n = 0;
  std::istringstream in(text);
  std::string line;
  while (std::getline(in, line))
    if (!line.empty() && line[0] != '#' && line != "\r") n++;
  return n;
}

bool same_bands(const Ideogram& a, const Ideogram& b, bool stains, std::string& detail) {
  if (a.bands.size() != b.bands.size()) { detail = std::to_string(a.bands.size()) + " vs " + std::to_string(b.bands.size()); return false; }
  for (std::size_t i = 0; i < a.bands.size(); i++) {
    const Band &x = a.bands[i], &y = b.bands[i];
    if (x.chromosome != y.chromosome || x.arm != y.arm || x.digits != y.digits || x.start != y.start || x.end != y.end || (stains && x.stain != y.stain)) {
      detail = x.full_name() + " [" + std::to_string(x.start) + "," + std::to_string(x.end) + ") " + std::string(stain_name(x.stain)) + " vs " + y.full_name() +
               " [" + std::to_string(y.start) + "," + std::to_string(y.end) + ") " + std::string(stain_name(y.stain));
      return false;
    }
  }
  return true;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- names --------------------------------------------------------------------------------------
  {
    s.require("a band's word of digits is its name without the point: 36.33 is 3633, 13 is 13, 13.3 is 133",
              band_digits("36.33") == std::string("3633") && band_digits("13") == std::string("13") && band_digits("13.3") == std::string("133"));
    s.require("names that are not bands are refused: 3.6, 36., 36.3.3, p13, an empty word",
              !band_digits("3.6") && !band_digits("36.") && !band_digits("36.3.3") && !band_digits("p13") && !band_digits(""));
    s.require("written back with the point after the second digit: p36.33, q21, q11.2",
              band_label('p', "3633") == "p36.33" && band_label('q', "21") == "q21" && band_label('q', "112") == "q11.2");
    s.require("chromosomes 1-22, X, Y are 0-23 and chr-prefixed names are read; 0, 23, 01 are not chromosomes",
              chromosome_index("1") == 0 && chromosome_index("22") == 21 && chromosome_index("X") == chr_x && chromosome_index("chrY") == chr_y &&
                  chromosome_index("0") < 0 && chromosome_index("23") < 0 && chromosome_index("01") < 0);
  }

  // ---- faults of a file -------------------------------------------------------------------------------
  {
    const std::string text =
        "#chromosome\tarm\tband\tiscn_start\tiscn_stop\tbp_start\tbp_stop\tstain\tdensity\n"
        "1\tp\t36.33\t0\t100\t1\t2300000\tgneg\t\n"
        "1\tp\t36.32\t100\t244\t2300001\t5300000\tgpos\t25\n"
        "1\tr\t36.31\t244\t344\t5300001\t7100000\tgneg\t\n"
        "1\tp\t36.31\t244\t344\t5300001\t7100000\tgblue\t\n"
        "23\tp\t11\t0\t1\t1\t2\tgneg\t\n"
        "1\tp\t36.23\t344\t459\t7200001\t9100000\tgpos\t25\n"
        "1\tp\t36\n";
    std::vector<Problem> problems;
    const Ideogram g = read_ideogram(text, 850, &problems);
    std::string lines;
    for (const auto& p : problems) lines += std::to_string(p.line) + ": " + p.message + "; ";
    s.require("faults reported with their lines: an arm r, a stain gblue, a chromosome 23, fewer than eight columns, and the gap they leave",
              problems.size() == 5 && problems[0].line == 4 && problems[1].line == 5 && problems[2].line == 6 && problems[3].line == 8 &&
                  problems[4].message.find("gap") != std::string::npos,
              lines);
    s.require("the bands that could be read are kept: three of chromosome 1", g.bands.size() == 3 && g.end_of(0) - g.begin_of(0) == 3);
  }

  const char* env = std::getenv("SCIENCE_CACHE");
  const std::filesystem::path cache = (env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE)) / "ncbi-ideogram";
  const auto file = [&](const std::string& name) { return cache / name; };
  if (!std::filesystem::exists(file("ideogram_9606_GCF_000001305.16_850_V1.tsv"))) {
    s.skip("the ideogram is not in the data cache (python data/acquisition/fetch.py ncbi-ideogram): its checks are not run");
    return s.done();
  }

  // ---- GRCh38 at three levels ----------------------------------------------------------------------------
  Ideogram levels[3];
  const int level_names[3] = {400, 550, 850};
  for (int k = 0; k < 3; k++) {
    const std::string text = read_file(file("ideogram_9606_GCF_000001305.16_" + std::to_string(level_names[k]) + "_V1.tsv"));
    std::vector<Problem> problems;
    levels[k] = read_ideogram(text, level_names[k], &problems);
    const Ideogram& g = levels[k];
    s.require("GRCh38 at " + std::to_string(level_names[k]) + " bands: every one of the file's " + std::to_string(data_lines(text)) +
                  " bands read, with no problem",
              g.bands.size() == data_lines(text) && problems.empty(), problems.empty() ? "" : problems[0].message);
    bool cover = true, order = true, acen = true;
    std::string where;
    for (int c = 0; c < chromosome_count; c++) {
      if (!g.has(c)) { cover = false; where = chromosome_name(c) + " missing"; continue; }
      std::uint64_t at = 0;
      for (std::size_t b = g.begin_of(c); b < g.end_of(c); b++) {
        const Band& x = g.bands[b];
        if (x.start != at || x.end <= x.start) { cover = false; where = x.full_name(); }
        at = x.end;
        if (b + 1 < g.end_of(c)) {
          const Band& y = g.bands[b + 1];
          // Outward from the centromere: on p the words fall toward the centromere, on q they rise.
          if (x.arm == 'p' && y.arm == 'p' && !(x.digits > y.digits)) { order = false; where = x.full_name() + " " + y.full_name(); }
          if (x.arm == 'q' && y.arm == 'q' && !(x.digits < y.digits)) { order = false; where = x.full_name() + " " + y.full_name(); }
          if (x.arm == 'q' && y.arm == 'p') { order = false; where = x.full_name(); }
          if (x.arm == 'p' && y.arm == 'q' && (x.stain != Stain::acen || y.stain != Stain::acen || g.centromere[c] != y.start)) { acen = false; where = x.full_name(); }
        }
      }
      if (at != g.length[c]) cover = false;
    }
    s.require("GRCh38 at " + std::to_string(level_names[k]) + ": each of the 24 chromosomes covered from its first base to its last without gaps or overlaps", cover, where);
    s.require("GRCh38 at " + std::to_string(level_names[k]) + ": the bands run from pter to qter, numbered outward from the centromere on each arm", order, where);
    s.require("GRCh38 at " + std::to_string(level_names[k]) + ": the centromere lies between two acen bands, the last of p and the first of q", acen, where);
  }
  const Ideogram& g = levels[2];
  {
    // The coarser levels are prefix classes of the finer, with their coordinates.
    for (int k = 0; k < 2; k++) {
      bool ok = true;
      std::string where;
      for (const Band& b : levels[k].bands) {
        const auto r = g.find(b.chromosome, b.arm, b.digits);
        if (!r || r->start != b.start || r->end != b.end) { ok = false; where = b.full_name(); break; }
      }
      s.require("every band of the " + std::to_string(level_names[k]) + " level is the run of 850-level bands whose names it begins, with the same coordinates", ok, where);
    }
    bool runs = true;
    std::string where;
    for (int c = 0; c < chromosome_count && runs; c++)
      for (std::size_t b = g.begin_of(c); b < g.end_of(c) && runs; b++)
        for (std::size_t k = 1; k <= g.bands[b].digits.size(); k++) {
          const std::string prefix = g.bands[b].digits.substr(0, k);
          const auto r = g.find(c, g.bands[b].arm, prefix);
          for (std::size_t i = r->first; i <= r->last; i++)
            if (g.bands[i].arm != g.bands[b].arm || g.bands[i].digits.compare(0, k, prefix) != 0) { runs = false; where = chromosome_name(c) + band_label(g.bands[b].arm, prefix); }
        }
    s.require("every prefix of a band's name (region, band, sub-band) names one run of consecutive bands: the hierarchy is the prefix order", runs, where);
    s.require("the arms: 1p from 0 to the centromere at 123,400,000, and 1q after it (NCBI's 1p11.1 ends at base 123,400,000)",
              g.arm(0, 'p')->start == 0 && g.arm(0, 'p')->end == g.centromere[0] && g.arm(0, 'q')->start == g.centromere[0] && g.centromere[0] == 123400000);
  }

  // ---- NCBI's 850 level is UCSC's hg38 table -------------------------------------------------------------
  {
    std::vector<Problem> problems;
    const Ideogram u = read_cytoband(read_file(file("ucsc-hg38-cytoBand.txt")), 850, &problems);
    std::string detail;
    s.require("UCSC's hg38 cytoBand holds NCBI's 862 bands of the 850 level with the same coordinates and stains", same_bands(g, u, true, detail) && problems.empty(), detail);
  }

  // ---- positions and names --------------------------------------------------------------------------------
  {
    s.require("base 0 of chromosome 1 lies in 1p36.33, and its last base in 1q44", g.bands[g.band_at(0, 0)].full_name() == "1p36.33" &&
                                                                                       g.bands[g.band_at(0, g.length[0] - 1)].full_name() == "1q44");
    // Conlin et al. (2025), Clin. Chem. 71, 45-55, Table 1 and text: arr[GRCh38] 1q21.1q21.2(147,105,904_147,917,509)x1 and
    // arr[GRCh38] 7q11.23(73,265,388_74,726,130)x1, the bands of published GRCh38 coordinates.
    s.require("the bands of 1:147,105,904-147,917,509 (GRCh38) are 1q21.1q21.2, as Conlin et al. write them", g.range_name(0, 147105903, 147917509) == "q21.1q21.2",
              g.range_name(0, 147105903, 147917509));
    s.require("the bands of 7:73,265,388-74,726,130 (GRCh38) are 7q11.23, as Conlin et al. write them", g.range_name(6, 73265387, 74726130) == "q11.23",
              g.range_name(6, 73265387, 74726130));
    s.require("the smallest band holding 1q21.1-1q21.2 is 1q21, and none holds an interval across the centromere",
              g.smallest_band(0, 147105903, 147917509) == std::string("q21") && !g.smallest_band(0, 100000000, 150000000));
    const auto r = g.find(0, "q2");
    s.require("region 1q2 runs from 1q21.1 to 1q25.3", r && g.bands[r->first].full_name() == "1q21.1" && g.bands[r->last].full_name() == "1q25.3");
    s.require("no band 5q37 and no band 21p14", !g.find(4, "q37") && !g.find(20, "p14"));
  }

  // ---- GRCh37 --------------------------------------------------------------------------------------------
  if (std::filesystem::exists(file("ucsc-hg19-cytoBand.txt"))) {
    std::vector<Problem> pu, p37;
    const Ideogram u19 = read_cytoband(read_file(file("ucsc-hg19-cytoBand.txt")), 850, &pu, "GRCh37");
    const Ideogram n19 = read_ideogram(read_file(file("ideogram_9606_GCF_000001305.13_850_V1.tsv")), 850, &p37, "GRCh37");
    std::string detail;
    s.require("GRCh37: NCBI's 850 level, whose starts are the previous band's stops, read as UCSC's hg19 table holds it, band for band (stains aside: NCBI's carry no density)",
              same_bands(u19, n19, false, detail) && pu.empty() && p37.empty(), detail);
    bool gaps = true;
    for (int lv : {400, 550}) {
      std::vector<Problem> pr;
      const Ideogram x = read_ideogram(read_file(file("ideogram_9606_GCF_000001305.13_" + std::to_string(lv) + "_V1.tsv")), lv, &pr, "GRCh37");
      gaps = gaps && pr.empty() && x.bands.size() > 300;
      for (const Band& b : x.bands) {
        const auto r = u19.find(b.chromosome, b.arm, b.digits);
        if (!r || r->start != b.start || r->end != b.end) { gaps = false; detail = b.full_name(); break; }
      }
    }
    s.require("GRCh37 at 400 and 550 bands, read whatever the convention of their starts, are prefix classes of the 850 level with its coordinates", gaps, detail);
    s.require("chromosome 1 of GRCh37 ends at 249,250,621 and that of GRCh38 at 248,956,422 (the files' own last bands)",
              u19.length[0] == 249250621 && g.length[0] == 248956422);
    // Ordulu et al. (2014): arr[hg19] 12p12.1(21,340,001–23,500,000)×1.
    s.require("the bands of 12:21,340,001-23,500,000 on GRCh37 are 12p12.1, as Ordulu et al. write them", u19.range_name(11, 21340000, 23500000) == "p12.1",
              u19.range_name(11, 21340000, 23500000));
  } else {
    s.skip("UCSC's hg19 cytoBand is not in the cache: GRCh37 is not checked");
  }
  return s.done();
}
