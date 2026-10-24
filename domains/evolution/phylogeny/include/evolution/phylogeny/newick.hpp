// The Newick format: trees read from text and written to it.
//
// The grammar is Felsenstein's specification ("The Newick tree format", the PHYLIP
// documentation, after the meeting at Newick's restaurant in 1986; Felsenstein 2004,
// Inferring Phylogenies, on the Newick format):
//
//   tree      ::= subtree [":" length] ";"
//   subtree   ::= label | "(" branch {"," branch} ")" [label]
//   branch    ::= subtree [":" length]
//   label     ::= unquoted | "'" {char | "''"} "'"
//
// An unquoted label is a run of characters other than blanks and ( ) [ ] ' : ; , in which an
// underscore stands for a blank; a quoted label may hold any character, a quote written twice.
// Comments in square brackets may stand wherever a blank may, and nest; a comment [&R] or [&U]
// before a tree says whether it is rooted (the convention of PAUP* and BEAST). Blanks and line
// breaks between tokens are ignored. An internal node's label that is a number is read as the
// support of the branch above it (the convention of RAxML, IQ-TREE and PhyML for bootstrap
// values), unless the options say otherwise.
//
// Reading never throws on content: it returns the trees it could read and a problem for each
// fault, with its line and column; a tree with a fault is skipped to its ';'. Without [&R] or
// [&U], a tree whose root has two children is read as rooted, one whose root has three or more
// as unrooted. Writing quotes a label when it must and gives lengths to `precision`
// significant digits, so that reading back gives the same tree. Both are O(length of text).

#pragma once

#include <evolution/phylogeny/tree.hpp>

#include <genomics/sequence/fasta.hpp>

#include <charconv>
#include <cmath>
#include <cstdio>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace evolution::phylogeny {

using genomics::sequence::Problem;

struct NewickOptions {
  bool underscores_are_blanks = true;      // Felsenstein's rule for unquoted labels
  bool numeric_labels_are_supports = true;  // an internal label that is a number is a support
};

struct NewickTrees {
  std::vector<Tree> trees;
  std::vector<Problem> problems;
};

namespace detail {

inline bool newick_special(char c) {
  return c == '(' || c == ')' || c == '[' || c == ']' || c == '\'' || c == ':' || c == ';' || c == ',' || c == ' ' || c == '\t' ||
         c == '\n' || c == '\r';
}

inline bool parse_number(std::string_view s, double& out) {
  if (s.empty()) return false;
  if (s.front() == '+') s.remove_prefix(1);
  const auto r = std::from_chars(s.data(), s.data() + s.size(), out);
  return r.ec == std::errc() && r.ptr == s.data() + s.size() && std::isfinite(out);
}

class NewickReader {
 public:
  NewickReader(std::string_view text, const NewickOptions& o) : s_(text), o_(o) {}

  NewickTrees read() {
    NewickTrees out;
    for (;;) {
      std::string open;
      rooted_mark_ = 0;
      if (!skip_blank(&open)) {
        out.problems.push_back({line_of(fault_), "column " + std::to_string(column_of(fault_)) + ": " + open});
        break;
      }
      if (at_end()) break;
      const std::size_t start = i_;
      Tree t;
      std::string message;
      if (tree(t, message)) {
        out.trees.push_back(std::move(t));
      } else {
        out.problems.push_back({line_of(fault_), "column " + std::to_string(column_of(fault_)) + ": " + message});
        // Skip to the end of this tree.
        while (!at_end() && s_[i_] != ';') i_++;
        if (!at_end()) i_++;
        if (i_ == start) i_++;
      }
    }
    return out;
  }

 private:
  std::string_view s_;
  NewickOptions o_;
  std::size_t i_ = 0, fault_ = 0;
  int rooted_mark_ = 0;  // +1 for [&R], −1 for [&U]

  bool at_end() const { return i_ >= s_.size(); }
  std::size_t line_of(std::size_t p) const {
    std::size_t l = 1;
    for (std::size_t k = 0; k < p && k < s_.size(); k++) l += s_[k] == '\n';
    return l;
  }
  std::size_t column_of(std::size_t p) const {
    std::size_t c = 1;
    for (std::size_t k = 0; k < p && k < s_.size(); k++) c = s_[k] == '\n' ? 1 : c + 1;
    return c;
  }
  bool fail(std::string& message, std::string m) {
    fault_ = i_;
    message = std::move(m);
    return false;
  }

  // Blanks and comments; a comment [&R] or [&U] marks the tree rooted or unrooted.
  bool skip_blank(std::string* message = nullptr) {
    for (;;) {
      while (!at_end() && (s_[i_] == ' ' || s_[i_] == '\t' || s_[i_] == '\n' || s_[i_] == '\r')) i_++;
      if (at_end() || s_[i_] != '[') return true;
      const std::size_t open = i_;
      int depth = 0;
      std::size_t body = i_ + 1;
      do {
        if (s_[i_] == '[') depth++;
        else if (s_[i_] == ']') depth--;
        i_++;
      } while (!at_end() && depth > 0);
      if (depth > 0) {
        if (message) {
          fault_ = open;
          *message = "a comment is not closed";
        }
        return false;
      }
      const std::string_view c = s_.substr(body, i_ - 1 - body);
      if (c == "&R" || c == "&r") rooted_mark_ = 1;
      if (c == "&U" || c == "&u") rooted_mark_ = -1;
    }
  }

  bool label(std::string& out, std::string& message) {
    out.clear();
    if (!at_end() && s_[i_] == '\'') {
      i_++;
      for (;;) {
        if (at_end()) return fail(message, "a quoted label is not closed");
        if (s_[i_] == '\'') {
          if (i_ + 1 < s_.size() && s_[i_ + 1] == '\'') {
            out += '\'';
            i_ += 2;
            continue;
          }
          i_++;
          return true;
        }
        out += s_[i_++];
      }
    }
    while (!at_end() && !newick_special(s_[i_])) {
      out += (o_.underscores_are_blanks && s_[i_] == '_') ? ' ' : s_[i_];
      i_++;
    }
    return true;
  }

  bool tree(Tree& t, std::string& message) {
    if (!skip_blank(&message)) return false;
    int cur = t.add(-1);
    std::vector<bool> named;  // whether a node has its label (by index)
    named.push_back(false);
    std::vector<bool> closed{false};  // whether a node's subtree has been closed by ')'
    for (;;) {
      if (!skip_blank(&message)) return false;
      if (at_end()) return fail(message, "the tree has no ';'");
      const char c = s_[i_];
      if (c == '(') {
        if (named[static_cast<std::size_t>(cur)] || !t.is_leaf(cur) || closed[static_cast<std::size_t>(cur)])
          return fail(message, "'(' after a label or a subtree");
        i_++;
        cur = t.add(cur);
        named.push_back(false);
        closed.push_back(false);
      } else if (c == ',') {
        const int p = t[cur].parent;
        if (p < 0) return fail(message, "',' outside every '('");
        i_++;
        cur = t.add(p);
        named.push_back(false);
        closed.push_back(false);
      } else if (c == ')') {
        const int p = t[cur].parent;
        if (p < 0) return fail(message, "')' without its '('");
        i_++;
        cur = p;
        closed[static_cast<std::size_t>(cur)] = true;
      } else if (c == ':') {
        i_++;
        if (!skip_blank(&message)) return false;
        const std::size_t b = i_;
        while (!at_end() && !newick_special(s_[i_])) i_++;
        double x = 0;
        if (!parse_number(s_.substr(b, i_ - b), x)) {
          i_ = b;
          return fail(message, "a branch length is not a number");
        }
        if (!std::isnan(t[cur].length)) return fail(message, "a branch has two lengths");
        t[cur].length = x;
      } else if (c == ';') {
        if (t[cur].parent >= 0) return fail(message, "a '(' is not closed");
        i_++;
        break;
      } else if (c == ']') {
        return fail(message, "']' without its '['");
      } else {
        if (named[static_cast<std::size_t>(cur)]) return fail(message, "a node has two labels");
        std::string l;
        if (!label(l, message)) return false;
        named[static_cast<std::size_t>(cur)] = true;
        t[cur].label = std::move(l);
      }
    }
    // Numeric internal labels are supports of the branch above.
    if (o_.numeric_labels_are_supports)
      for (std::size_t v = 0; v < t.size(); v++) {
        auto& n = t.nodes[v];
        double x = 0;
        if (!n.children.empty() && static_cast<int>(v) != t.root && parse_number(n.label, x)) {
          n.support = x;
          n.label.clear();
        }
      }
    t.rooted = rooted_mark_ > 0 || (rooted_mark_ == 0 && t[t.root].children.size() == 2);
    return true;
  }
};

}  // namespace detail

inline NewickTrees read_newick(std::string_view text, const NewickOptions& options = {}) {
  return detail::NewickReader(text, options).read();
}

struct NewickWriting {
  bool lengths = true;
  bool supports = true;          // an internal node's support, written as its label when it has none
  bool rooted_mark = false;      // [&R] or [&U] before the tree
  int precision = 12;            // significant digits of lengths
  bool underscores_are_blanks = true;
};

namespace detail {

inline std::string newick_label(const std::string& l, bool underscores) {
  bool quote = false, blank = false;
  for (char c : l) {
    if (c == ' ') blank = true;
    else if (newick_special(c) || (underscores && c == '_')) quote = true;
  }
  if (blank && !underscores) quote = true;
  if (!quote) {
    std::string out = l;
    for (char& c : out)
      if (c == ' ') c = '_';
    return out;
  }
  std::string out = "'";
  for (char c : l) {
    if (c == '\'') out += "''";
    else out += c;
  }
  return out + "'";
}

inline std::string newick_number(double x, int precision) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.*g", precision, x);
  return buf;
}

}  // namespace detail

inline std::string write_newick(const Tree& t, const NewickWriting& w = {}) {
  std::string out;
  if (t.root < 0) return ";";
  if (w.rooted_mark) out += t.rooted ? "[&R] " : "[&U] ";
  // Iterative: a stack of (node, next child to write).
  std::vector<std::pair<int, std::size_t>> stack{{t.root, 0}};
  while (!stack.empty()) {
    auto& [v, k] = stack.back();
    const auto& n = t[v];
    if (k == 0 && !n.children.empty()) out += '(';
    if (k < n.children.size()) {
      if (k > 0) out += ',';
      const int c = n.children[k++];
      stack.push_back({c, 0});
      continue;
    }
    if (!n.children.empty()) out += ')';
    if (!n.label.empty()) out += detail::newick_label(n.label, w.underscores_are_blanks);
    else if (w.supports && !n.children.empty() && !std::isnan(n.support)) out += detail::newick_number(n.support, w.precision);
    if (w.lengths && v != t.root && !std::isnan(n.length)) out += ':' + detail::newick_number(n.length, w.precision);
    stack.pop_back();
  }
  return out + ';';
}

}  // namespace evolution::phylogeny
