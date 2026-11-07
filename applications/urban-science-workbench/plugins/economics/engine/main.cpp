// The engine of Quantitative economics (ADR 0080): the toolbox of quantitative economics (api/cpp/science/toolbox.hpp), asked one
// request at a time, a request a line and an answer a line ({"ok": true, "result": …} or {"ok": false, "error": …}): with --serve until
// its input ends, as the workbench's server keeps it; otherwise the one request given as an argument. --data names the directory of
// the data library's products, the only files it reads.
//
//   tool    {tool, columns?, roles, params, by_row?, places?, each?}   a tool of the toolbox on the columns the page sends
//   status  {}                                                         the engine answers
//
// The page holds the tables and sends only the columns a tool asks for. The tools on the world's trade read their product, whose
// path this engine gives them.

#include <science/json.hpp>
#include <science/toolbox.hpp>

#include <cstring>
#include <iostream>
#include <string>

namespace {

using science::Json;

std::string& data_dir() {
  static std::string dir = "data/cache/products";
  return dir;
}

Json answer(const Json& r) {
  const std::string op = r["op"].text();
  if (op == "status") return Json::object({{"ok", true}});
  if (op != "tool") throw std::invalid_argument("no operation \"" + op + "\"");
  const Json* tool = r.find("tool");
  if (tool && tool->is_string() && tool->as_string().rfind("trade-", 0) == 0) {
    Json request = r;
    request.set("flows", data_dir() + "/trade-flows.json");
    request.set("sizes", data_dir() + "/tables/economies.json");
    return science::toolbox::run(request);
  }
  return science::toolbox::run(r);
}

std::string respond(const std::string& line) {
  Json reply = Json::object();
  try {
    const Json r = Json::parse(line);
    if (r.has("id")) reply.set("id", r["id"]);
    reply.set("ok", true);
    reply.set("result", answer(r));
  } catch (const std::exception& e) {
    reply = Json::object({{"ok", false}, {"error", e.what()}});
  }
  return reply.dump();
}

}  // namespace

int main(int argc, char** argv) {
  std::ios::sync_with_stdio(false);
  bool serve = false;
  std::string request;
  for (int i = 1; i < argc; i++) {
    if (std::strcmp(argv[i], "--serve") == 0) serve = true;
    else if (std::strcmp(argv[i], "--data") == 0 && i + 1 < argc) data_dir() = argv[++i];
    else request = argv[i];
  }
  if (serve) {
    std::string line;
    while (std::getline(std::cin, line)) {
      if (line.empty()) continue;
      std::cout << respond(line) << "\n" << std::flush;
    }
    return 0;
  }
  if (request.empty()) std::getline(std::cin, request);
  std::cout << respond(request) << "\n";
  return 0;
}
