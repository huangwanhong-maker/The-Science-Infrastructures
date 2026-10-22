// The HTTP of the API's clients in C and C++ (api/README.md §3.5): a client of
// HTTP/1.1 over TCP, on the standard library and the platform's sockets (POSIX;
// Winsock on Windows), with what the API needs and no more.
//
// - Plain http:// only. HTTPS needs a TLS library, which this header does not
//   bring. Across a network that is not trusted, put the service behind a TLS
//   proxy and run the program where the hop to the proxy is trusted (the proxy's
//   own computer, or a tunnel), or use the Python client, whose standard
//   library speaks TLS. A token sent over plain HTTP crosses the network in the
//   clear.
// - One request per connection (Connection: close): a client holds no
//   connection between its calls, and so outlives a network that comes and
//   goes, as a device's does.
// - Bodies framed by Content-Length, by the chunked transfer coding (RFC 9112
//   §7.1), or by the end of the connection, decoded as they arrive, so that a
//   stream is read while it lasts. Redirects are not followed.
// - A time limit on connecting and on each wait for bytes (not on a whole
//   answer, which may be long, as a large grid is); names resolved by
//   getaddrinfo, IPv4 and IPv6.
//
// The events of a text/event-stream are parsed by EventStream, after the WHATWG
// HTML Standard, §9.2.6 "Interpreting an event stream".

#pragma once

#include "json.hpp"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <winsock2.h>
#include <ws2tcpip.h>
#ifdef _MSC_VER
#pragma comment(lib, "Ws2_32.lib")
#endif
#else
#include <fcntl.h>
#include <netdb.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <poll.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <unistd.h>
#endif

namespace science {

// An answer that is an error: its HTTP status, or 0 when no answer came (the
// service could not be reached, did not answer in time, broke the connection or
// answered what is not HTTP), its message, the address asked, and, for a run of
// a plug-in that failed, its job.
class ApiError : public std::runtime_error {
 public:
  ApiError(int status, std::string message, std::string url = {}, Json job = {})
      : std::runtime_error((status ? std::to_string(status) : std::string("unreachable")) + ": " + message +
                           (url.empty() ? std::string() : " (" + url + ")")),
        status(status), message(std::move(message)), url(std::move(url)), job(std::move(job)) {}
  int status;
  std::string message, url;
  Json job;
};

namespace http {

// A wait for bytes that ran out of time: told apart, so that a stream reconnects.
class Timeout : public ApiError {
 public:
  using ApiError::ApiError;
};

inline std::string seconds_text(double seconds) { return Json::number_text(seconds); }

// ---- addresses ----------------------------------------------------------------------------------

struct Url {
  std::string host, port = "80", path;  // path: empty, or from its first '/'
};

// The parts of an http:// address; std::invalid_argument for any other.
inline Url parse_url(const std::string& url) {
  auto refuse = [&](const std::string& why) { return std::invalid_argument(why + ": " + url); };
  if (url.compare(0, 8, "https://") == 0)
    throw refuse("this client speaks plain http:// only: reach an HTTPS service through a TLS proxy on a trusted hop, "
                 "or use the Python client (api/README.md §3.5)");
  const std::string scheme = "http://";
  if (url.compare(0, scheme.size(), scheme) != 0) throw refuse("an address of the API begins with http://");
  std::string rest = url.substr(scheme.size());
  std::size_t slash = rest.find_first_of("/?#");
  std::string authority = rest.substr(0, slash);
  Url u;
  if (slash != std::string::npos) u.path = rest.substr(slash);
  if (authority.find('@') != std::string::npos) throw refuse("an address with a user name is not taken: give the token instead");
  std::string port;
  if (!authority.empty() && authority[0] == '[') {  // an IPv6 literal, [::1]:8800
    std::size_t close = authority.find(']');
    if (close == std::string::npos) throw refuse("an IPv6 address without its closing bracket");
    u.host = authority.substr(1, close - 1);
    if (close + 1 < authority.size()) {
      if (authority[close + 1] != ':') throw refuse("an address that is not host:port");
      port = authority.substr(close + 2);
    }
  } else {
    std::size_t colon = authority.rfind(':');
    u.host = authority.substr(0, colon);
    if (colon != std::string::npos) port = authority.substr(colon + 1);
  }
  if (!port.empty()) {
    if (port.size() > 5 || !std::all_of(port.begin(), port.end(), [](unsigned char c) { return std::isdigit(c); }))
      throw refuse("the port is not a number");
    u.port = port;
  }
  if (u.host.empty()) throw refuse("an address without a host");
  return u;
}

// Percent-encoding of a path segment or a query's key or value: every byte but
// the unreserved characters of RFC 3986 §2.3, so that '+' in now+6h stays a plus.
inline std::string encode(std::string_view s) {
  static const char* hex = "0123456789ABCDEF";
  std::string out;
  out.reserve(s.size());
  for (unsigned char c : s) {
    if (std::isalnum(c) || c == '-' || c == '.' || c == '_' || c == '~') {
      out += static_cast<char>(c);
    } else {
      out += '%';
      out += hex[c >> 4];
      out += hex[c & 15];
    }
  }
  return out;
}

inline std::string lower(std::string s) {
  for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

inline std::string trim(std::string_view s) {
  std::size_t a = 0, b = s.size();
  while (a < b && (s[a] == ' ' || s[a] == '\t')) a++;
  while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t')) b--;
  return std::string(s.substr(a, b - a));
}

// ---- the bytes of an answer -------------------------------------------------------------------------

// The source of an answer's bytes: a socket, or, in the checks, a text.
class Input {
 public:
  virtual ~Input() = default;
  // At most n bytes into out, waiting at most `seconds` for the first; 0 at the
  // end of the input; Timeout when none came in time.
  virtual std::size_t read(char* out, std::size_t n, double seconds) = 0;
};

// A TCP connection, non-blocking underneath, every wait bounded in time.
class Socket : public Input {
 public:
  Socket() = default;
  Socket(const Socket&) = delete;
  Socket& operator=(const Socket&) = delete;
  ~Socket() override { close(); }

  void connect(const std::string& host, const std::string& port, double seconds) {
    close();
    startup();
    addrinfo hints{};
    hints.ai_family = AF_UNSPEC;
    hints.ai_socktype = SOCK_STREAM;
    hints.ai_protocol = IPPROTO_TCP;
    addrinfo* found = nullptr;
    int rc = ::getaddrinfo(host.c_str(), port.c_str(), &hints, &found);
    if (rc != 0 || !found) throw ApiError(0, "the name " + host + " cannot be resolved: " + resolve_text(rc));
    std::string why = "it has no address";
    for (addrinfo* a = found; a; a = a->ai_next) {
      Native s = ::socket(a->ai_family, a->ai_socktype, a->ai_protocol);
      if (s == NONE) {
        why = error_text(last_error());
        continue;
      }
      nonblocking(s);
#ifdef _WIN32
      int err = ::connect(s, a->ai_addr, static_cast<int>(a->ai_addrlen)) == 0 ? 0 : last_error();
#else
      int err = ::connect(s, a->ai_addr, a->ai_addrlen) == 0 ? 0 : last_error();
#endif
      if (err != 0 && in_progress(err)) {
        if (!wait(s, true, seconds)) {
          err = -1;
        } else {
          int e = 0;
          socklen_t len = sizeof e;
          ::getsockopt(s, SOL_SOCKET, SO_ERROR, reinterpret_cast<char*>(&e), &len);
          err = e;
        }
      }
      if (err == 0) {
        options(s);
        fd_ = s;
        ::freeaddrinfo(found);
        return;
      }
      why = err == -1 ? "no answer within " + seconds_text(seconds) + " s" : error_text(err);
      close_native(s);
    }
    ::freeaddrinfo(found);
    throw ApiError(0, "the service at " + host + ":" + port + " cannot be reached: " + why);
  }

  void write(std::string_view data, double seconds) {
    while (!data.empty()) {
      if (!wait(fd_, true, seconds)) throw Timeout(0, "the service took nothing for " + seconds_text(seconds) + " s");
      int n = static_cast<int>(std::min<std::size_t>(data.size(), 1 << 20));
      auto sent = ::send(fd_, data.data(), n, SEND_FLAGS);
      if (sent < 0) {
        int err = last_error();
        if (again(err)) continue;
        throw ApiError(0, "the connection broke while sending: " + error_text(err));
      }
      data.remove_prefix(static_cast<std::size_t>(sent));
    }
  }

  std::size_t read(char* out, std::size_t n, double seconds) override {
    for (;;) {
      if (!wait(fd_, false, seconds)) throw Timeout(0, "no answer within " + seconds_text(seconds) + " s");
      auto got = ::recv(fd_, out, static_cast<int>(std::min<std::size_t>(n, INT_MAX)), 0);
      if (got >= 0) return static_cast<std::size_t>(got);
      int err = last_error();
      if (again(err)) continue;
      throw ApiError(0, "the connection broke: " + error_text(err));
    }
  }

  void close() {
    if (fd_ != NONE) close_native(fd_);
    fd_ = NONE;
  }

 private:
#ifdef _WIN32
  using Native = SOCKET;
  static constexpr Native NONE = INVALID_SOCKET;
  static constexpr int SEND_FLAGS = 0;
  static void startup() {
    static const bool started = [] {
      WSADATA data;
      return ::WSAStartup(MAKEWORD(2, 2), &data) == 0;
    }();
    if (!started) throw ApiError(0, "Winsock could not be started");
  }
  static int last_error() { return ::WSAGetLastError(); }
  static bool in_progress(int err) { return err == WSAEWOULDBLOCK || err == WSAEINPROGRESS; }
  static bool again(int err) { return err == WSAEWOULDBLOCK || err == WSAEINTR; }
  static void nonblocking(Native s) {
    u_long on = 1;
    ::ioctlsocket(s, FIONBIO, &on);
  }
  static void close_native(Native s) { ::closesocket(s); }
  static std::string resolve_text(int rc) { return std::system_category().message(rc); }
  // True when the socket is ready (or has failed: the next call says how); false when the time ran out.
  static bool wait(Native s, bool writing, double seconds) {
    fd_set ready, failed;
    FD_ZERO(&ready);
    FD_ZERO(&failed);
    FD_SET(s, &ready);
    FD_SET(s, &failed);
    timeval tv;
    double whole = std::floor(std::max(0.0, seconds));
    tv.tv_sec = static_cast<long>(std::min(whole, 2e9));
    tv.tv_usec = static_cast<long>((std::max(0.0, seconds) - whole) * 1e6);
    int r = ::select(0, writing ? nullptr : &ready, writing ? &ready : nullptr, &failed, &tv);
    if (r < 0) throw ApiError(0, "waiting on the connection failed: " + error_text(last_error()));
    return r > 0;
  }
#else
  using Native = int;
  static constexpr Native NONE = -1;
#ifdef MSG_NOSIGNAL
  static constexpr int SEND_FLAGS = MSG_NOSIGNAL;  // a closed peer is an error, not SIGPIPE
#else
  static constexpr int SEND_FLAGS = 0;
#endif
  static void startup() {}
  static int last_error() { return errno; }
  static bool in_progress(int err) { return err == EINPROGRESS || err == EINTR; }
  static bool again(int err) { return err == EAGAIN || err == EWOULDBLOCK || err == EINTR; }
  static void nonblocking(Native s) { ::fcntl(s, F_SETFL, ::fcntl(s, F_GETFL, 0) | O_NONBLOCK); }
  static void close_native(Native s) { ::close(s); }
  static std::string resolve_text(int rc) { return ::gai_strerror(rc); }
  static bool wait(Native s, bool writing, double seconds) {
    double ms = std::ceil(std::max(0.0, seconds) * 1000.0);
    pollfd p{s, static_cast<short>(writing ? POLLOUT : POLLIN), 0};
    for (;;) {
      int r = ::poll(&p, 1, ms > INT_MAX ? INT_MAX : static_cast<int>(ms));
      if (r < 0 && errno == EINTR) continue;
      if (r < 0) throw ApiError(0, "waiting on the connection failed: " + error_text(errno));
      return r > 0;
    }
  }
#endif
  static std::string error_text(int code) { return std::system_category().message(code); }
  static void options(Native s) {
    int on = 1;
    ::setsockopt(s, IPPROTO_TCP, TCP_NODELAY, reinterpret_cast<const char*>(&on), sizeof on);
#ifdef SO_NOSIGPIPE
    ::setsockopt(s, SOL_SOCKET, SO_NOSIGPIPE, reinterpret_cast<const char*>(&on), sizeof on);
#endif
  }
  Native fd_ = NONE;
};

// An answer: its status line and headers, then its body, decoded from its
// framing as it arrives.
class Message {
 public:
  explicit Message(Input& in) : in_(in) {}

  int status = 0;
  std::string reason;
  std::vector<std::pair<std::string, std::string>> headers;

  // A header's value by its name, in any case; empty when it is absent.
  std::string header(std::string_view name) const {
    std::string want = lower(std::string(name));
    for (const auto& h : headers)
      if (lower(h.first) == want) return h.second;
    return {};
  }

  // The status line and the headers (after any interim 1xx answer), and the framing of the body.
  void read_head(double seconds) {
    for (;;) {
      std::string first = line(seconds);
      if (first.compare(0, 5, "HTTP/") != 0) throw ApiError(0, "the answer is not HTTP: " + first.substr(0, 80));
      std::size_t sp = first.find(' ');
      if (sp == std::string::npos || sp + 4 > first.size() || !std::isdigit(static_cast<unsigned char>(first[sp + 1])) ||
          !std::isdigit(static_cast<unsigned char>(first[sp + 2])) || !std::isdigit(static_cast<unsigned char>(first[sp + 3])))
        throw ApiError(0, "an HTTP status line without its status: " + first.substr(0, 80));
      status = std::stoi(first.substr(sp + 1, 3));
      reason = sp + 5 <= first.size() ? first.substr(sp + 5) : std::string();
      headers.clear();
      for (;;) {
        std::string h = line(seconds);
        if (h.empty()) break;
        if (headers.size() > 200) throw ApiError(0, "an answer with more than 200 headers");
        if ((h[0] == ' ' || h[0] == '\t') && !headers.empty()) {  // an obsolete folded line continues the last
          headers.back().second += ' ' + trim(h);
          continue;
        }
        std::size_t colon = h.find(':');
        if (colon != std::string::npos) headers.emplace_back(trim(std::string_view(h).substr(0, colon)), trim(std::string_view(h).substr(colon + 1)));
      }
      if (status >= 200 || status < 100) break;  // 1xx is interim: the final answer follows
    }
    std::string coding = lower(header("Transfer-Encoding"));
    std::string length = header("Content-Length");
    if (status == 204 || status == 304) {
      framing_ = Framing::none;
    } else if (coding.find("chunked") != std::string::npos) {
      framing_ = Framing::chunked;
    } else if (!length.empty()) {
      if (!std::all_of(length.begin(), length.end(), [](unsigned char c) { return std::isdigit(c); }) || length.size() > 15)
        throw ApiError(0, "a Content-Length that is not a number: " + length);
      framing_ = Framing::length;
      remaining_ = std::stoll(length);
    } else {
      framing_ = Framing::close;
    }
  }

  // Appends the next bytes of the body to out, waiting at most `seconds` for
  // them; false, with nothing appended, once the body has ended.
  bool read_body(std::string& out, double seconds) {
    for (;;) {
      switch (framing_) {
        case Framing::none:
          return false;
        case Framing::close:
          if (take(out, SIZE_MAX, seconds)) return true;
          framing_ = Framing::none;
          return false;
        case Framing::length:
          if (remaining_ == 0) return false;
          if (!take(out, static_cast<std::size_t>(remaining_), seconds))
            throw ApiError(0, "the connection closed " + std::to_string(remaining_) + " bytes before the end of the answer");
          remaining_ -= static_cast<long long>(taken_);
          return true;
        case Framing::chunked:
          if (remaining_ > 0) {
            if (!take(out, static_cast<std::size_t>(remaining_), seconds))
              throw ApiError(0, "the connection closed inside a chunk of the answer");
            remaining_ -= static_cast<long long>(taken_);
            if (remaining_ == 0) after_chunk_ = true;
            return true;
          }
          if (after_chunk_) {  // the CR LF that ends a chunk's data
            if (!line(seconds).empty()) throw ApiError(0, "a chunk of the answer longer than it said");
            after_chunk_ = false;
          }
          remaining_ = chunk_size(line(seconds));
          if (remaining_ == 0) {
            while (!line(seconds).empty()) {}  // the trailer's fields, unused
            framing_ = Framing::none;
            return false;
          }
          continue;
      }
    }
  }

  // The whole body.
  std::string read_all(double seconds) {
    std::string out;
    if (framing_ == Framing::length) out.reserve(static_cast<std::size_t>(remaining_));
    while (read_body(out, seconds)) {}
    return out;
  }

 private:
  enum class Framing { none, length, chunked, close };
  static constexpr std::size_t MAX_LINE = 1 << 16;
  Input& in_;
  std::string buf_;
  std::size_t pos_ = 0, taken_ = 0;
  Framing framing_ = Framing::none;
  long long remaining_ = 0;
  bool after_chunk_ = false;

  // At most `most` bytes of the body into out, from what is buffered or else
  // straight from the input; false at the end of the input.
  bool take(std::string& out, std::size_t most, double seconds) {
    taken_ = 0;
    if (pos_ < buf_.size()) {
      taken_ = std::min(most, buf_.size() - pos_);
      out.append(buf_, pos_, taken_);
      pos_ += taken_;
      return true;
    }
    std::size_t want = std::min<std::size_t>(most, 1 << 16), old = out.size();
    out.resize(old + want);
    std::size_t got = in_.read(&out[old], want, seconds);
    out.resize(old + got);
    taken_ = got;
    return got > 0;
  }

  bool fill(double seconds) {
    if (pos_ == buf_.size()) {
      buf_.clear();
      pos_ = 0;
    } else if (pos_ > MAX_LINE) {
      buf_.erase(0, pos_);
      pos_ = 0;
    }
    std::size_t old = buf_.size();
    buf_.resize(old + 4096);
    std::size_t got = in_.read(&buf_[old], 4096, seconds);
    buf_.resize(old + got);
    return got > 0;
  }

  // A line of the head or of the chunks' framing, without its CR LF.
  std::string line(double seconds) {
    for (;;) {
      std::size_t nl = buf_.find('\n', pos_);
      if (nl != std::string::npos) {
        std::string out = buf_.substr(pos_, nl - pos_);
        pos_ = nl + 1;
        if (!out.empty() && out.back() == '\r') out.pop_back();
        return out;
      }
      if (buf_.size() - pos_ > MAX_LINE) throw ApiError(0, "a line of the answer longer than 64 KiB");
      if (!fill(seconds)) throw ApiError(0, "the connection closed in the middle of the answer");
    }
  }

  static long long chunk_size(const std::string& text) {
    long long n = 0;
    std::size_t k = 0;
    for (; k < text.size() && std::isxdigit(static_cast<unsigned char>(text[k])); k++) {
      if (n > (1LL << 40)) throw ApiError(0, "a chunk of the answer too large to be true");
      char c = static_cast<char>(std::tolower(static_cast<unsigned char>(text[k])));
      n = n * 16 + (c <= '9' ? c - '0' : c - 'a' + 10);
    }
    if (k == 0 || (k < text.size() && text[k] != ';' && text[k] != ' ' && text[k] != '\t'))
      throw ApiError(0, "a chunk's size that is not hexadecimal: " + text.substr(0, 40));
    return n;  // what follows ';' are chunk extensions, unused
  }
};

// ---- Server-Sent Events ---------------------------------------------------------------------------

// One event of a stream: its type (event:), its data (data: lines joined by
// LF), and the last event id seen when it was dispatched.
struct ServerEvent {
  std::string type, data, id;
};

// Events from the bytes of a text/event-stream, as the WHATWG HTML Standard
// §9.2.6 interprets them: a line beginning with ':' is a comment (the service's
// heartbeats); an empty line dispatches the event gathered, when it has data;
// the id persists from event to event, as the Last-Event-ID to reconnect with;
// retry: sets the delay before reconnecting. Lines end in LF or CR LF; a lone
// CR, which the Standard also allows and the service never sends, is not taken
// as an end.
class EventStream {
 public:
  std::string last_id;  // the Last-Event-ID to reconnect with
  long retry_ms = -1;   // the delay the service asked for, or -1

  // Feeds bytes. For each event dispatched calls on_event(const ServerEvent&),
  // for each comment on_comment(std::string_view); either returning false stops,
  // and then feed returns false.
  template <class OnEvent, class OnComment>
  bool feed(std::string_view bytes, OnEvent&& on_event, OnComment&& on_comment) {
    pending_.append(bytes.data(), bytes.size());
    std::size_t start = 0;
    bool go = true;
    for (;;) {
      std::size_t nl = pending_.find('\n', start);
      if (nl == std::string::npos) break;
      std::string_view line(pending_.data() + start, nl - start);
      start = nl + 1;
      if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
      if (!take(line, on_event, on_comment)) {
        go = false;
        break;
      }
    }
    pending_.erase(0, start);
    if (pending_.size() + data_.size() > (16u << 20)) throw ApiError(0, "an event of the stream longer than 16 MiB");
    return go;
  }

  // Forgets a partial event, as after a reconnection; the last id is kept.
  void reset() {
    pending_.clear();
    data_.clear();
    type_.clear();
    has_data_ = false;
  }

 private:
  std::string pending_, data_, type_;
  bool has_data_ = false;

  template <class OnEvent, class OnComment>
  bool take(std::string_view line, OnEvent& on_event, OnComment& on_comment) {
    if (line.empty()) {
      if (!has_data_) {
        type_.clear();
        return true;
      }
      ServerEvent e{type_.empty() ? std::string("message") : type_, data_, last_id};
      if (!e.data.empty()) e.data.pop_back();  // the LF after the last data line
      data_.clear();
      type_.clear();
      has_data_ = false;
      return on_event(static_cast<const ServerEvent&>(e));
    }
    if (line[0] == ':') return on_comment(line.substr(1));
    std::size_t colon = line.find(':');
    std::string_view field = line.substr(0, colon);
    std::string_view value = colon == std::string_view::npos ? std::string_view() : line.substr(colon + 1);
    if (!value.empty() && value[0] == ' ') value.remove_prefix(1);
    if (field == "data") {
      data_.append(value.data(), value.size());
      data_ += '\n';
      has_data_ = true;
    } else if (field == "event") {
      type_ = std::string(value);
    } else if (field == "id") {
      if (value.find('\0') == std::string_view::npos) last_id = std::string(value);
    } else if (field == "retry") {
      if (!value.empty() && value.size() < 10 &&
          std::all_of(value.begin(), value.end(), [](unsigned char c) { return std::isdigit(c); }))
        retry_ms = std::stol(std::string(value));
    }
    return true;
  }
};

}  // namespace http
}  // namespace science
