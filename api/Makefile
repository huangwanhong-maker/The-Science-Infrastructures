# The C and C++ clients of the API (README.md §3.5): the C interface as
# libscience.a and libscience.so (c/science.cpp over the header cpp/science.hpp),
# the examples of examples/c/ and examples/cpp/, and the check programs of tests/.
#
#   make -C api                         the libraries and the examples, in api/build/
#   make -C api check                   the check programs (tests/run-tests.mjs builds and runs them)
#   make -C api CXX=clang++ CC=clang    another compiler; CXXSTD=c++20 for C++20
#   make -C api BUILD=/tmp/science      another place for what is built
#   make -C api clean
#
# The header needs a C++17 compiler, its standard library and the system's
# sockets, nothing else; a C program links the library with the C++ runtime, as
# the rules for the C examples here do. GNU make, on Linux and macOS; on Windows
# with MinGW, the sockets come from -lws2_32 and the shared library is a DLL.
# Paths here are relative to api/, since make cannot hold a path with a space
# in it, and the repository's own has one.

BUILD    ?= build
CXXSTD   ?= c++17
CSTD     ?= c99
CXXFLAGS ?= -O2 -Wall -Wextra -pedantic
CFLAGS   ?= -O2 -Wall -Wextra -pedantic
THREADS  ?= -pthread

ifeq ($(OS),Windows_NT)
  NET    := -lws2_32
  SHARED := $(BUILD)/libscience.dll
  EXE    := .exe
else
  NET    :=
  SHARED := $(BUILD)/libscience.so
  EXE    :=
endif

HEADERS      := cpp/science.hpp cpp/science/json.hpp cpp/science/http.hpp
CPP_EXAMPLES := $(patsubst examples/cpp/%.cpp,$(BUILD)/cpp-%$(EXE),$(wildcard examples/cpp/*.cpp))
C_EXAMPLES   := $(patsubst examples/c/%.c,$(BUILD)/c-%$(EXE),$(wildcard examples/c/*.c))
CHECKS       := $(BUILD)/check_cpp$(EXE) $(BUILD)/check_cpp20$(EXE) $(BUILD)/check_c$(EXE)

.PHONY: all libraries examples check clean
all: libraries examples
libraries: $(BUILD)/libscience.a $(SHARED)
examples: $(CPP_EXAMPLES) $(C_EXAMPLES)
check: $(CHECKS)

$(BUILD) $(BUILD)/o:
	mkdir -p $@

# The C interface: one object, archived and shared.
$(BUILD)/o/science.o: c/science.cpp c/science.h $(HEADERS) | $(BUILD)/o
	$(CXX) -std=$(CXXSTD) $(CXXFLAGS) -fPIC -c $< -o $@

$(BUILD)/libscience.a: $(BUILD)/o/science.o | $(BUILD)
	$(AR) rcs $@ $^

$(SHARED): $(BUILD)/o/science.o | $(BUILD)
	$(CXX) -shared $^ -o $@ $(NET)

# The C++ examples include the header alone; the C ones link the library.
$(BUILD)/cpp-%$(EXE): examples/cpp/%.cpp $(HEADERS) | $(BUILD)
	$(CXX) -std=$(CXXSTD) $(CXXFLAGS) -I cpp $< -o $@ $(THREADS) $(NET)

$(BUILD)/o/c-%.o: examples/c/%.c c/science.h | $(BUILD)/o
	$(CC) -std=$(CSTD) $(CFLAGS) -I c -c $< -o $@

$(BUILD)/c-%$(EXE): $(BUILD)/o/c-%.o $(BUILD)/libscience.a
	$(CXX) $^ -o $@ $(THREADS) -lm $(NET)

# The checks: the C++ one as C++17 and, to hold the header to it, as C++20.
$(BUILD)/check_cpp$(EXE): tests/check_cpp.cpp $(HEADERS) | $(BUILD)
	$(CXX) -std=$(CXXSTD) $(CXXFLAGS) -I cpp $< -o $@ $(THREADS) $(NET)

$(BUILD)/check_cpp20$(EXE): tests/check_cpp.cpp $(HEADERS) | $(BUILD)
	$(CXX) -std=c++20 $(CXXFLAGS) -I cpp $< -o $@ $(THREADS) $(NET)

$(BUILD)/o/check_c.o: tests/check_c.c c/science.h | $(BUILD)/o
	$(CC) -std=$(CSTD) $(CFLAGS) -I c -c $< -o $@

$(BUILD)/check_c$(EXE): $(BUILD)/o/check_c.o $(BUILD)/libscience.a
	$(CXX) $^ -o $@ $(THREADS) -lm $(NET)

clean:
	rm -rf $(BUILD)
