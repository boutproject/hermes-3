//
// Test the speed of passing state between components 
//
// Result: (Intel(R) Core(TM) i5-7200U CPU @ 2.50GHz)
//
// Global state: 4.58632e-05s
// Options state: 6.65891e-05s
// Options state with resets: 6.48178e-05s

// -> Overhead of around 20 us, re-using state doesn't seem to improve time

#include <bout.hxx>
#include <field3d.hxx>
#include <options.hxx>

#include <chrono>

////////////////////////////////////////////////
// 1. Global state
//
// This is a baseline case with minimal/no overhead

namespace globalstate {
  
  Field3D a, b, c, d, e, f, g, h, i, j, k, l;
  BoutReal x, y, z;
  
  void component1() {
    a = 1.0;
    b = 2.0;
    c = 3.0;
    d = 4.0;

    x = -5.2;
  }

  void component2() {
    e = 2 * a + b;
    f = 5.0 * x;
    g = d - c;
    h = 6.0;

    y = 42;
  }

  void component3() {
    i = 7.0;
    j = d + f + g;
    l = h - a * y;
    k = a + e + i;
    z = 32;
  }

  void component4() {
    l = b + f + j * y;
  }

  /// Run components in order
  void run() {
    component1();
    component2();
    component3();
    component4();
  }
  
} // globalstate


////////////////////////////////////////////////
// Options class, with some nesting

#include "../include/component.hxx"

namespace optionstate {
  Options state;
  
  void component1(Options &state) {
    state["a"] = Field3D(1.0);
    state["b"] = Field3D(2.0);
    state["c"] = Field3D(3.0);
    state["d"] = Field3D(4.0);

    state["x"] = -5.2;
  }

  void component2(Options &state) {
    state["e"] = 2 * get<Field3D>(state["a"]) + get<Field3D>(state["b"]);
    state["f"] = 5.0 * get<BoutReal>(state["x"]);
    state["g"] = get<Field3D>(state["d"]) - get<Field3D>(state["c"]);
    state["h"] = 6.0;

    state["y"] = 42.0;
  }

  void component3(Options &state) {
    state["i"] = 7.0;
    state["j"] = get<Field3D>(state["d"]) + get<Field3D>(state["f"]) + get<Field3D>(state["g"]);
    state["l"] = get<Field3D>(state["h"]) - get<Field3D>(state["a"]) * get<BoutReal>(state["y"]);
    state["k"] = get<Field3D>(state["a"]) + get<Field3D>(state["e"]) + get<Field3D>(state["i"]);
    state["z"] = 32.0;
  }

  void component4(Options &state) {
    state["l"] = get<Field3D>(state["b"]) + get<Field3D>(state["f"]) + get<Field3D>(state["j"]) * get<BoutReal>(state["y"]);
  }

  /// Run components in order
  void run() {
    component1(state);
    component2(state);
    component3(state);
    component4(state);
  }
} // optionstate

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  int N = 10000;
  
  globalstate::run();

  {
    auto start = std::chrono::steady_clock::now();
    for(int i = 0; i < N; i++) {
      globalstate::run();
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    output << "Global state: " << elapsed_seconds.count() / N << "s\n";
  }

  optionstate::run();
  {
    auto start = std::chrono::steady_clock::now();
    for(int i = 0; i < N; i++) {
      optionstate::run();
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    output << "Options state: " << elapsed_seconds.count() / N << "s\n";
  }

  // Try resetting the state between runs
  {
    auto start = std::chrono::steady_clock::now();
    for(int i = 0; i < N; i++) {
      optionstate::run();
      optionstate::state = Options{}; // reset state
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    output << "Options state with resets: " << elapsed_seconds.count() / N << "s\n";
  }
  
  BoutFinalise();
}
