#include <cassert>
#include <cstdlib>

#include <iostream>
#include <string>
#include <tuple>

using std::get;
using std::size_t;
using std::string;

// Represents the value of a cell (stored in the second value) and its
// immediate neighborhood (stored in the leftmost and rightmost values
// respectively).
//
using neighborhood = std::tuple<int, int, int>;

// Returns the `n`-th bit of an integer `x`.
// 
inline int get_bit(int x, int n) {
  return ((x >> n) & 1);
}

// Returns the index to the left of a given position in an array of size
// `size`.
//
inline int index_left(int idx, int size) {
  return (((idx - 1) % size) + size) % size;
}

// Returns the index to the right of a given position in an array of size
// `size`.
//
inline int index_right(int idx, int size) {
  return (idx + 1) % size;
}

// Applies a Time-Reversible Cellular Automata Rule for a single cell given
// the current state of that cell and its immediate neighborhood, and the
// cell's past state.
//
int apply_reverse_rule(int rule_num, neighborhood state_cur, int state_past) {
  int idx = get<0>(state_cur) * 4 + get<1>(state_cur) * 2 + get<2>(state_cur);
  return get_bit(rule_num, idx) ^ state_past;
}

// Unit tests for routines provided in this file.
//
void run_tests() {
  std::cout << "Running tests..." << std::endl;
  assert(apply_reverse_rule(110,  neighborhood(1, 0, 1), 1) == 0);
  assert(apply_reverse_rule(30, neighborhood(1, 0, 0), 0) == 1);
  std::cout << "All tests passed." << std::endl;
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    std::cerr << "Error: usage is " << argv[0] << "<space width> "
              << "<number of steps>" << std::endl;
  }
  string str_width(argv[1]);
  string str_steps(argv[2]);
  if (str_width.empty() ||
      str_width.find_first_not_of("0123456789") != std::string::npos ||
      str_steps.empty() ||
      str_steps.find_first_not_of("0123456789") != std::string::npos) {
    std::cerr << "Error: command line args <space width> and <number of steps>"
              << " must be integers." << std::endl; 
  }
  int input_width = std::stoi(str_width);
  int num_steps = std::stoi(str_steps);
  run_tests();
}