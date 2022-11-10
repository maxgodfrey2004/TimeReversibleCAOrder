#include <cassert>
#include <cstdlib>

#include <bitset>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

using std::bitset;
using std::get;
using std::size_t;
using std::string;
using std::vector;

// Represents the value of a cell (stored in the second value) and its
// immediate neighborhood (stored in the leftmost and rightmost values
// respectively).
//
using neighborhood = std::tuple<int, int, int>;

// Returns the `n`-th bit of an integer `x`.
// 
inline bool get_bit(int x, int n) {
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

// Applies a Cellular Automaton rule for a single cell given the state of that
// cell and its immediate neighborhood.
//
bool apply_forward_rule(int rule_num, neighborhood state_cur) {
  int idx = get<0>(state_cur) * 4 + get<1>(state_cur) * 2 + get<2>(state_cur);
  return get_bit(rule_num, idx);
}

// Applies a Time-Reversible Cellular Automaton for a single cell given the
// current state of that cell and its immediate neighborhood, and the cell's
// past state.
//
bool apply_reverse_rule(int rule_num, neighborhood state_cur, bool state_past) {
  int idx = get<0>(state_cur) * 4 + get<1>(state_cur) * 2 + get<2>(state_cur);
  return get_bit(rule_num, idx) ^ state_past;
}

// Runs a Cellular Automaton for `N` cells. Returns the state of all cells at
// each iteration of the rule application, including the starting state.
//
template <size_t N>
vector<bitset<N>> ca_forward_rule(int rule_num, int steps, bitset<N> m0) {
  vector<bitset<N>> states(steps);
  states[0] = m0;
  for (int step = 0; step + 1 < steps; ++step) {
    for (size_t i = 0; i < N; ++i) {
      states[step + 1][i] = apply_forward_rule(
        rule_num,
        neighborhood(
          states[step][index_left(i, N)], 
          states[step][i],
          states[step][index_right(i, N)]
        )
      );
    }
  }
  return states;
}

// Runs a Time-Reversible Cellular Automaton for `N` cells. Returns the state
// of all cells at each iteration of the rule application, including the
// starting states.
//
// Note: this implementation uses wrapping, such that the first and last cells
//       are considered adjacent. 
//
template <size_t N>
vector<bitset<N>> ca_reverse_rule(
  int rule_num,
  int steps,
  bitset<N> m0,  // The initial past state of the cells.
  bitset<N> m1   // The initial current state of the cells.
) {
  vector<bitset<N>> states(steps);
  states[0] = m0;
  states[1] = m1;
  for (int step = 1; step + 1 < steps; ++step) {
    for (size_t i = 0; i < N; ++i) {
      states[step + 1][i] = apply_reverse_rule(
        rule_num,
        neighborhood(
          states[step][index_left(i, N)],
          states[step][i],
          states[step][index_right(i, N)]
        ),
        states[step - 1][i]
      );
    }
  }
  return states;
}

// Writes a bitset to standard output, replacing the ones with solid characters
// and outputting spaces wherever there is a zero.
//
template <size_t N>
void write_bitset(bitset<N> b) {
  for (size_t i = 0; i < N; ++i) {
    std::cout << (b[i] ? '#' : ' ');
  }
  std::cout << std::endl;
}

// Unit tests for routines provided in this file.
//
void run_tests() {
  std::cout << "Running tests..." << std::endl;
  assert(index_left(0, 3) == 2);
  assert(index_left(2, 3) == 1);
  assert(index_right(0, 3) == 1);
  assert(index_right(2, 3) == 0);
  std::cout << ">  Index tests passed." << std::endl;
  assert(apply_reverse_rule(110,  neighborhood(1, 0, 1), 1) == 0);
  assert(apply_reverse_rule(30, neighborhood(1, 0, 0), 0) == 1);
  std::cout << ">  Reverse rule tests passed." << std::endl;
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
  
  // run_tests();
  
  // bitset<18> m0;
  // bitset<18> m1("111110011000111111");
  // std::cout << "\nFORWARD:" << std::endl;
  // for (auto line : ca_forward_rule(37, 18, m1)) {
  //   write_bitset(line);
  // }
  // std::cout << "\nREVERSE:" << std::endl;
  // for (auto line : ca_reverse_rule(37, 18, m0, m1)) {
  //   write_bitset(line);
  // }
}