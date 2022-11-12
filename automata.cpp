#include <cassert>
#include <cstdlib>

#include <bitset>
#include <fstream>
#include <iomanip>
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

// The number of input states used to analyze a specific Cellular Automaton.
const int NUM_SAMPLES = 100;

// Represents the result of finding a cycle for a Cellular Automaton.
//
template <size_t N>
struct forward_cycle_result {
  int steps;
  int total_hamming_distance;
  bitset<N> state_cur;

  forward_cycle_result(int s, int thd, bitset<N> sc)
    : steps(s), total_hamming_distance(thd), state_cur(sc) {}

  string to_string() {
    return "{STEPS=" + std::to_string(steps) + ", THD="
      + std::to_string(total_hamming_distance) + ", SC="
      + state_cur.to_string() + "}";
  }
};

// Represents the result of finding a cycle for a Time-Reversible Cellular
// Automaton.
//
template <size_t N>
struct reverse_cycle_result {
  int steps;
  int total_hamming_distance;
  bitset<N> state_prev;
  bitset<N> state_cur;

  reverse_cycle_result(int s, int thd, bitset<N> sp, bitset<N> sc)
    : steps(s), total_hamming_distance(thd), state_prev(sp), state_cur(sc) {}
  
  string to_string() {
    return "{STEPS=" + std::to_string(steps) + ", THD="
      + std::to_string(total_hamming_distance) + ", SP="
      + state_prev.to_string() + ", SC=" + state_cur.to_string() + "}";
  }
};

// Represents the result of analyzing a cycle.
//
struct analysis_result {
  double avg_cycle_length;
  double avg_hamming_distance;

  analysis_result() : avg_cycle_length(0), avg_hamming_distance(0) {}
  analysis_result(double cl, double hd)
    : avg_cycle_length(cl), avg_hamming_distance(hd) {}
  
  string to_string() {
    return "{AVG_CYC_LEN=" + std::to_string(avg_cycle_length) + ", AVG_HD="
      + std::to_string(avg_hamming_distance) + "}";
  }
};

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

// Returns the Hamming Distance between two bitsets of the same length.
//
template <size_t N>
int hamming_distance(bitset<N> a, bitset<N> b) {
  int dist = 0;
  for (size_t i = 0; i < N; ++i) {
    dist += (a[i] ^ b[i]);
  }
  return dist;
}

// Generates a random bitset of size N.
//
template <size_t N>
bitset<N> random_bitset() {
  bitset<N> res;
  for (size_t i = 0; i < N; ++i) {
    res[i] = (rand() % 2);
  }
  return res;
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

// Returns the number of steps, total hamming distance, and repeated state
// found when searching a Cellular Automaton for a cycle given some input
// state.
//
template <size_t N>
forward_cycle_result<N> find_forward_cycle(int rule_num, bitset<N> state_cur) {
  bitset<1 << N> seen_states;
  bitset<N> state_next;
  int steps = 0;
  int total_hamming_distance = 0;
  while (!seen_states[state_cur.to_ulong()]) {
    seen_states[state_cur.to_ulong()] = 1;
    ++steps;
    for (size_t i = 0; i < N; ++i) {
      state_next[i] = apply_forward_rule(
        rule_num,
        neighborhood(
          state_cur[index_left(i, N)],
          state_cur[i],
          state_cur[index_right(i, N)]
        )
      );
    }
    total_hamming_distance += hamming_distance(state_cur, state_next);
    state_cur = state_next;
  }
  return forward_cycle_result<N>(steps, total_hamming_distance, state_cur);
}

// Returns the number of steps, total hamming distance, and repeated state
// found when searching a Time-Reversible Cellular Automaton for a cycle given
// some input state.
//
template <size_t N>
reverse_cycle_result<N> find_reverse_cycle(
  int rule_num,
  bitset<N> state_prev,
  bitset<N> state_cur
) {
  // Each state comprises the values of the cells in the previous step and the
  // current values of the cells. This means that there are 2 * (2^N) possible
  // states.
  // We can map each state to a unique integer since the values of all cells
  // can be expressed as a bitset of length N. If A and B represent the
  // integers represented by the current and past bitsets respectively, then we
  // can assign this state the integer ((A << N) + B).
  bitset<1 << (N * 2)> seen_states;
  bitset<N> state_next;
  int steps = 0;
  int total_hamming_distance = 0;
  while (!seen_states[(state_cur.to_ulong() << N) + state_prev.to_ulong()]) {
    seen_states[(state_cur.to_ulong() << N) + state_prev.to_ulong()] = 1;
    ++steps;
    for (size_t i = 0; i < N; ++i) {
      state_next[i] = apply_reverse_rule(
        rule_num,
        neighborhood(
          state_cur[index_left(i, N)],
          state_cur[i],
          state_cur[index_right(i, N)]
        ),
        state_prev[i]
      );
    }
    total_hamming_distance += hamming_distance(state_cur, state_next);
    state_prev = state_cur;
    state_cur = state_next;
  }
  return reverse_cycle_result<N>(
    steps,
    total_hamming_distance,
    state_prev,
    state_cur
  );
}

// Determines the cycle length and average hamming distance for the repeated
// state tended to by a given Cellular Automaton and input state.
//
template <size_t N>
struct analysis_result analyze_forward_cycle(
  int rule_num,
  bitset<N> state_cur
) {
  // Determine a state within the cycle.
  struct forward_cycle_result<N> rep = find_forward_cycle(rule_num, state_cur);
  // Determine information for the cycle.
  struct forward_cycle_result<N> cyc = find_forward_cycle(rule_num, rep.state_cur);
  double avg_cycle_length = double(cyc.steps) / double(1);
  double avg_hamming_distance = double(cyc.total_hamming_distance) / double(cyc.steps);
  return analysis_result(avg_cycle_length, avg_hamming_distance);
}

// Determines the cycle length and average hamming distance for the repeated
// state tended to by a given Time-Reversible Cellular Automaton and input
// state.
//
template <size_t N>
struct analysis_result analyze_reverse_cycle(
  int rule_num,
  bitset<N> state_prev,
  bitset<N> state_cur
) {
  // Determine a state within the cycle.
  struct reverse_cycle_result<N> rep = find_reverse_cycle(rule_num, state_prev, state_cur);
  // Determine information for the cycle.
  struct reverse_cycle_result<N> cyc = find_reverse_cycle(rule_num, rep.state_prev, rep.state_cur);
  double avg_cycle_length = double(cyc.steps) / double(1);
  double avg_hamming_distance = double(cyc.total_hamming_distance) / double(cyc.steps);
  return analysis_result(avg_cycle_length, avg_hamming_distance);
}

// Determines the average cycle lengths and average hamming distances for a
// given input size `N` and rule number and its Time-Reversible cognate.
//
template <size_t N>
void analyze_rule(int rule_num) {
  // Compute the forward result.
  struct analysis_result result_forward;
  for (int i = 0; i < NUM_SAMPLES; ++i) {
    bitset<N> m1 = random_bitset<N>();
    struct analysis_result res = analyze_forward_cycle(rule_num, m1);
    result_forward.avg_cycle_length += res.avg_cycle_length;
    result_forward.avg_hamming_distance += res.avg_hamming_distance;
  }
  result_forward.avg_cycle_length /= double(NUM_SAMPLES);
  result_forward.avg_hamming_distance /= double(NUM_SAMPLES);
  // Compute the reverse result.
  struct analysis_result result_reverse;
  for (int i = 0; i < NUM_SAMPLES; ++i) {
    bitset<N> m0 = random_bitset<N>();
    bitset<N> m1 = random_bitset<N>();
    struct analysis_result res = analyze_reverse_cycle(rule_num, m0, m1);
    result_reverse.avg_cycle_length += res.avg_cycle_length;
    result_reverse.avg_hamming_distance += res.avg_hamming_distance;
  }
  result_reverse.avg_cycle_length /= double(NUM_SAMPLES);
  result_reverse.avg_hamming_distance /= double(NUM_SAMPLES);
  // Write the results to the console. TR stands for "Time-Reversed" and DC
  // stands for "Deterministic Cognate".
  std::cout << "TR: " << result_reverse.to_string() << "\t DC: "
            << result_forward.to_string() << std::endl;
}

// Writes an analysis of all 256 1-Dimensional Automaton (Time-Reversible and
// their deterministic cognates) to a specified output file in CSV format.
//
template <size_t N>
void full_analysis(string outfile_name) {
  std::ofstream outfile(outfile_name);
  outfile << std::fixed << std::setprecision(2);
  outfile << "FD ACL,FD AHD,TR ACL,TR AHD" << std::endl;
  for (int rule_num = 0; rule_num < 256; ++rule_num) {
    // Compute the forward result.
    struct analysis_result result_forward;
    for (int i = 0; i < NUM_SAMPLES; ++i) {
      std::cout << "Analyzing Rule " << rule_num << ", FD Sample " << i << "\r"
                << std::flush;
      bitset<N> m1 = random_bitset<N>();
      struct analysis_result res = analyze_forward_cycle(rule_num, m1);
      result_forward.avg_cycle_length += res.avg_cycle_length;
      result_forward.avg_hamming_distance += res.avg_hamming_distance;
    }
    result_forward.avg_cycle_length /= double(NUM_SAMPLES);
    result_forward.avg_hamming_distance /= double(NUM_SAMPLES);
    // Compute the reverse result.
    struct analysis_result result_reverse;
    for (int i = 0; i < NUM_SAMPLES; ++i) {
      std::cout << "Analyzing Rule " << rule_num << ", TR Sample " << i << "\r"
                << std::flush;
      bitset<N> m0 = random_bitset<N>();
      bitset<N> m1 = random_bitset<N>();
      struct analysis_result res = analyze_reverse_cycle(rule_num, m0, m1);
      result_reverse.avg_cycle_length += res.avg_cycle_length;
      result_reverse.avg_hamming_distance += res.avg_hamming_distance;
    }
    result_reverse.avg_cycle_length /= double(NUM_SAMPLES);
    result_reverse.avg_hamming_distance /= double(NUM_SAMPLES);
    outfile << result_forward.avg_cycle_length << ","
            << result_forward.avg_hamming_distance << ","
            << result_reverse.avg_cycle_length << ","
            << result_reverse.avg_hamming_distance << std::endl;
  }
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
  srand(100);  // Produce consistent output every time.

  // Testing rule applications and wraparound index calculation.
  //
  // run_tests();
  
  // Confirming that forward and reverse rule computations are done correctly.
  //
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

  // Confirming that cycles are found and analyzed correctly.
  //
  // bitset<12> m0("101101010110");
  // bitset<12> m1("110110111010");
  // std::cout << find_reverse_cycle(37, m0, m1).to_string() << std::endl;
  // std::cout << analyze_reverse_cycle(37, m0, m1).to_string() << std::endl;
  // std::cout << find_forward_cycle(37, m1).to_string() << std::endl;
  // std::cout << analyze_forward_cycle(37, m1).to_string() << std::endl;

  // Analyze every rule.
  full_analysis<12>("results.csv");
}