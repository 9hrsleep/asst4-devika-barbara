/**
 * Parallel VLSI Wire Routing via MPI
 * Name 1(andrew_id 1), Name 2(andrew_id 2)
 */

#ifndef __WIREOPT_H__
#define __WIREOPT_H__

#include <cstdint>
#include <vector>
#include <array>

#define MAX_PTS_PER_WIRE 5
#define COST_REPORT_DEPTH 10

/** README(student):
 We provide a way to validate consistency between wire layout
 and occupancy. Within the program, a wire checker is provided 
 by:
 wr_checker Checker(wires, occupancy);
 and its validate() method can be called to validate the consistency.
 
 The struct below is the standard format for wires used by the wire checker.
 It contains a buffer that holds up to MAX_PTS_PER_WIRE points, and a num_pts
 field that specifies the number of points. Regardless of what representation you use for your wires, you should
 implement the Wire::to_validate_format method to convert your Wire
 to a validate_wire_t if you wish to use the checker.
*/

/* validate_wire_t is a format to represent wires by key points

For example a wire with 3 bend has 5 key poinst:
start, bend1, bend2, bend3, end.
Notice this implies two consecutive keypoints share at least same x or same y

num_pts is number of keypoints
p is an array store the xy coordinates of each keypoints

Requires (see validate.cpp validate_wire_t::cleanup() for how is this
checked）： (1) two consecutive keypoints share at least same x or same y, (2) 2
<= num_pts <= MAX_PTS_PER_WIRE.
*/
struct validate_wire_t {
  uint8_t num_pts;
  struct {
    uint16_t x;
    uint16_t y;
  } p[MAX_PTS_PER_WIRE];
  validate_wire_t &cleanup(void);
  void print_wire(void) const;
};

/* TODO (student): Define the data structure for wire here.
Feel free to redefine this, notice you also need to change how bend is read in
inside main of wireroute.cpp.

A better readable way of define this is:
int start_x, start_y, end_x, end_y, bend1_x, bend1_y, bend2_x, bend2_y, bend3_x,
bend3_y; but this might not be the most efficient way to define the solution
space for a wire with <= 3 bends.
*/
struct Point {
  int x, y;
};

struct Wire {
  int start_x, start_y;
  int end_x, end_y;
  std::array<Point, 3> bends{};
  int num_bends = 0;
  validate_wire_t to_validate_format(void) const;
};

// Definition of the wire checker
struct wr_checker {
  std::vector<Wire> wires;
  std::vector<std::vector<int>> occupancies;
  const int nwires;
  const int dim_x;
  const int dim_y;
  wr_checker(std::vector<Wire> &wires,
             std::vector<std::vector<int>> &occupancies)
      : wires(wires), occupancies(occupancies), nwires(wires.size()),
        dim_x(occupancies[0].size()), dim_y(occupancies.size()) {}
  void validate() const;
};

const char *get_option_string(const char *option_name,
                              const char *default_value);
int get_option_int(const char *option_name, int default_value);
float get_option_float(const char *option_name, float default_value);

#endif
