/**
 * Parallel VLSI Wire Routing via MPI
 * Name 1(andrew_id 1), Name 2(andrew_id 2)
 */

#include "wireroute.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <mpi.h>
#include <unistd.h>

void print_stats(const std::vector<std::vector<int>> &occupancy) {
  int max_occupancy = 0;
  long long total_cost = 0;

  for (const auto &row : occupancy) {
    for (const int count : row) {
      max_occupancy = std::max(max_occupancy, count);
      total_cost += count * count;
    }
  }

  std::cout << "Max occupancy: " << max_occupancy << '\n';
  std::cout << "Total cost: " << total_cost << '\n';
}

/* This function write the output into 2 files
(1) It write occupancy grids into a file
(2) It convert wires from Wire to validate_wire_t by to_validate_format
(2) It write wires into another file
*/
void write_output(
    const std::vector<Wire> &wires, const int num_wires,
    const std::vector<std::vector<int>> &occupancy, const int dim_x,
    const int dim_y,
    std::string wires_output_file_path = "outputs/wire_output.txt",
    std::string occupancy_output_file_path = "outputs/occ_output.txt") {

  std::ofstream out_occupancy(occupancy_output_file_path, std::fstream::out);
  if (!out_occupancy) {
    std::cerr << "Unable to open file: " << occupancy_output_file_path << '\n';
    exit(EXIT_FAILURE);
  }
  out_occupancy << dim_x << ' ' << dim_y << '\n';

  for (const auto &row : occupancy) {
    for (size_t i = 0; i < row.size(); ++i)
      out_occupancy << row[i] << (i == row.size() - 1 ? "" : " ");
    out_occupancy << '\n';
  }
  out_occupancy.close();

  std::ofstream out_wires(wires_output_file_path, std::fstream::out);
  if (!out_wires) {
    std::cerr << "Unable to open file: " << wires_output_file_path << '\n';
    exit(EXIT_FAILURE);
  }

  out_wires << dim_x << ' ' << dim_y << '\n';
  out_wires << num_wires << '\n';

  for (const auto &wire : wires) {
    // NOTICE: we convert to keypoint representation here, using
    // to_validate_format which need to be defined in the bottom of this file
    validate_wire_t keypoints = wire.to_validate_format();
    for (int i = 0; i < keypoints.num_pts; ++i) {
      out_wires << keypoints.p[i].x << ' ' << keypoints.p[i].y;
      if (i < keypoints.num_pts - 1)
        out_wires << ' ';
    }
    out_wires << '\n';
  }

  out_wires.close();
}


/** from a chosen wire, get the point coordinates of starting point (at index 0),
or a bend points (index 1, 2, 3) or end point (index 4)
*/
Point get_point_i(const Wire &wire, int i) {
  if (i == 0){
    return {wire.start_x, wire.start_y};
  } else if (i == wire.num_bends + 1) {
    return {wire.end_x, wire.end_y};
  } else {
    return wire.bends[i-1];
  }
}

// update both copies of one square. across wires needs atomics because
// different threads can update the same square at the same time
inline void update_occupancy_cell(std::vector<std::vector<int>> &occupancy,
                           std::vector<std::vector<int>> &occupancy_T,
                           int x, int y, int delta, bool atomic_updates) {
  occupancy[y][x] += delta;
  occupancy_T[x][y] += delta;
}

/** this updates the occupancy matrix by adding or subtracting 1 in the entire
 wire's path on the matrix. to add a wire, use delta = 1 and to subtract
 a wire, use delta = -1. we also update occupancy_T (the transposed copy,
 occupancy_T[x][y] == occupancy[y][x]) so the two agree once all updates finish.
*/
void update_occupancy_with_wire(std::vector<std::vector<int>> &occupancy,
                                std::vector<std::vector<int>> &occupancy_T,
                                const Wire &wire, int delta,
                                bool atomic_updates = false) {
  Point curr_pt = get_point_i(wire, 0);
  update_occupancy_cell(occupancy, occupancy_T, curr_pt.x, curr_pt.y,
                        delta, atomic_updates);

  // loop through all line segments and add the delta (+1 or -1) on all points you traverse
  const int num_segments = wire.num_bends + 1;
  for (int s = 0; s < num_segments; s++) {
    // get the end point where this current segment will end
    Point next_pt = get_point_i(wire, s + 1);
    // figure out if we are adding to the x coordinate or subtracting
    int step_x = 0;
    if (next_pt.x > curr_pt.x) {
      step_x = 1;
    } else if (next_pt.x < curr_pt.x) {
      step_x = -1;
    }
    // figure out if we are adding to the y coordinate or subtracting
    int step_y = 0;
    if (next_pt.y > curr_pt.y) {
      step_y = 1;
    } else if (next_pt.y < curr_pt.y) {
      step_y = -1;
    }

    // move through from curr_pt to next_pt, updating both grids
    while (curr_pt.x != next_pt.x || curr_pt.y != next_pt.y) {
      curr_pt.x += step_x;
      curr_pt.y += step_y;
      update_occupancy_cell(occupancy, occupancy_T, curr_pt.x, curr_pt.y,
                            delta, atomic_updates);
    }
  }
}

// record the old route's segments once so we can check overlap once per
// candidate segment, instead of checking every square against the old wire
struct OldRouteSegment {
  bool horizontal;
  int fixed_coord, low, high;
};

struct OldRouteSegments {
  int num_segments;
  OldRouteSegment segments[4];
};

void build_old_route_segments(const Wire &wire, OldRouteSegments &old_route_segments) {
  old_route_segments.num_segments = wire.num_bends + 1;
  Point curr_pt = get_point_i(wire, 0);
  for (int s = 0; s < old_route_segments.num_segments; s++) {
    Point next_pt = get_point_i(wire, s + 1);
    OldRouteSegment &segment = old_route_segments.segments[s];
    segment.horizontal = (curr_pt.y == next_pt.y);
    if (segment.horizontal) {
      segment.fixed_coord = curr_pt.y;
      segment.low = std::min(curr_pt.x, next_pt.x);
      segment.high = std::max(curr_pt.x, next_pt.x);
    } else {
      segment.fixed_coord = curr_pt.x;
      segment.low = std::min(curr_pt.y, next_pt.y);
      segment.high = std::max(curr_pt.y, next_pt.y);
    }
    curr_pt = next_pt;
  }
}

/** gets cost of laying a route on the occupancy matrix without
 calculating all square costs by only looking at the val at each square the route
 passes through. add (occupancy + 1)^2 to that square and keep traversing
*/
int calculate_route_cost(const std::vector<std::vector<int>> &occupancy,
                         const std::vector<std::vector<int>> &occupancy_T,
                         const Wire &route) {
  Point curr_pt = get_point_i(route, 0);
  int new_count = occupancy[curr_pt.y][curr_pt.x] + 1;
  int route_cost = new_count * new_count;

  // loop through all segments and add the cost of every point
  const int num_segments = route.num_bends + 1;
  for (int s = 0; s < num_segments; s++) {
    // get end point for curr segment
    Point next_pt = get_point_i(route, s + 1);

    if (next_pt.y == curr_pt.y) {
      // horizontal segment: walk along one row of occupancy
      const int *row = occupancy[curr_pt.y].data();
      int step_x = (next_pt.x > curr_pt.x) ? 1 : -1;
      for (int x = curr_pt.x + step_x; x != next_pt.x + step_x; x += step_x) {
        new_count = row[x] + 1;
        route_cost += new_count * new_count;
      }
    } else {
      // vertical segment: walk along one row of occupancy_T (= one column of occupancy)
      const int *col = occupancy_T[curr_pt.x].data();
      int step_y = (next_pt.y > curr_pt.y) ? 1 : -1;
      for (int y = curr_pt.y + step_y; y != next_pt.y + step_y; y += step_y) {
        new_count = col[y] + 1;
        route_cost += new_count * new_count;
      }
    }
    // this segment's end is the next segment's start
    curr_pt = next_pt;
  }
  return route_cost;
}

// walk every square in this part of the segment. adjustment is 0 where
// the old wire is already present, and 1 if new wire goes through but old wire doesn't
inline int calculate_segment_cost_across(const int *line, int low, int high,
                                          int adjustment) {
  int cost = 0;
  for (int p = low; p <= high; p++) {
    int count;
    #pragma omp atomic read
    count = line[p]; // atomic read
    int new_count = count + adjustment;
    cost += new_count * new_count;
    // cost is just the square of the wire count (after adjustment)
  }
  return cost;
}

// same scoring rule for across wires. find overlap with the old route once
// per segment, then still read and score every square along the candidate
int calculate_route_cost_across(const std::vector<std::vector<int>> &occupancy,
                              const std::vector<std::vector<int>> &occupancy_T,
                              const Wire &route,
                              const OldRouteSegments &old_route_segments) {
  Point curr_pt = get_point_i(route, 0);
  int count;
  #pragma omp atomic read
  count = occupancy[curr_pt.y][curr_pt.x];
  // the start is always on the old wire, so removing then adding it cancels
  int route_cost = count * count;

  const int num_segments = route.num_bends + 1;
  // LOOP 1 GOAL: walk the wire one straight piece at a time
  for (int s = 0; s < num_segments; s++) {
    Point next_pt = get_point_i(route, s + 1);
    bool horizontal = (curr_pt.y == next_pt.y);
    int fixed_coord, low, high;
    const int *line;
    if (horizontal) {
      fixed_coord = curr_pt.y;
      // line: specific row or col that the segment is part of, along with start and end points
      // low and high are the segments end x or y
      // it doesnt have to work out occupancy[y][x] at every square
      // its for cache
      line = occupancy[curr_pt.y].data();
      low = std::min(curr_pt.x, next_pt.x);
      high = std::max(curr_pt.x, next_pt.x);
      // don't count the segment's starting point again
      if (curr_pt.x <= next_pt.x) {
        low++;
      } else {
        high--;
      }
    } else {
      fixed_coord = curr_pt.x;
      line = occupancy_T[curr_pt.x].data();
      low = std::min(curr_pt.y, next_pt.y);
      high = std::max(curr_pt.y, next_pt.y);
      if (curr_pt.y <= next_pt.y) {
        low++;
      } else {
        high--;
      }
    }

    // legal routes move monotonically towards the end so old wire
    // intersection with row or column is a whole interval
    int old_low = std::numeric_limits<int>::max();
    int old_high = std::numeric_limits<int>::min();
    // LOOP 2 GOAL: checking if old wire is on this line (for segemnt currently investigating)
    // goes through EVERY SEGMENT of OLD WIRE
    // if it is, no +1
    // this gets you the range 
    for (int j = 0; j < old_route_segments.num_segments; j++) {
      const OldRouteSegment &old_segment = old_route_segments.segments[j];
      if (old_segment.horizontal == horizontal) {
        // fixed_coord is the row or col number it stays on
        // so if thats the same as our current one, then the segments MIGHT overlap
        if (old_segment.fixed_coord == fixed_coord) {
          old_low = std::min(old_low, old_segment.low);
          old_high = std::max(old_high, old_segment.high);
        }
      } else if (fixed_coord >= old_segment.low && fixed_coord <= old_segment.high) {
        old_low = std::min(old_low, old_segment.fixed_coord);
        old_high = std::max(old_high, old_segment.fixed_coord);
      }
    }
    // segment intersection code - end goal is to get cost
    int overlap_low = std::max(low, old_low);
    int overlap_high = std::min(high, old_high);
    if (overlap_low <= overlap_high) {
      // outside the old route, add 1; inside it, subtract 1 then add 1
      route_cost += calculate_segment_cost_across(line, low, overlap_low - 1, 1);
      route_cost += calculate_segment_cost_across(line, overlap_low, overlap_high, 0);
      route_cost += calculate_segment_cost_across(line, overlap_high + 1, high, 1);
    } else {
      route_cost += calculate_segment_cost_across(line, low, high, 1);
    }
    curr_pt = next_pt;
  }
  return route_cost;
}

// we are building the route based on the index, bounding box 
Wire build_route(const Wire &wire, int route_idx) {
  Wire route = wire;
  const int start_x = wire.start_x;
  const int start_y = wire.start_y;
  const int end_x = wire.end_x;
  const int end_y = wire.end_y;
  const int dx = std::abs(end_x - start_x);
  const int dy = std::abs(end_y - start_y);
  // direction the start pt points towards the end pt if we think about
  // it directionally
  const int x_dir = (end_x > start_x) ? 1 : -1;
  const int y_dir = (end_y > start_y) ? 1 : -1;
  const int routes_per_three_bend_shape = (dx - 1) * (dy - 1);

  // low id chooses horizontal first options: turn at column x1
  if (route_idx < dx) {
    int x1 = start_x + x_dir * (route_idx + 1);
    if (x1 == end_x) {
      route.num_bends = 1;
      route.bends[0] = {end_x, start_y};
    } else {
      route.num_bends = 2;
      route.bends[0] = {x1, start_y};
      route.bends[1] = {x1, end_y};
    }
  // go down vertical first options with v height = 1, 2 ... turn at row y1 into horizontal
  } else if (route_idx < dx + dy) {
    int y1 = start_y + y_dir * (route_idx - dx + 1);
    if (y1 == end_y) {
      route.num_bends = 1;
      route.bends[0] = {start_x, end_y};
    } else {
      route.num_bends = 2;
      route.bends[0] = {start_x, y1};
      route.bends[1] = {end_x, y1};
    }
  } else {
    // pick the middle bend location randomly within the "bounding box"
    // and then set the other two bends to make the route work
    // k = which 3 bend im on, subtract dx + dy because get rid of 1 bend and 2 bend

    // 3-bend routes
    // position within the 3 bend section
    int three_bend_offset = route_idx - (dx + dy);

    // first half of the section is horizontal first, second half is vertical first
    bool starts_horizontal = (three_bend_offset < routes_per_three_bend_shape);
    if (!starts_horizontal) {
      three_bend_offset -= routes_per_three_bend_shape;
    }

    // unflatten the offset
    int turn_column_choice = three_bend_offset / (dy - 1);
    int turn_row_choice    = three_bend_offset % (dy - 1);
    int turn_x = start_x + x_dir * (turn_column_choice + 1);
    int turn_y = start_y + y_dir * (turn_row_choice + 1);

    route.num_bends = 3;
    if (starts_horizontal) {
      route.bends[0] = {turn_x, start_y};
      route.bends[1] = {turn_x, turn_y};
      route.bends[2] = {end_x, turn_y};
    } else {
      route.bends[0] = {start_x, turn_y};
      route.bends[1] = {turn_x, turn_y};
      route.bends[2] = {turn_x, end_y};
    }
  }
  return route;
}

void route_within_wires(std::vector<Wire> &wires, std::vector<std::vector<int>> &occupancy,
     std::vector<std::vector<int>> &occupancy_T,
     int SA_iters, double SA_prob, int num_threads) {
  // random number generator for simulated annealing
  std::mt19937 rng(676);
  std::bernoulli_distribution choose_random_route(SA_prob);

  // do iterative algorithm
  for (int iters = 0; iters < SA_iters; iters++) {
    for (auto &wire : wires) { // no parallelization across wires
      // straight line conditions
      if (wire.num_bends == 0) {
        continue;
      }
      // remove original route placed in greedy placement so it doesn't count the cost of that layout
      update_occupancy_with_wire(occupancy, occupancy_T, wire, -1);

      int dx = std::abs(wire.end_x - wire.start_x);
      int dy = std::abs(wire.end_y - wire.start_y);
      int num_routes = dx + dy + 2 * (dx - 1) * (dy - 1);

      Wire best_wire = wire;
      // either pick a random route or do a full search for the best route
      if (choose_random_route(rng)) {
        std::uniform_int_distribution<int> pick_random_route(0, num_routes - 1);
        best_wire = build_route(wire, pick_random_route(rng));
      } 
      else {
        // initialize worst cost to highest number so better costs can compare if <
        // https://en.cppreference.com/cpp/types/numeric_limits/max 
        int best_cost = std::numeric_limits<int>::max();
        int best_route_idx = num_routes;

        // parallelize the candidate search across threads (threads might have to deal with more than 1 route)
        // threads score every possible route in parallel, and each one updates its private best variables
        // here theres a fork (overhead preventing perfect speedup)
        // tells the compiler: “Divide the iterations of the following loop among the threads in this parallel region.”
        #pragma omp parallel num_threads(num_threads)
        {
          int local_best_cost = std::numeric_limits<int>::max();
          int local_best_idx = num_routes;

          #pragma omp for schedule(static) nowait
          for (int route_idx = 0; route_idx < num_routes; route_idx++) {

            // calculate a route given the enumeration strat we have with route_idx
            Wire new_route = build_route(wire, route_idx);
            int new_cost = calculate_route_cost(occupancy, occupancy_T, new_route);

            if (new_cost < local_best_cost || (new_cost == local_best_cost &&
                route_idx < local_best_idx)) { // ties go to the lowest index
              local_best_cost = new_cost;
              local_best_idx = route_idx;
            }
          }

          // now we have best, best result
          // among all the thread's private local bests pick the lowest cost, so we have the best amongst all threads
          #pragma omp critical
          {
            if (local_best_cost < best_cost ||
                (local_best_cost == best_cost &&
                local_best_idx < best_route_idx)) {
              best_cost = local_best_cost;
              best_route_idx = local_best_idx;
            }
          }
        }
        best_wire = build_route(wire, best_route_idx);
      }
      // actually update the occupancy matrix with the new best wire
      wire = best_wire;
      update_occupancy_with_wire(occupancy, occupancy_T, wire, 1);
    }
  }
}

/** counts how many legal routes a bent wire has (all 1, 2, and 3 bend routes) */
int count_routes(const Wire &wire) {
  int dx = std::abs(wire.end_x - wire.start_x);
  int dy = std::abs(wire.end_y - wire.start_y);
  return dx + dy + 2 * (dx - 1) * (dy - 1);
}

/** sequentially try every route for one wire and returns the index of the
 cheapest one (for across wires approach) */
int find_best_route_idx(const std::vector<std::vector<int>> &occupancy,
                        const std::vector<std::vector<int>> &occupancy_T,
                        const Wire &wire, int num_routes,
                        const OldRouteSegments *old_route_segments = nullptr) {
  // found from https://en.cppreference.com/cpp/types/numeric_limits/max
  // initialize worst cost to highest number so better costs can compare if <
  int best_cost = std::numeric_limits<int>::max();
  int best_route_idx = 0;
  for (int route_idx = 0; route_idx < num_routes; route_idx++) {
    Wire new_route = build_route(wire, route_idx);
    int new_cost;
    new_cost = calculate_route_cost(occupancy, occupancy_T , new_route);
    
    if (new_cost < best_cost) { // < ties go to the lowest index
      best_cost = new_cost;
      best_route_idx = route_idx;
    }
  }
  return best_route_idx;
}

/** across wires: wires are processed in batches of batch_size.
 each thread handles its own batch of wires sequentially.
 logic:
 - each thread makes the simulated annealing random choices for its batch
  up front to prevent wasting path traversal work
 - the thread searches for best routes without changing occupancy yet
 - after all choices are ready, the thread writes every wire in its batch
  back on the grid using atomic updates
*/
void route_across_wires_MPI(std::vector<Wire> &wires, std::vector<std::vector<int>> &occupancy,
  std::vector<std::vector<int>> &occupancy_T,
  int SA_iters, double SA_prob, int num_threads, int batch_size) {
  const int num_wires = wires.size();

  // do iterative algorithm
  for (int iters = 0; iters < SA_iters; iters++) {
    for (int batch_start = 0; batch_start < num_wires; batch_start += batch_size) {
      // make sure final batch is right size
      int batch_end = std::min(batch_start + batch_size, num_wires);

      // random number generator for simulated annealing, private to this batch
      std::seed_seq seed{676, iters, batch_start};
      std::mt19937 rng(seed);
      std::bernoulli_distribution choose_random_route(SA_prob);
      // decide randomness up front but leave old routes on the grid
      // wires in the same batch don't see each other's new routes until all
      // choices are ready, so there can be staler information
      for (int w = batch_start; w < batch_end; w++) {
        // batch_idx is wire's position in current batch
        int batch_idx = w - batch_start;
        // reset to false
        search_curr_wire[batch_idx] = false;
        if (wires[w].num_bends == 0) {
          continue;
        }
        Wire old_wire = wires[w];
        update_occupancy_with_wire(occupancy, occupancy_T, old_wire, -1);
        int num_routes = count_routes(old_wire);
        int chosen_idx;
        
        // stimulated annealing
        if (choose_random_route(rng)) {
          std::uniform_int_distribution<int> pick(0, num_routes - 1);
          chosen_idx = pick(rng);
        } else {
          chosen_idx = find_best_route_idx(occupancy, occupancy_T, old_wire, num_routes);
        }

        // place best "chosen" wire
        wires[w] = build_route(old_wire, chosen_idx);
        update_occupancy_with_wire(occupancy, occupancy_T, wires[w], +1);
      }
    }
  }
}

////////////// END HELPERS /////////////////////

int main(int argc, char *argv[]) {
  const auto init_start = std::chrono::steady_clock::now();
  int pid;
  int nproc;

  // Initialize MPI
  MPI_Init(&argc, &argv);
  // Get process rank
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  // Get total number of processes  
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  std::string input_filename;
  double SA_prob = 0.1;
  int SA_iters = 5;
  int batch_size = 1;

  // Read command line arguments
  int opt;
  while ((opt = getopt(argc, argv, "f:p:i:m:b:")) != -1) {
    switch (opt) {
      case 'f':
        input_filename = optarg;
        break;
      case 'p':
        SA_prob = atof(optarg);
        break;
      case 'i':
        SA_iters = atoi(optarg);
        break;
      case 'b':
        batch_size = atoi(optarg);
        break;
      default:
        if (pid == 0) {
          std::cerr << "Usage: " << argv[0] << " -f input_filename [-p SA_prob] [-i SA_iters] -b batch_size\n";
        }

        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
  }

  // Check if required options are provided
  if (empty(input_filename) || SA_iters <= 0 || batch_size <= 0) {
    if (pid == 0) {
      std::cerr << "Usage: " << argv[0] << " -f input_filename [-p SA_prob] [-i SA_iters] -b batch_size\n";
    }

    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  if (pid == 0) {
    std::cout << "Number of processes: " << nproc << '\n';
    std::cout << "Simulated annealing probability parameter: " << SA_prob << '\n';
    std::cout << "Simulated annealing iterations: " << SA_iters << '\n';
    std::cout << "Input file: " << input_filename << '\n';
    std::cout << "Batch size: " << batch_size << '\n';
  }

  int dim_x, dim_y, num_wires;
  std::vector<Wire> wires;
  std::vector<std::vector<int>> occupancy;

  if (pid == 0) {
      std::ifstream fin(input_filename);

      if (!fin) {
        std::cerr << "Unable to open file: " << input_filename << ".\n";
        exit(EXIT_FAILURE);
      }

      /* Read the grid dimension and wire information from file */
      fin >> dim_x >> dim_y >> num_wires;

      wires.resize(num_wires);
      for (auto &wire : wires) {
        fin >> wire.start_x >> wire.start_y >> wire.end_x >> wire.end_y;
        wire.move_x_start = true;
        wire.move_x_end = false;
        wire.mid_x = wire.end_x;
        wire.mid_y = wire.start_y;
      }
  }

  /* Initialize any additional data structures needed in the algorithm */

  if (pid == 0) {
    const double init_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - init_start).count();
    std::cout << "Initialization time (sec): " << std::fixed << std::setprecision(10) << init_time << '\n';
  }

  const auto compute_start = std::chrono::steady_clock::now();

  /** 
   * (TODO)
   * Implement the wire routing algorithm here
   * Feel free to structure the algorithm into different functions
   * Use MPI to parallelize the algorithm. 
   */
  
  // TODO WE HAVE TO PARSE THE INPUT WIRE START AND END AND POPULATE THE
  // WIRES list, making sure we 0 initiate all parameters
  // place initial route using greedy heuristic
  for (auto &wire : wires){ // OPENMP PARALLELIZE: each wire is independent for this initial laying
    // straight line conditions
    if (wire.start_x == wire.end_x) {
      wire.num_bends = 0;
      // draw straight vertical line
      update_occupancy_with_wire(occupancy, occupancy_T, wire, 1);
    } else if (wire.start_y == wire.end_y) {
      // draw straight horizontal line
      wire.num_bends = 0;
      update_occupancy_with_wire(occupancy, occupancy_T, wire, 1);
    } else {
      // naive so draw 1 bend line (just to lay smth down)
      wire.bends[0] = {wire.start_x, wire.end_y};
      wire.num_bends = 1;
      update_occupancy_with_wire(occupancy, occupancy_T, wire, 1);
    }
  }
  // initialize wires
  // Within wires
  route_across_wires(wires, occupancy, occupancy_T, SA_iters, SA_prob, num_threads, batch_size);


  if (pid == 0) {
    const double compute_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - compute_start).count();
    std::cout << "Computation time (sec): " << std::fixed << std::setprecision(10) << compute_time << '\n';
  }

  if (pid == 0) {
    wr_checker checker(wires, occupancy);
    checker.validate();

    /* Write wires and occupancy matrix to files */
    print_stats(occupancy);
    write_output(wires, num_wires, occupancy, dim_x, dim_y);
  }

  // Cleanup
  MPI_Finalize();
}


/* TODO (student): implement to_validate_format to convert Wire to
  validate_wire_t keypoint representation in order to run checker and
  write output
*/
validate_wire_t Wire::to_validate_format(void) const {
  validate_wire_t w;
  
  return w;
}
