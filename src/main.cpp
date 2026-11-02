#include <Eigen/Dense>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef HAVE_GUROBI
#include <gurobi_c++.h>
#endif

namespace {

struct Options {
    std::string input;
    std::string mode = "joint";
    std::string bias_mode = "direct";
    std::string core_solver = "bsb";
    std::string row_set = "free";
    std::string preset = "paper";
    std::string tie_mode = "zero";
    std::string field_schedule = "coupled";
    std::string dynamics = "ballistic";
    std::string normalization = "full";
    int n_free = 4;
    int rounds = 5;
    int max_partitions = 1000;
    int max_steps = 2000;
    int sample_period = 0;
    int energy_window = 0;
    int local_sweeps = 0;
    int restarts = 1;
    int refine_restarts = 0;
    int refine_sweeps = 20;
    int swap_restarts = 0;
    int swap_sweeps = 10;
    int anneal_restarts = 0;
    int anneal_steps = 0;
    float anneal_temperature = 0.1F;
    int row_sa_restarts = 0;
    int row_sa_steps = 0;
    float row_sa_temperature = 0.1F;
    float field_scale = 1.0F;
    float separate_slack = 0.0F;
    float oracle_time = 0.1F;
    int greedy_sweeps = 0;
    unsigned seed = 1;
    float epsilon = 1.0e-8F;
    float time_step = 0.1F;
    float pressure_slope = 0.01F;
    float initial_scale = 1.0F;
    float momentum_scale = 1.0F;
    // The paper omits the exact normalization used around the standard
    // random-matrix c0 rule.  0.3 reproduces its Table-1 quality/time tradeoff.
    float coupling_scale = 0.3F;
    bool quiet = false;
    bool paper_plus = false;
};

struct Partition {
    std::vector<int> free;
    std::vector<int> bound;
};

struct MatrixView {
    int rows = 0;
    int cols = 0;
    std::vector<int> row_of_address;
    std::vector<int> col_of_address;
    std::vector<int> address_of_cell;
};

struct ColumnSolution {
    std::vector<uint8_t> v1;
    std::vector<uint8_t> v2;
    std::vector<uint8_t> type;
    double objective = std::numeric_limits<double>::infinity();
    int steps = 0;
};

std::unordered_map<std::string, std::string> parse_arguments(int argc, char ** argv) {
    std::unordered_map<std::string, std::string> values;
    for (int i = 1; i < argc; ++i) {
        std::string key(argv[i]);
        if (key == "--quiet") {
            values[key] = "1";
            continue;
        }
        if (key.rfind("--", 0) != 0 || i + 1 >= argc)
            throw std::invalid_argument("expected --key value, got: " + key);
        values[key] = argv[++i];
    }
    return values;
}
Options make_options(int argc, char ** argv) {
    const auto args = parse_arguments(argc, argv);
    Options opt;
    auto set_string = [&](const char * key, std::string & target) {
        if (auto it = args.find(key); it != args.end()) target = it->second;
    };
    auto set_int = [&](const char * key, int & target) {
        if (auto it = args.find(key); it != args.end()) target = std::stoi(it->second);
    };
    auto set_uint = [&](const char * key, unsigned & target) {
        if (auto it = args.find(key); it != args.end()) target =
            static_cast<unsigned>(std::stoul(it->second));
    };
    auto set_float = [&](const char * key, float & target) {
        if (auto it = args.find(key); it != args.end()) target = std::stof(it->second);
    };

    // The paper preset uses a fixed field calibration for each formulation,
    // so parse mode before applying preset defaults. Explicit CLI parameters
    // below can still override every preset choice.
    set_string("--input", opt.input);
    set_string("--mode", opt.mode);
    set_string("--preset", opt.preset);
    if (opt.preset == "paper" || opt.preset == "paper-plus") {
        opt.row_set = "free";
        if (opt.mode == "separate") {
            opt.field_schedule = "sqrt";
            opt.field_scale = 0.7F;
        } else {
            opt.field_schedule = "coupled-linear";
            opt.field_scale = 0.329F;
        }
        opt.normalization = "full";
        opt.seed = 1;
        if (opt.preset == "paper-plus") {
            opt.paper_plus = true;
        }
    } else {
        throw std::invalid_argument(
            "--preset must be paper or paper-plus");
    }
    set_string("--bias-mode", opt.bias_mode);
    set_string("--core-solver", opt.core_solver);
    set_string("--row-set", opt.row_set);
    set_string("--tie-mode", opt.tie_mode);
    set_string("--field-schedule", opt.field_schedule);
    set_string("--dynamics", opt.dynamics);
    set_string("--normalization", opt.normalization);
    set_int("--n-free", opt.n_free);
    set_int("--rounds", opt.rounds);
    set_int("--partitions", opt.max_partitions);
    set_int("--max-steps", opt.max_steps);
    set_int("--sample-period", opt.sample_period);
    set_int("--energy-window", opt.energy_window);
    set_int("--local-sweeps", opt.local_sweeps);
    set_int("--restarts", opt.restarts);
    set_int("--refine-restarts", opt.refine_restarts);
    set_int("--refine-sweeps", opt.refine_sweeps);
    set_int("--swap-restarts", opt.swap_restarts);
    set_int("--swap-sweeps", opt.swap_sweeps);
    set_int("--anneal-restarts", opt.anneal_restarts);
    set_int("--anneal-steps", opt.anneal_steps);
    set_float("--anneal-temperature", opt.anneal_temperature);
    set_int("--row-sa-restarts", opt.row_sa_restarts);
    set_int("--row-sa-steps", opt.row_sa_steps);
    set_float("--row-sa-temperature", opt.row_sa_temperature);
    set_float("--field-scale", opt.field_scale);
    set_float("--separate-slack", opt.separate_slack);
    set_float("--oracle-time", opt.oracle_time);
    set_int("--greedy-sweeps", opt.greedy_sweeps);
    set_uint("--seed", opt.seed);
    set_float("--epsilon", opt.epsilon);
    set_float("--time-step", opt.time_step);
    set_float("--pressure-slope", opt.pressure_slope);
    set_float("--initial-scale", opt.initial_scale);
    set_float("--momentum-scale", opt.momentum_scale);
    set_float("--coupling-scale", opt.coupling_scale);
    opt.quiet = args.contains("--quiet");

    if (opt.input.empty()) throw std::invalid_argument("--input is required");
    if (opt.mode != "separate" && opt.mode != "joint")
        throw std::invalid_argument("--mode must be separate or joint");
    if (opt.bias_mode != "direct" && opt.bias_mode != "auxiliary")
        throw std::invalid_argument("--bias-mode must be direct or auxiliary");
    if (opt.core_solver != "bsb" && opt.core_solver != "gurobi")
        throw std::invalid_argument("--core-solver must be bsb or gurobi");
    if (opt.row_set != "free" && opt.row_set != "bound")
        throw std::invalid_argument("--row-set must be free or bound");
    if (opt.tie_mode != "zero" && opt.tie_mode != "preserve" &&
        opt.tie_mode != "random")
        throw std::invalid_argument("--tie-mode must be zero, preserve, or random");
    if (opt.field_schedule != "coupled" &&
        opt.field_schedule != "coupled-linear" &&
        opt.field_schedule != "linear" &&
        opt.field_schedule != "sqrt")
        throw std::invalid_argument(
            "--field-schedule must be coupled, coupled-linear, linear, or sqrt");
    if (opt.dynamics != "ballistic" && opt.dynamics != "discrete")
        throw std::invalid_argument("--dynamics must be ballistic or discrete");
    if (opt.normalization != "full" && opt.normalization != "pairwise")
        throw std::invalid_argument("--normalization must be full or pairwise");
#ifndef HAVE_GUROBI
    if (opt.core_solver == "gurobi")
        throw std::invalid_argument("this build has no Gurobi support");
#endif
    return opt;
}

std::vector<uint64_t> load_table(const std::string & path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open input table: " + path);
    std::vector<uint64_t> table;
    int address = 0;
    uint64_t value = 0;
    while (input >> address >> value) {
        if (address != static_cast<int>(table.size()))
            throw std::runtime_error("truth-table addresses must be contiguous from zero");
        table.push_back(value);
    }
    if (table.empty() || (table.size() & (table.size() - 1)) != 0)
        throw std::runtime_error("truth-table length must be a nonzero power of two");
    return table;
}

int exact_log2(std::size_t value) {
    int result = 0;
    while ((std::size_t{1} << result) < value) ++result;
    if ((std::size_t{1} << result) != value)
        throw std::invalid_argument("value is not a power of two");
    return result;
}

int output_bits(const std::vector<uint64_t> & table) {
    const auto maximum = *std::max_element(table.begin(), table.end());
    int bits = 1;
    while ((uint64_t{1} << bits) <= maximum) ++bits;
    return bits;
}

void generate_combinations(int n, int k, int start, std::vector<int> & current,
                           std::vector<Partition> & result) {
    if (static_cast<int>(current.size()) == k) {
        Partition partition;
        partition.free = current;
        std::vector<uint8_t> selected(n, 0);
        for (int value : current) selected[value] = 1;
        for (int value = 0; value < n; ++value)
            if (!selected[value]) partition.bound.push_back(value);
        result.push_back(std::move(partition));
        return;
    }
    for (int value = start; value <= n - (k - static_cast<int>(current.size())); ++value) {
        current.push_back(value);
        generate_combinations(n, k, value + 1, current, result);
        current.pop_back();
    }
}

std::vector<Partition> make_partitions(int n, int n_free, int maximum,
                                       std::mt19937 & random) {
    std::vector<Partition> all;
    std::vector<int> current;
    generate_combinations(n, n_free, 0, current, all);
    if (maximum <= 0 || static_cast<int>(all.size()) <= maximum) return all;
    std::shuffle(all.begin(), all.end(), random);
    all.resize(maximum);
    return all;
}

MatrixView make_matrix_view(int n_inputs, const Partition & partition,
                            const std::string & row_set) {
    MatrixView view;
    const std::vector<int> & row_variables =
        row_set == "free" ? partition.free : partition.bound;
    const std::vector<int> & col_variables =
        row_set == "free" ? partition.bound : partition.free;
    view.rows = 1 << static_cast<int>(row_variables.size());
    view.cols = 1 << static_cast<int>(col_variables.size());
    const int addresses = 1 << n_inputs;
    view.row_of_address.resize(addresses);
    view.col_of_address.resize(addresses);
    view.address_of_cell.resize(addresses);
    for (int address = 0; address < addresses; ++address) {
        int row = 0;
        int col = 0;
        for (int index = 0; index < static_cast<int>(row_variables.size()); ++index)
            row |= ((address >> row_variables[index]) & 1) << index;
        for (int index = 0; index < static_cast<int>(col_variables.size()); ++index)
            col |= ((address >> col_variables[index]) & 1) << index;
        view.row_of_address[address] = row;
        view.col_of_address[address] = col;
        view.address_of_cell[row * view.cols + col] = address;
    }
    return view;
}

double mean_error_distance(const std::vector<uint64_t> & exact,
                           const std::vector<uint64_t> & approximate) {
    double sum = 0.0;
    for (std::size_t i = 0; i < exact.size(); ++i) {
        const auto lhs = static_cast<int64_t>(exact[i]);
        const auto rhs = static_cast<int64_t>(approximate[i]);
        sum += static_cast<double>(std::llabs(lhs - rhs));
    }
    return sum / static_cast<double>(exact.size());
}

double candidate_objective(const std::string & mode, int bit,
                           const MatrixView & view,
                           const std::vector<uint64_t> & exact,
                           const std::vector<uint64_t> & approximate_without_bit,
                           const std::vector<uint8_t> & v1,
                           const std::vector<uint8_t> & v2,
                           std::vector<uint8_t> * optimal_type,
                           const std::vector<uint8_t> * tie_type = nullptr) {
    std::vector<double> error1(view.cols, 0.0);
    std::vector<double> error2(view.cols, 0.0);
    const uint64_t weight = uint64_t{1} << bit;
    for (int address = 0; address < static_cast<int>(exact.size()); ++address) {
        const int row = view.row_of_address[address];
        const int col = view.col_of_address[address];
        if (mode == "separate") {
            const int target = static_cast<int>((exact[address] >> bit) & 1U);
            error1[col] += v1[row] != target;
            error2[col] += v2[row] != target;
        } else {
            const int64_t target = static_cast<int64_t>(exact[address]);
            const int64_t base = static_cast<int64_t>(approximate_without_bit[address]);
            error1[col] += std::llabs(target - (base + static_cast<int64_t>(weight * v1[row])));
            error2[col] += std::llabs(target - (base + static_cast<int64_t>(weight * v2[row])));
        }
    }
    double total = 0.0;
    if (optimal_type) optimal_type->assign(view.cols, 0);
    for (int col = 0; col < view.cols; ++col) {
        const bool tied = std::abs(error2[col] - error1[col]) < 1.0e-12;
        const bool use_second = error2[col] < error1[col] ||
            (tied && tie_type && (*tie_type)[col]);
        if (optimal_type) (*optimal_type)[col] = static_cast<uint8_t>(use_second);
        total += use_second ? error2[col] : error1[col];
    }
    return total / static_cast<double>(exact.size());
}

double polish_patterns(const std::string & mode, int bit,
                       const MatrixView & view,
                       const std::vector<uint64_t> & exact,
                       const std::vector<uint64_t> & approximate_without_bit,
                       int sweeps, std::vector<uint8_t> & v1,
                       std::vector<uint8_t> & v2,
                       std::vector<uint8_t> & type) {
    double objective = candidate_objective(mode, bit, view, exact,
                                           approximate_without_bit, v1, v2, &type);
    const int64_t weight = int64_t{1} << bit;
    for (int sweep = 0; sweep < sweeps; ++sweep) {
        const double previous = objective;
        for (int row = 0; row < view.rows; ++row) {
            double first[2] = {0.0, 0.0};
            double second[2] = {0.0, 0.0};
            for (int col = 0; col < view.cols; ++col) {
                const int address = view.address_of_cell[row * view.cols + col];
                double * errors = type[col] ? second : first;
                if (mode == "separate") {
                    const int target = static_cast<int>((exact[address] >> bit) & 1U);
                    errors[0] += target != 0;
                    errors[1] += target != 1;
                } else {
                    const int64_t target = static_cast<int64_t>(exact[address]);
                    const int64_t base = static_cast<int64_t>(approximate_without_bit[address]);
                    errors[0] += std::llabs(target - base);
                    errors[1] += std::llabs(target - (base + weight));
                }
            }
            v1[row] = static_cast<uint8_t>(first[1] < first[0]);
            v2[row] = static_cast<uint8_t>(second[1] < second[0]);
        }
        objective = candidate_objective(mode, bit, view, exact,
                                         approximate_without_bit, v1, v2, &type);
        if (objective >= previous - 1.0e-12) break;
    }
    return objective;
}

void optimal_patterns_for_type(
    const std::string & mode, int bit, const MatrixView & view,
    const std::vector<uint64_t> & exact,
    const std::vector<uint64_t> & approximate_without_bit,
    const std::vector<uint8_t> & type,
    std::vector<uint8_t> & v1, std::vector<uint8_t> & v2) {
    v1.assign(view.rows, 0);
    v2.assign(view.rows, 0);
    const int64_t weight = int64_t{1} << bit;
    for (int row = 0; row < view.rows; ++row) {
        double cost[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
        for (int col = 0; col < view.cols; ++col) {
            const int address = view.address_of_cell[row * view.cols + col];
            const int cluster = type[col] ? 1 : 0;
            if (mode == "separate") {
                const int target = static_cast<int>((exact[address] >> bit) & 1U);
                cost[cluster][0] += target != 0;
                cost[cluster][1] += target != 1;
            } else {
                const int64_t target = static_cast<int64_t>(exact[address]);
                const int64_t base = static_cast<int64_t>(approximate_without_bit[address]);
                cost[cluster][0] += std::llabs(target - base);
                cost[cluster][1] += std::llabs(target - (base + weight));
            }
        }
        v1[row] = static_cast<uint8_t>(cost[0][1] < cost[0][0]);
        v2[row] = static_cast<uint8_t>(cost[1][1] < cost[1][0]);
    }
}

double polish_type_swaps(
    const std::string & mode, int bit, const MatrixView & view,
    const std::vector<uint64_t> & exact,
    const std::vector<uint64_t> & approximate_without_bit,
    int sweeps, std::vector<uint8_t> & v1, std::vector<uint8_t> & v2,
    std::vector<uint8_t> & type) {
    optimal_patterns_for_type(mode, bit, view, exact, approximate_without_bit,
                              type, v1, v2);
    double objective = candidate_objective(mode, bit, view, exact,
                                            approximate_without_bit,
                                            v1, v2, &type);
    for (int sweep = 0; sweep < sweeps; ++sweep) {
        bool improved = false;
        for (int col = 0; col < view.cols; ++col) {
            std::vector<uint8_t> trial_type = type;
            trial_type[col] ^= 1U;
            std::vector<uint8_t> trial_v1, trial_v2, optimal_type;
            optimal_patterns_for_type(mode, bit, view, exact,
                                      approximate_without_bit, trial_type,
                                      trial_v1, trial_v2);
            const double trial_objective = candidate_objective(
                mode, bit, view, exact, approximate_without_bit,
                trial_v1, trial_v2, &optimal_type);
            if (trial_objective + 1.0e-12 < objective) {
                objective = trial_objective;
                v1 = std::move(trial_v1);
                v2 = std::move(trial_v2);
                type = std::move(optimal_type);
                improved = true;
            }
        }
        if (!improved) break;
    }
    return objective;
}

double solution_from_row_pattern(
    const std::string & mode, int bit, const MatrixView & view,
    const std::vector<uint64_t> & exact,
    const std::vector<uint64_t> & approximate_without_bit,
    const std::vector<uint8_t> & pattern,
    std::vector<uint8_t> & v1, std::vector<uint8_t> & v2) {
    v1.assign(view.rows, 0);
    v2.assign(view.rows, 0);
    const int64_t weight = int64_t{1} << bit;
    double total = 0.0;
    for (int row = 0; row < view.rows; ++row) {
        double cost[4] = {0.0, 0.0, 0.0, 0.0};
        for (int col = 0; col < view.cols; ++col) {
            const int address = view.address_of_cell[row * view.cols + col];
            double error0 = 0.0;
            double error1 = 0.0;
            if (mode == "separate") {
                const int target = static_cast<int>((exact[address] >> bit) & 1U);
                error0 = target != 0;
                error1 = target != 1;
            } else {
                const int64_t target = static_cast<int64_t>(exact[address]);
                const int64_t base = static_cast<int64_t>(approximate_without_bit[address]);
                error0 = std::llabs(target - base);
                error1 = std::llabs(target - (base + weight));
            }
            cost[0] += error0;
            cost[1] += error1;
            cost[2] += pattern[col] ? error1 : error0;
            cost[3] += pattern[col] ? error0 : error1;
        }
        const int type = static_cast<int>(
            std::min_element(cost, cost + 4) - cost);
        total += cost[type];
        if (type == 1) {
            v1[row] = 1;
            v2[row] = 1;
        } else if (type == 2) {
            v2[row] = 1;
        } else if (type == 3) {
            v1[row] = 1;
        }
    }
    return total / static_cast<double>(exact.size());
}

Eigen::MatrixXf build_coefficients(const std::string & mode, int bit,
                                   const MatrixView & view,
                                   const std::vector<uint64_t> & exact,
                                   const std::vector<uint64_t> & approximate_without_bit) {
    Eigen::MatrixXf coefficient(view.rows, view.cols);
    coefficient.setZero();
    const int64_t weight = int64_t{1} << bit;
    const float probability = 1.0F / static_cast<float>(exact.size());
    for (int address = 0; address < static_cast<int>(exact.size()); ++address) {
        const int row = view.row_of_address[address];
        const int col = view.col_of_address[address];
        if (mode == "separate") {
            const int target = static_cast<int>((exact[address] >> bit) & 1U);
            coefficient(row, col) = probability * static_cast<float>(1 - 2 * target);
        } else {
            const int64_t difference = static_cast<int64_t>(approximate_without_bit[address]) -
                                       static_cast<int64_t>(exact[address]);
            int64_t q = 0;
            if (-weight <= difference && difference <= 0)
                q = weight + 2 * difference;
            else
                q = weight * (difference >= 0 ? 1 : -1);
            coefficient(row, col) = probability * static_cast<float>(q);
        }
    }
    return coefficient;
}

double ising_energy(const Eigen::VectorXf & spin_v1,
                    const Eigen::VectorXf & spin_v2,
                    const Eigen::VectorXf & spin_type,
                    const Eigen::MatrixXf & coupling,
                    const Eigen::VectorXf & field_v) {
    return -static_cast<double>(field_v.dot(spin_v1) + field_v.dot(spin_v2))
           - static_cast<double>(spin_v1.dot(coupling * spin_type))
           + static_cast<double>(spin_v2.dot(coupling * spin_type));
}

ColumnSolution solve_bsb(const Options & opt, int bit, const MatrixView & view,
                         const std::vector<uint64_t> & exact,
                         const std::vector<uint64_t> & approximate_without_bit,
                         std::mt19937 & random) {
    const Eigen::MatrixXf coefficient =
        build_coefficients(opt.mode, bit, view, exact, approximate_without_bit);
    const Eigen::MatrixXf coupling = coefficient * 0.25F;
    const Eigen::VectorXf linear = coupling.rowwise().sum();
    const Eigen::VectorXf field_v = -linear;

    const bool auxiliary_bias = opt.bias_mode == "auxiliary";
    const int total_spins = 2 * view.rows + view.cols + (auxiliary_bias ? 1 : 0);
    // The zero-field bSB normalization uses sum(J_ij^2).  For an Ising
    // problem with linear fields, the fixed-spin embedding contributes
    // 2*sum(h_i^2).  There are two pattern vectors with the same field, hence
    // the 4*||field_v||^2 term below in both direct and auxiliary modes.
    double symmetric_j_squared =
        4.0 * static_cast<double>(coupling.squaredNorm());
    if (opt.normalization == "full")
        symmetric_j_squared += 4.0 * static_cast<double>(field_v.squaredNorm());
    const double rms = symmetric_j_squared > 0.0
        ? std::sqrt(symmetric_j_squared /
                    static_cast<double>(total_spins * std::max(1, total_spins - 1)))
        : 1.0;
    const float c0 = opt.coupling_scale * static_cast<float>(0.5 /
        (std::sqrt(static_cast<double>(total_spins)) * rms));

    std::uniform_real_distribution<float> initial_position(
        -opt.initial_scale, opt.initial_scale);
    std::uniform_real_distribution<float> initial_momentum(
        -opt.momentum_scale, opt.momentum_scale);
    Eigen::VectorXf x1(view.rows), x2(view.rows), xt(view.cols);
    Eigen::VectorXf y1(view.rows), y2(view.rows), yt(view.cols);
    float x_bias = 1.0F;
    float y_bias = 0.0F;
    if (auxiliary_bias) {
        x_bias = initial_position(random);
        y_bias = initial_momentum(random);
    }
    for (int i = 0; i < view.rows; ++i) {
        x1[i] = initial_position(random);
        x2[i] = initial_position(random);
        y1[i] = initial_momentum(random);
        y2[i] = initial_momentum(random);
    }
    for (int j = 0; j < view.cols; ++j) {
        xt[j] = initial_position(random);
        yt[j] = initial_momentum(random);
    }

    const int sample_period = opt.sample_period > 0 ? opt.sample_period
                                                    : (exact.size() <= 512 ? 20 : 10);
    const int energy_window = opt.energy_window > 0 ? opt.energy_window
                                                    : (exact.size() <= 512 ? 20 : 10);
    std::vector<double> energies;
    energies.reserve(energy_window);
    ColumnSolution best;
    Eigen::VectorXf force1(view.rows), force2(view.rows), force_t(view.cols);
    Eigen::VectorXf spin1(view.rows), spin2(view.rows), spin_t(view.cols);
    Eigen::VectorXf interaction1(view.rows), interaction2(view.rows);
    Eigen::VectorXf interaction_t(view.cols);
    std::uniform_int_distribution<int> tie_random_bit(0, 1);

    auto apply_walls = [](Eigen::VectorXf & position, Eigen::VectorXf & momentum) {
        for (Eigen::Index i = 0; i < position.size(); ++i) {
            if (std::abs(position[i]) > 1.0F) {
                position[i] = position[i] > 0.0F ? 1.0F : -1.0F;
                momentum[i] = 0.0F;
            }
        }
    };

    for (int step = 1; step <= opt.max_steps; ++step) {
        const float pressure = std::min(1.0F, opt.pressure_slope * static_cast<float>(step));
        if (opt.dynamics == "discrete") {
            for (int i = 0; i < view.rows; ++i) {
                interaction1[i] = x1[i] >= 0.0F ? 1.0F : -1.0F;
                interaction2[i] = x2[i] >= 0.0F ? 1.0F : -1.0F;
            }
            for (int j = 0; j < view.cols; ++j)
                interaction_t[j] = xt[j] >= 0.0F ? 1.0F : -1.0F;
        } else {
            interaction1 = x1;
            interaction2 = x2;
            interaction_t = xt;
        }
        force1.noalias() = coupling * interaction_t;
        force2.noalias() = -coupling * interaction_t;
        if (auxiliary_bias) {
            force1 += field_v * x_bias;
            force2 += field_v * x_bias;
        }
        force_t.noalias() = coupling.transpose() * (interaction1 - interaction2);
        const float force_bias = auxiliary_bias
            ? field_v.dot(interaction1 + interaction2) : 0.0F;

        float field_eta = c0;
        if (opt.field_schedule == "coupled-linear")
            field_eta = c0 * pressure;
        else if (opt.field_schedule == "linear")
            field_eta = pressure;
        else if (opt.field_schedule == "sqrt")
            field_eta = std::sqrt(pressure);
        field_eta *= opt.field_scale;
        y1 += (-(1.0F - pressure) * x1 + c0 * force1 +
               (auxiliary_bias ? 0.0F : field_eta) * field_v) * opt.time_step;
        y2 += (-(1.0F - pressure) * x2 + c0 * force2 +
               (auxiliary_bias ? 0.0F : field_eta) * field_v) * opt.time_step;
        yt += (-(1.0F - pressure) * xt + c0 * force_t) * opt.time_step;
        if (auxiliary_bias)
            y_bias += (-(1.0F - pressure) * x_bias + c0 * force_bias) * opt.time_step;
        x1 += y1 * opt.time_step;
        x2 += y2 * opt.time_step;
        xt += yt * opt.time_step;
        if (auxiliary_bias) x_bias += y_bias * opt.time_step;
        apply_walls(x1, y1);
        apply_walls(x2, y2);
        apply_walls(xt, yt);
        if (auxiliary_bias && std::abs(x_bias) > 1.0F) {
            x_bias = x_bias > 0.0F ? 1.0F : -1.0F;
            y_bias = 0.0F;
        }

        if (step % sample_period != 0) continue;

        std::vector<uint8_t> v1(view.rows), v2(view.rows), type;
        std::vector<uint8_t> tie_type(view.cols);
        const float bias_spin = auxiliary_bias && x_bias < 0.0F ? -1.0F : 1.0F;
        for (int i = 0; i < view.rows; ++i) {
            v1[i] = static_cast<uint8_t>(x1[i] * bias_spin >= 0.0F);
            v2[i] = static_cast<uint8_t>(x2[i] * bias_spin >= 0.0F);
            spin1[i] = v1[i] ? 1.0F : -1.0F;
            spin2[i] = v2[i] ? 1.0F : -1.0F;
        }
        for (int j = 0; j < view.cols; ++j) {
            if (opt.tie_mode == "random")
                tie_type[j] = static_cast<uint8_t>(tie_random_bit(random));
            else
                tie_type[j] = static_cast<uint8_t>(xt[j] * bias_spin >= 0.0F);
        }
        const double objective = opt.local_sweeps > 0
            ? polish_patterns(opt.mode, bit, view, exact, approximate_without_bit,
                              opt.local_sweeps, v1, v2, type)
            : candidate_objective(
                  opt.mode, bit, view, exact, approximate_without_bit, v1, v2,
                  &type, opt.tie_mode == "zero" ? nullptr : &tie_type);
        for (int j = 0; j < view.cols; ++j) {
            spin_t[j] = type[j] ? 1.0F : -1.0F;
            xt[j] = spin_t[j] * bias_spin;
            yt[j] = 0.0F;
        }
        if (objective < best.objective) {
            best.v1 = std::move(v1);
            best.v2 = std::move(v2);
            best.type = type;
            best.objective = objective;
            best.steps = step;
        }

        const double energy = ising_energy(spin1, spin2, spin_t, coupling, field_v);
        energies.push_back(energy);
        if (static_cast<int>(energies.size()) > energy_window) energies.erase(energies.begin());
        if (static_cast<int>(energies.size()) == energy_window) {
            const double mean = std::accumulate(energies.begin(), energies.end(), 0.0) /
                                static_cast<double>(energy_window);
            double variance = 0.0;
            for (double value : energies) variance += (value - mean) * (value - mean);
            variance /= static_cast<double>(energy_window);
            if (variance < opt.epsilon) break;
        }
    }
    for (int sweep = 0; sweep < opt.greedy_sweeps; ++sweep) {
        bool improved = false;
        for (int pattern = 0; pattern < 2; ++pattern) {
            auto & values = pattern == 0 ? best.v1 : best.v2;
            for (int row = 0; row < view.rows; ++row) {
                values[row] ^= 1U;
                std::vector<uint8_t> type;
                const double objective = candidate_objective(
                    opt.mode, bit, view, exact, approximate_without_bit,
                    best.v1, best.v2, &type);
                if (objective + 1.0e-12 < best.objective) {
                    best.objective = objective;
                    best.type = std::move(type);
                    improved = true;
                } else {
                    values[row] ^= 1U;
                }
            }
        }
        if (!improved) break;
    }
    std::uniform_int_distribution<int> random_bit(0, 1);
    for (int restart = 0; restart < opt.refine_restarts; ++restart) {
        std::vector<uint8_t> v1(view.rows), v2(view.rows), type;
        for (int row = 0; row < view.rows; ++row) {
            v1[row] = static_cast<uint8_t>(random_bit(random));
            v2[row] = static_cast<uint8_t>(random_bit(random));
        }
        const double objective = polish_patterns(
            opt.mode, bit, view, exact, approximate_without_bit,
            opt.refine_sweeps, v1, v2, type);
        if (objective < best.objective) {
            best.v1 = std::move(v1);
            best.v2 = std::move(v2);
            best.type = std::move(type);
            best.objective = objective;
            best.steps = -restart - 1;
        }
    }
    for (int restart = 0; restart < opt.swap_restarts; ++restart) {
        std::vector<uint8_t> type(view.cols);
        if (restart == 0 && !best.type.empty()) {
            type = best.type;
        } else {
            for (int col = 0; col < view.cols; ++col)
                type[col] = static_cast<uint8_t>(random_bit(random));
        }
        std::vector<uint8_t> v1, v2;
        const double objective = polish_type_swaps(
            opt.mode, bit, view, exact, approximate_without_bit,
            opt.swap_sweeps, v1, v2, type);
        if (objective < best.objective) {
            best.v1 = std::move(v1);
            best.v2 = std::move(v2);
            best.type = std::move(type);
            best.objective = objective;
            best.steps = -opt.refine_restarts - restart - 1;
        }
    }
    std::uniform_real_distribution<double> probability(0.0, 1.0);
    std::uniform_int_distribution<int> pattern_index(0, 2 * view.rows - 1);
    for (int restart = 0; restart < opt.anneal_restarts; ++restart) {
        std::vector<uint8_t> v1 = best.v1;
        std::vector<uint8_t> v2 = best.v2;
        if (restart > 0) {
            for (int row = 0; row < view.rows; ++row) {
                v1[row] = static_cast<uint8_t>(random_bit(random));
                v2[row] = static_cast<uint8_t>(random_bit(random));
            }
        }
        std::vector<uint8_t> type;
        double objective = candidate_objective(
            opt.mode, bit, view, exact, approximate_without_bit,
            v1, v2, &type);
        for (int step = 0; step < opt.anneal_steps; ++step) {
            const int index = pattern_index(random);
            auto & pattern = index < view.rows ? v1 : v2;
            const int row = index < view.rows ? index : index - view.rows;
            pattern[row] ^= 1U;
            std::vector<uint8_t> trial_type;
            const double trial_objective = candidate_objective(
                opt.mode, bit, view, exact, approximate_without_bit,
                v1, v2, &trial_type);
            const double progress = opt.anneal_steps > 1
                ? static_cast<double>(step) / static_cast<double>(opt.anneal_steps - 1)
                : 1.0;
            const double temperature = std::max(
                1.0e-9, static_cast<double>(opt.anneal_temperature) *
                std::pow(1.0e-4, progress));
            const bool accept = trial_objective <= objective ||
                probability(random) < std::exp((objective - trial_objective) / temperature);
            if (accept) {
                objective = trial_objective;
                type = std::move(trial_type);
                if (objective < best.objective) {
                    best.v1 = v1;
                    best.v2 = v2;
                    best.type = type;
                    best.objective = objective;
                    best.steps = -opt.refine_restarts - opt.swap_restarts - restart - 1;
                }
            } else {
                pattern[row] ^= 1U;
            }
        }
    }
    std::uniform_int_distribution<int> column_index(0, view.cols - 1);
    for (int restart = 0; restart < opt.row_sa_restarts; ++restart) {
        std::vector<uint8_t> pattern(view.cols);
        if (restart == 0 && !best.type.empty()) {
            pattern = best.type;
        } else {
            for (int col = 0; col < view.cols; ++col)
                pattern[col] = static_cast<uint8_t>(random_bit(random));
        }
        std::vector<uint8_t> v1, v2;
        double row_objective = solution_from_row_pattern(
            opt.mode, bit, view, exact, approximate_without_bit,
            pattern, v1, v2);
        for (int step = 0; step < opt.row_sa_steps; ++step) {
            const int col = column_index(random);
            pattern[col] ^= 1U;
            std::vector<uint8_t> trial_v1, trial_v2;
            const double trial_row_objective = solution_from_row_pattern(
                opt.mode, bit, view, exact, approximate_without_bit,
                pattern, trial_v1, trial_v2);
            const double progress = opt.row_sa_steps > 1
                ? static_cast<double>(step) / static_cast<double>(opt.row_sa_steps - 1)
                : 1.0;
            const double temperature = std::max(
                1.0e-9, static_cast<double>(opt.row_sa_temperature) *
                std::pow(1.0e-4, progress));
            const bool accept = trial_row_objective <= row_objective ||
                probability(random) <
                    std::exp((row_objective - trial_row_objective) / temperature);
            if (accept) {
                row_objective = trial_row_objective;
                v1 = std::move(trial_v1);
                v2 = std::move(trial_v2);
                std::vector<uint8_t> optimal_type;
                const double objective = candidate_objective(
                    opt.mode, bit, view, exact, approximate_without_bit,
                    v1, v2, &optimal_type);
                if (objective < best.objective) {
                    best.v1 = v1;
                    best.v2 = v2;
                    best.type = std::move(optimal_type);
                    best.objective = objective;
                    best.steps = -opt.refine_restarts - opt.swap_restarts -
                                 opt.anneal_restarts - restart - 1;
                }
            } else {
                pattern[col] ^= 1U;
            }
        }
    }
    return best;
}

#ifdef HAVE_GUROBI
GRBEnv & gurobi_environment() {
    static GRBEnv environment = [] {
        GRBEnv value(true);
        value.set(GRB_IntParam_OutputFlag, 0);
        value.start();
        return value;
    }();
    return environment;
}

ColumnSolution solve_gurobi(const Options & opt, int bit, const MatrixView & view,
                            const std::vector<uint64_t> & exact,
                            const std::vector<uint64_t> & approximate_without_bit) {
    const Eigen::MatrixXf coefficient =
        build_coefficients(opt.mode, bit, view, exact, approximate_without_bit);
    GRBModel model(gurobi_environment());
    model.set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_DoubleParam_TimeLimit, opt.oracle_time);
    model.set(GRB_DoubleParam_MIPGap, 0.0);
    std::vector<GRBVar> v1(view.rows), v2(view.rows), type(view.cols);
    std::vector<GRBVar> product1(view.rows * view.cols);
    std::vector<GRBVar> product2(view.rows * view.cols);
    for (int row = 0; row < view.rows; ++row) {
        v1[row] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        v2[row] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    }
    for (int col = 0; col < view.cols; ++col)
        type[col] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    GRBLinExpr objective = 0.0;
    for (int row = 0; row < view.rows; ++row) {
        for (int col = 0; col < view.cols; ++col) {
            const int index = row * view.cols + col;
            product1[index] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            product2[index] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            model.addConstr(product1[index] <= v1[row]);
            model.addConstr(product1[index] <= type[col]);
            model.addConstr(product1[index] >= v1[row] + type[col] - 1.0);
            model.addConstr(product2[index] <= v2[row]);
            model.addConstr(product2[index] <= type[col]);
            model.addConstr(product2[index] >= v2[row] + type[col] - 1.0);
            objective += static_cast<double>(coefficient(row, col)) *
                (v1[row] - product1[index] + product2[index]);
        }
    }
    model.setObjective(objective, GRB_MINIMIZE);
    model.optimize();
    if (model.get(GRB_IntAttr_SolCount) == 0)
        throw std::runtime_error("Gurobi core oracle found no feasible solution");
    ColumnSolution solution;
    solution.v1.resize(view.rows);
    solution.v2.resize(view.rows);
    solution.type.resize(view.cols);
    for (int row = 0; row < view.rows; ++row) {
        solution.v1[row] = static_cast<uint8_t>(v1[row].get(GRB_DoubleAttr_X) > 0.5);
        solution.v2[row] = static_cast<uint8_t>(v2[row].get(GRB_DoubleAttr_X) > 0.5);
    }
    solution.objective = candidate_objective(
        opt.mode, bit, view, exact, approximate_without_bit,
        solution.v1, solution.v2, &solution.type);
    return solution;
}
#endif

void apply_solution(int bit, const MatrixView & view, const ColumnSolution & solution,
                    std::vector<uint64_t> & approximate) {
    const uint64_t mask = uint64_t{1} << bit;
    for (int address = 0; address < static_cast<int>(approximate.size()); ++address) {
        const int row = view.row_of_address[address];
        const int col = view.col_of_address[address];
        const bool value = solution.type[col] ? solution.v2[row] : solution.v1[row];
        approximate[address] &= ~mask;
        if (value) approximate[address] |= mask;
    }
}

std::vector<uint64_t> run_decomposition(
    const Options & opt, const std::vector<uint64_t> & exact, int n_inputs,
    int n_outputs, const std::vector<Partition> & partitions) {
    std::vector<uint64_t> approximate = exact;
    std::mt19937 random(opt.seed);
    for (int round = 0; round < opt.rounds; ++round) {
        for (int bit = n_outputs - 1; bit >= 0; --bit) {
            const uint64_t mask = uint64_t{1} << bit;
            std::vector<uint64_t> approximate_without_bit = approximate;
            for (auto & value : approximate_without_bit) value &= ~mask;

            ColumnSolution best;
            int best_partition = -1;
            MatrixView best_view;
            double best_secondary = std::numeric_limits<double>::infinity();
            struct RankedCandidate {
                ColumnSolution solution;
                MatrixView view;
                double secondary = 0.0;
            };
            std::vector<RankedCandidate> ranked_candidates;
            for (int partition_id = 0;
                 partition_id < static_cast<int>(partitions.size()); ++partition_id) {
                const MatrixView view = make_matrix_view(
                    n_inputs, partitions[partition_id], opt.row_set);
                for (int restart = 0; restart < opt.restarts; ++restart) {
                    ColumnSolution candidate;
#ifdef HAVE_GUROBI
                    if (opt.core_solver == "gurobi")
                        candidate = solve_gurobi(
                            opt, bit, view, exact, approximate_without_bit);
                    else
#endif
                        candidate = solve_bsb(
                            opt, bit, view, exact, approximate_without_bit, random);
                    double secondary = 0.0;
                    if (opt.mode == "separate") {
                        std::vector<uint64_t> trial_approximate = approximate;
                        apply_solution(bit, view, candidate, trial_approximate);
                        secondary = mean_error_distance(exact, trial_approximate);
                        if (opt.separate_slack > 0.0F)
                            ranked_candidates.push_back(
                                {candidate, view, secondary});
                    }
                    const bool primary_better =
                        candidate.objective + 1.0e-12 < best.objective;
                    const bool primary_tied =
                        std::abs(candidate.objective - best.objective) <= 1.0e-12;
                    if (primary_better ||
                        (opt.mode == "separate" && primary_tied &&
                         secondary < best_secondary)) {
                        best = std::move(candidate);
                        best_partition = partition_id;
                        best_view = view;
                        best_secondary = secondary;
                    }
                }
            }
            if (opt.mode == "separate" && opt.separate_slack > 0.0F &&
                !ranked_candidates.empty()) {
                const double minimum_primary = std::min_element(
                    ranked_candidates.begin(), ranked_candidates.end(),
                    [](const RankedCandidate & lhs, const RankedCandidate & rhs) {
                        return lhs.solution.objective < rhs.solution.objective;
                    })->solution.objective;
                best_secondary = std::numeric_limits<double>::infinity();
                for (RankedCandidate & ranked : ranked_candidates) {
                    if (ranked.solution.objective <=
                            minimum_primary + opt.separate_slack + 1.0e-12 &&
                        ranked.secondary < best_secondary) {
                        best = std::move(ranked.solution);
                        best_view = std::move(ranked.view);
                        best_secondary = ranked.secondary;
                        best_partition = 0;
                    }
                }
            }
            if (best_partition < 0)
                throw std::runtime_error("solver produced no solution");
            apply_solution(bit, best_view, best, approximate);
            if (!opt.quiet) {
                std::cout << "round=" << round << " bit=" << bit
                          << " partition=" << best_partition
                          << " core_objective=" << std::setprecision(10)
                          << best.objective
                          << " med=" << mean_error_distance(exact, approximate)
                          << " best_step=" << best.steps << '\n';
            }
        }
    }
    return approximate;
}

}  // namespace

int main(int argc, char ** argv) {
    try {
        const Options opt = make_options(argc, argv);
        const std::vector<uint64_t> exact = load_table(opt.input);
        const int n_inputs = exact_log2(exact.size());
        const int n_outputs = output_bits(exact);
        if (opt.n_free <= 0 || opt.n_free >= n_inputs)
            throw std::invalid_argument("--n-free must be in [1, n_inputs-1]");

        std::mt19937 random(opt.seed);
        const std::vector<Partition> partitions =
            make_partitions(n_inputs, opt.n_free, opt.max_partitions, random);
        std::cout << "mode=" << opt.mode << " inputs=" << n_inputs
                  << " outputs=" << n_outputs << " free=" << opt.n_free
                  << " bound=" << (n_inputs - opt.n_free)
                  << " row_set=" << opt.row_set
                  << " partitions=" << partitions.size() << " seed=" << opt.seed << '\n';

        const auto start = std::chrono::steady_clock::now();
        std::vector<uint64_t> approximate;
        if (!opt.paper_plus) {
            approximate = run_decomposition(
                opt, exact, n_inputs, n_outputs, partitions);
        } else {
            std::vector<Options> variants;
            if (opt.mode == "separate") {
                variants.assign(2, opt);
                variants[0].field_schedule = "sqrt";
                variants[0].field_scale = 0.7F;
                variants[1].field_schedule = "coupled-linear";
                variants[1].field_scale = 1.0F;
                for (Options & variant : variants)
                    variant.normalization = "full";
            } else {
                variants.assign(3, opt);
                variants[0].field_schedule = "coupled-linear";
                variants[0].field_scale = 0.329F;
                variants[1].field_schedule = "linear";
                variants[1].field_scale = 1.15F;
                variants[2].field_schedule = "coupled";
                variants[2].field_scale = 1.0F;
                variants[0].normalization = "full";
                variants[1].normalization = "full";
                variants[2].normalization = "pairwise";
            }
            double best_med = std::numeric_limits<double>::infinity();
            std::vector<std::future<std::vector<uint64_t>>> futures;
            futures.reserve(variants.size());
            for (Options & variant : variants) {
                variant.paper_plus = false;
                variant.quiet = true;
                futures.push_back(std::async(
                    std::launch::async,
                    [variant, &exact, n_inputs, n_outputs, &partitions] {
                        return run_decomposition(
                            variant, exact, n_inputs, n_outputs, partitions);
                    }));
            }
            for (auto & future : futures) {
                std::vector<uint64_t> trial = future.get();
                const double med = mean_error_distance(exact, trial);
                if (med < best_med) {
                    best_med = med;
                    approximate = std::move(trial);
                }
            }
        }
        const auto stop = std::chrono::steady_clock::now();
        const double seconds = std::chrono::duration<double>(stop - start).count();
        std::cout << "RESULT mode=" << opt.mode
                  << " med=" << std::setprecision(10)
                  << mean_error_distance(exact, approximate)
                  << " time_s=" << seconds << '\n';
        return 0;
    } catch (const std::exception & error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
