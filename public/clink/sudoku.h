#pragma once

#include <algorithm>

#include "./adapt.h"
#include "./details.h"
#include "./parallel_r1cs.h"
#include "circuit/sudoku_gadget.h"

namespace clink {
struct Sudoku {
  // using Policy = groth09::OrdinaryPolicy;
  using Policy = groth09::SuccinctPolicy;
  using R1cs = typename clink::ParallelR1cs<Policy>;
  using HyraxA = typename Policy::HyraxA;
  using Sec51 = typename Policy::Sec51;
  using Sec53 = typename Policy::Sec53;
  using Sec43 = typename Policy::Sec43;

  struct Proof {
    clink::ParallelR1cs<Policy>::Proof r1cs_proof;
    hyrax::A4::Proof adapt_proof;
    G1 com_x;
    std::vector<G1> com_w;
    bool operator==(Proof const& b) const {
      return r1cs_proof == b.r1cs_proof && com_x == b.com_x && com_w == b.com_w;
    }

    bool operator!=(Proof const& b) const { return !(*this == b); }

    template <typename Ar>
    void serialize(Ar& ar) const {
      ar& YAS_OBJECT_NVP("sudoku.p", ("rp", r1cs_proof), ("rp", adapt_proof),
                         ("cx", com_x), ("cw", com_w));
    }
    template <typename Ar>
    void serialize(Ar& ar) {
      ar& YAS_OBJECT_NVP("sudoku.p", ("rp", r1cs_proof), ("rp", adapt_proof),
                         ("cx", com_x), ("cw", com_w));
    }
  };

  struct ProveInput {
    size_t d;
    size_t D;
    std::vector<Fr> const& x;
    G1 const& com_x;
    Fr const& com_x_r;
    std::vector<size_t> const& open_positions;
    GetRefG1 const& get_g;

    std::vector<Fr> open_values;
    std::unique_ptr<R1csInfo> r1cs_info;
    int64_t s;
    std::vector<std::vector<Fr>> mutable w;
    std::vector<std::vector<Fr>> row_col_cells;
    std::vector<std::vector<size_t>> row_col_cells_pos;

    Fr const& get_x(size_t i, size_t j) const { return x[i * D + j]; }
    ProveInput(std::vector<Fr> const& x, G1 const& com_x, Fr const& com_x_r,
               std::vector<size_t> const& open_positions, GetRefG1 const& get_g)
        : x(x),
          com_x(com_x),
          com_x_r(com_x_r),
          open_positions(open_positions),
          get_g(get_g) {
      D = (size_t)(std::sqrt(x.size()));
      if (D * D != x.size()) throw std::invalid_argument("invalid dimension");
      d = (size_t)std::sqrt(D);
      if (d * d != D) throw std::invalid_argument("invalid dimension");

      libsnark::protoboard<Fr> pb;
      circuit::SudokuGadget gadget(pb, D);

      int64_t const primary_input_size = D;
      pb.set_input_sizes(primary_input_size);
      r1cs_info.reset(new R1csInfo(pb));
      s = r1cs_info->num_variables;
      w.resize(s);
      for (auto& i : w) i.resize(D * 3);

      row_col_cells.resize(D * 3);

      // every row
      for (size_t i = 0; i < D; ++i) {
        std::vector<Fr>& row = row_col_cells[i];
        row.resize(D);
        for (size_t j = 0; j < D; ++j) {
          row[j] = get_x(i, j);
        }
      }
      // every col
      for (size_t i = 0; i < D; ++i) {
        std::vector<Fr>& col = row_col_cells[i + D];
        col.resize(D);
        for (size_t j = 0; j < D; ++j) {
          col[j] = get_x(j, i);
        }
      }
      // every cell
      for (size_t i = 0; i < D; ++i) {
        std::vector<Fr>& cell = row_col_cells[i + D * 2];
        cell.resize(D);
        for (size_t j = 0; j < D; ++j) {
          size_t ii = (i / d) * d + j / d;
          size_t jj = (i % d) * d + j % d;
          cell[j] = get_x(ii, jj);
        }
      }

      row_col_cells_pos.resize(D * 3);
      // every row
      for (size_t i = 0; i < D; ++i) {
        auto& pos = row_col_cells_pos[i];
        pos.resize(D);
        for (size_t j = 0; j < D; ++j) {
          pos[j] = i * D + j;
        }
      }
      // every col
      for (size_t i = 0; i < D; ++i) {
        auto& pos = row_col_cells_pos[i + D];
        pos.resize(D);
        for (size_t j = 0; j < D; ++j) {
          pos[j] = j * D + i;
        }
      }
      // every cell
      for (size_t i = 0; i < D; ++i) {
        auto& pos = row_col_cells_pos[i + D * 2];
        pos.resize(D);
        for (size_t j = 0; j < D; ++j) {
          size_t ii = (i / d) * d + j / d;
          size_t jj = (i % d) * d + j % d;
          pos[j] = ii * D + jj;
        }
      }

      for (size_t j = 0; j < D * 3; ++j) {
        gadget.Assign(row_col_cells[j]);
        assert(pb.is_satisfied());
        auto v = pb.full_variable_assignment();
        for (int64_t i = 0; i < s; ++i) {
          w[i][j] = v[i];
        }
      }

      open_values.resize(open_positions.size());
      for (size_t i = 0; i < open_positions.size(); ++i) {
        open_values[i] = x[open_positions[i]];
      }
    }
  };

  static void AdaptUpdateSeed(h256_t& seed, G1 const& com_x,
                              std::vector<G1> const& com_w) {
    CryptoPP::Keccak_256 hash;
    HashUpdate(hash, seed);
    HashUpdate(hash, com_x);
    for (auto const& i : com_w) {
      HashUpdate(hash, i);
    }
    hash.Final(seed.data());
  }

  static void AdaptComputeFst(h256_t const& seed, std::string const& tag,
                              std::vector<Fr>& e) {
    std::string salt = "sudoku adapt " + std::to_string(e.size()) + tag;
    ComputeFst(seed, salt, e);
  }

  static void Prove(Proof& proof, h256_t seed, ProveInput const& input) {
    Tick tick(__FN__);
    proof.com_x = input.com_x;

    proof.com_w.resize(input.s);
    std::vector<Fr> com_w_r(input.s);

    // public input
    G1 sum_g = pc::ComputeSigmaG(0, input.D * 3);
    for (size_t i = 0; i < input.D; ++i) {
      assert(input.w[i].size() == input.D * 3);
      for (size_t j = 0; j < input.D * 3; ++j) {
        assert(input.w[i][j] == i + 1);
      }
      com_w_r[i] = FrZero();
      proof.com_w[i] = sum_g * (i + 1);
    }

    {
      auto parallel_f = [&com_w_r, &proof, &input](int64_t i) {
        com_w_r[i] = FrRand();
        proof.com_w[i] = pc::ComputeCom(input.get_g, input.w[i], com_w_r[i]);
      };
      parallel::For<int64_t>(input.D, input.s, parallel_f);
    }

    // prove adapt
    std::vector<AdaptProveItem> adapt_items(input.D + 1);
    std::vector<std::vector<Fr>> e(input.D + 1);
    for (auto& i : e) i.reserve(input.D * input.D);
    std::vector<std::vector<Fr>> ee(input.D + 1);
    for (auto& i : ee) i.resize(input.D * input.D);

    // consistent with open values
    AdaptUpdateSeed(seed, input.com_x, proof.com_w);
    // std::cout << "prove seed: " <<  misc::HexToStr(seed) << "\n";

    e[0].resize(input.open_positions.size());
    AdaptComputeFst(seed, "open", e[0]);
    for (auto& i : ee[0]) i = FrZero();
    for (size_t i = 0; i < input.open_positions.size(); ++i) {
      ee[0][input.open_positions[i]] = e[0][i];
    }
    Fr z = InnerProduct(e[0], input.open_values);

    AdaptProveItem& adapt_open = adapt_items[0];
    adapt_open.Init(1, "open", z);
    adapt_open.x[0] = input.x;
    adapt_open.a[0] = std::move(ee[0]);
    adapt_open.cx[0] = input.com_x;
    adapt_open.rx[0] = input.com_x_r;
    // assert(adapt_open.CheckData());

    // consistent with com_w[N~2N)
    for (size_t i = 0; i < input.D; ++i) {
      auto& this_e = e[i + 1];
      auto& this_ee = ee[i + 1];
      this_e.resize(input.D * 3);
      AdaptComputeFst(seed, "open", this_e);
      std::vector<Fr> row(input.D * 3);
      for (size_t j = 0; j < input.D * 3; ++j) {
        row[j] = input.row_col_cells[j][i];
      }
      for (auto& item : this_ee) item = FrZero();
      for (size_t j = 0; j < input.D * 3; ++j) {
        auto pos = input.row_col_cells_pos[j][i];
        this_ee[pos] += -this_e[j];
      }
      AdaptProveItem& adapt_cell = adapt_items[i + 1];
      adapt_cell.Init(2, "cell" + std::to_string(i), FrZero());
      adapt_cell.x[0] = input.x;
      adapt_cell.a[0] = std::move(this_ee);
      adapt_cell.cx[0] = input.com_x;
      adapt_cell.rx[0] = input.com_x_r;
      adapt_cell.x[1] = std::move(row);
      adapt_cell.a[1] = std::move(this_e);
      adapt_cell.cx[1] = proof.com_w[input.D + i];
      adapt_cell.rx[1] = com_w_r[input.D + i];
      // assert(adapt_cell.CheckData());
    }

    std::array<parallel::VoidTask, 2> tasks;

    tasks[0] = [&seed,&adapt_items, &proof]() {
      AdaptProve(seed, std::move(adapt_items), proof.adapt_proof);
    };

    tasks[1] = [&seed, &proof, &input, &com_w_r]() {
      // prove permutation
      typename R1cs::ProveInput r1cs_input(*input.r1cs_info, "sudoku",
                                           std::move(input.w), proof.com_w,
                                           com_w_r, input.get_g);
      R1cs::Prove(proof.r1cs_proof, seed, std::move(r1cs_input));
    };

    parallel::Invoke(tasks);
  }

  struct VerifyInput {
    size_t d;
    size_t D;
    std::vector<size_t> const& open_positions;
    std::vector<Fr> const& open_values;
    GetRefG1 const& get_g;

    std::unique_ptr<R1csInfo> r1cs_info;
    int64_t m;
    int64_t s;
    std::vector<std::vector<Fr>> public_w;
    std::vector<std::vector<size_t>> row_col_cells_pos;

    VerifyInput(size_t d, std::vector<size_t> const& open_positions,
                std::vector<Fr> const& open_values, GetRefG1 const& get_g)
        : d(d),
          D(d * d),
          open_positions(open_positions),
          open_values(open_values),
          get_g(get_g) {
      libsnark::protoboard<Fr> pb;
      circuit::SudokuGadget gadget(pb, D);
      int64_t const primary_input_size = D;
      pb.set_input_sizes(primary_input_size);
      r1cs_info.reset(new R1csInfo(pb));
      m = r1cs_info->num_constraints;
      s = r1cs_info->num_variables;

      // public input
      public_w.resize(D);
      for (size_t i = 0; i < D; ++i) {
        auto& row = public_w[i];
        row.resize(D * 3);
        for (size_t j = 0; j < D * 3; ++j) {
          row[j] = i + 1;
        }
      }

      row_col_cells_pos.resize(D * 3);
      // every row
      for (size_t i = 0; i < D; ++i) {
        auto& pos = row_col_cells_pos[i];
        pos.resize(D);
        for (size_t j = 0; j < D; ++j) {
          pos[j] = i * D + j;
        }
      }
      // every col
      for (size_t i = 0; i < D; ++i) {
        auto& pos = row_col_cells_pos[i + D];
        pos.resize(D);
        for (size_t j = 0; j < D; ++j) {
          pos[j] = j * D + i;
        }
      }
      // every cell
      for (size_t i = 0; i < D; ++i) {
        auto& pos = row_col_cells_pos[i + D * 2];
        pos.resize(D);
        for (size_t j = 0; j < D; ++j) {
          size_t ii = (i / d) * d + j / d;
          size_t jj = (i % d) * d + j % d;
          pos[j] = ii * D + jj;
        }
      }
    }
  };

  static bool Verify(Proof const& proof, h256_t seed,
                     VerifyInput const& input) {
    Tick tick(__FN__);
    if ((int64_t)proof.com_w.size() != input.s) {
      assert(false);
      return false;
    }

    // public input
    G1 sum_g = pc::ComputeSigmaG(0, input.D * 3);
    for (size_t i = 0; i < input.D; ++i) {
      if (proof.com_w[i] != sum_g * (i + 1)) {
        assert(false);
        return false;
      }
    }

    AdaptUpdateSeed(seed, proof.com_x, proof.com_w);
    
    // verify adapt
    std::vector<AdaptVerifyItem> adapt_items(input.D + 1);
    std::vector<std::vector<Fr>> e(input.D + 1);
    for (auto& i : e) i.reserve(input.D * input.D);
    std::vector<std::vector<Fr>> ee(input.D + 1);
    for (auto& i : ee) i.resize(input.D * input.D);

    // consistent with open values
    e[0].resize(input.open_positions.size());
    AdaptComputeFst(seed, "open", e[0]);
    for (auto& i : ee[0]) i = FrZero();
    for (size_t i = 0; i < input.open_positions.size(); ++i) {
      ee[0][input.open_positions[i]] = e[0][i];
    }
    Fr z = InnerProduct(e[0], input.open_values);

    AdaptVerifyItem& adapt_open = adapt_items[0];
    adapt_open.Init(1, "open", z);
    adapt_open.a[0] = std::move(ee[0]);
    adapt_open.cx[0] = proof.com_x;

    // consistent with com_w[N~2N)
    for (size_t i = 0; i < input.D; ++i) {
      auto& this_e = e[i + 1];
      auto& this_ee = ee[i + 1];
      this_e.resize(input.D * 3);
      AdaptComputeFst(seed, "open", this_e);
      for (auto& item : this_ee) item = FrZero();
      for (size_t j = 0; j < input.D * 3; ++j) {
        auto pos = input.row_col_cells_pos[j][i];
        this_ee[pos] += -this_e[j];
      }
      AdaptVerifyItem& adapt_cell = adapt_items[i + 1];
      adapt_cell.Init(2, "cell" + std::to_string(i), FrZero());
      adapt_cell.a[0] = std::move(this_ee);
      adapt_cell.cx[0] = proof.com_x;
      adapt_cell.a[1] = std::move(this_e);
      adapt_cell.cx[1] = proof.com_w[input.D + i];
    }

    std::array<std::atomic<bool>, 2> rets;
    std::array<parallel::VoidTask, 2> tasks;
    tasks[0] = [&input, &proof, &seed, &rets]() {
      typename ParallelR1cs<Policy>::VerifyInput pr_input(
          input.D * 3, *input.r1cs_info, "sudoku", proof.com_w, input.public_w,
          input.get_g);
      rets[0] = ParallelR1cs<Policy>::Verify(proof.r1cs_proof, seed, pr_input);
    };

    tasks[1] = [&adapt_items, &seed, &proof,&rets]() {
      rets[1] = AdaptVerify(seed, std::move(adapt_items), proof.adapt_proof);
    };

    parallel::Invoke(tasks);

    if (!rets[0] || !rets[1]) {
      assert(false);
      return false;
    }
    return true;
  }

  static std::vector<Fr> GeneratePuzzle(size_t d) {    
    size_t D = d * d;
    std::vector<std::vector<size_t>> matrix;
    matrix.reserve(D);

    std::vector<size_t> base(D);
    for (size_t i = 0; i < D; ++i) {
      base[i] = i + 1;
    }
    std::vector<size_t> rnd(d);
    for (size_t i = 0; i < d; ++i) {
      rnd[i] = i;
    }

    std::random_device rng;
    std::mt19937 urng(rng());
    auto buf = base;
    std::shuffle(buf.begin(), buf.end(), urng);

    for (size_t i = 0; i < d; ++i) {
      std::rotate(buf.begin(), buf.begin() + 1, buf.end());
      matrix.push_back(buf);
      for (size_t j = 1; j < d; ++j) {
        std::rotate(buf.begin(), buf.begin() + d, buf.end());
        matrix.push_back(buf);
      }
    }

    // shuffle
    for (size_t k = 0; k < d; ++k) {
      auto begin = matrix.begin() + k * d;
      auto end = begin + d;
      std::shuffle(begin, end, urng);
    }

    // rotate 90
    std::vector<std::vector<size_t>> dup(D);
    for (size_t i = 0; i < D; ++i) {
      dup[i].resize(D);
      for (size_t j = 0; j < D; ++j) {
        dup[i][j] = matrix[j][i];
      }
    }
    matrix = std::move(dup);

    // shuffle
    for (size_t k = 0; k < d; ++k) {
      auto begin = matrix.begin() + k * d;
      auto end = begin + d;
      std::shuffle(begin, end, urng);
    }

    // matrix to vector
    std::vector<size_t> ret(D * D);
    for (size_t i = 0; i < D; ++i) {
      for (size_t j = 0; j < D; ++j) {
        ret[i * D + j] = matrix[i][j];
      }
    }

    // check row
    for (size_t i = 0; i < D; ++i) {
      buf.resize(0);
      for (size_t j = 0; j < D; ++j) {
        buf.push_back(ret[i * D + j]);
      }
      std::sort(buf.begin(), buf.end());
      assert(buf == base);
    }

    // check col
    for (size_t j = 0; j < D; ++j) {
      buf.resize(0);
      for (size_t i = 0; i < D; ++i) {
        buf.push_back(ret[i * D + j]);
      }
      std::sort(buf.begin(), buf.end());
      assert(buf == base);
    }

    // check cell
    for (size_t i = 0; i < D; ++i) {
      buf.resize(0);
      for (size_t j = 0; j < D; ++j) {
        size_t ii = (i / d) * d + j / d;
        size_t jj = (i % d) * d + j % d;
        buf.push_back(ret[ii * D + jj]);
      }
      std::sort(buf.begin(), buf.end());
      assert(buf == base);
    }

    for (size_t i = 0; i < D; ++i) {
      for (size_t j = 0; j < D; ++j) {
        std::cout << std::right << std::setw(4) << std::setfill(' ')
                  << ret[i * D + j];
      }
      std::cout << "\n";
    }
    std::cout << "\n";

    std::vector<Fr> fr_ret(ret.size());
    for (size_t i = 0; i < ret.size(); ++i) {
      fr_ret[i] = ret[i];
    }
    return fr_ret;
  }

  static bool Test(size_t d) {
    size_t D = d * d;
    auto x = GeneratePuzzle(d);
    assert(x.size() == D * D);

    Fr com_x_r = FrRand();
    G1 com_x = pc::ComputeCom(x, com_x_r);

    std::vector<size_t> open_positions(D * D / 3);
    for (auto& i : open_positions) i = rand() % (D * D);
    std::sort(open_positions.begin(), open_positions.end());
    open_positions.erase(
        std::unique(open_positions.begin(), open_positions.end()),
        open_positions.end());

    std::vector<Fr> open_values(open_positions.size());
    for (size_t i = 0; i < open_positions.size(); ++i) {
      open_values[i] = x[open_positions[i]];
    }

    ProveInput prove_input(x, com_x, com_x_r, open_positions, pc::kGetRefG1);
    h256_t seed = misc::RandH256();
    Proof proof;
    Prove(proof, seed, prove_input);

    VerifyInput verify_input(d, open_positions, open_values, pc::kGetRefG1);
    bool success = Verify(proof, seed, verify_input);
    std::cout << __FILE__ << " " << __FN__ << ": " << success << "\n\n\n\n\n\n";
    return success;
  }
};
}  // namespace clink