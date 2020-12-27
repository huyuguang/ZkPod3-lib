#pragma once

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

    std::vector<G1> com_w(input.s);
    std::vector<Fr> com_w_r(input.s);

    // public input
    G1 sum_g = pc::ComputeSigmaG(0, input.D * 3);
    for (size_t i = 0; i < input.D; ++i) {
      assert(input.w[i].size() == input.D * 3);
      for (size_t j = 0; j < input.D * 3; ++j) {
        assert(input.w[i][j] == i + 1);
      }
      com_w_r[i] = FrZero();
      com_w[i] = sum_g * (i + 1);
    }

    {
      auto parallel_f = [&com_w_r, &com_w, &input](int64_t i) {
        com_w_r[i] = FrRand();
        com_w[i] = pc::ComputeCom(input.get_g, input.w[i], com_w_r[i]);
      };
      parallel::For<int64_t>(input.D, input.s, parallel_f);
    }

    // prove permutation
    typename R1cs::ProveInput r1cs_input(*input.r1cs_info, "sudoku",
                                         std::move(input.w), com_w, com_w_r,
                                         input.get_g);
    R1cs::Prove(proof.r1cs_proof, seed, std::move(r1cs_input));
    proof.com_w = std::move(com_w);

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

    AdaptProve(seed, std::move(adapt_items), proof.adapt_proof);
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

    typename ParallelR1cs<Policy>::VerifyInput pr_input(
        input.D * 3, *input.r1cs_info, "sudoku", proof.com_w, input.public_w,
        input.get_g);
    if (!ParallelR1cs<Policy>::Verify(proof.r1cs_proof, seed, pr_input))
      return false;

    AdaptUpdateSeed(seed, proof.com_x, proof.com_w);
    // std::cout << "verify seed: " <<  misc::HexToStr(seed) << "\n";

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

    if (!AdaptVerify(seed, std::move(adapt_items), proof.adapt_proof)) {
      assert(false);
      return false;
    }
    return true;
  }

  static bool Test() {
    const static size_t d = 3;
    //const static size_t D = d * d;
    std::vector<Fr> x{
      7, 5, 2, 8, 3, 9, 6, 1, 4,
      3, 9, 1, 4, 5, 6, 2, 8, 7,
      6, 8, 4, 1, 7, 2, 9, 5, 3,
      2, 1, 7, 9, 6, 4, 5, 3, 8,
      5, 4, 9, 3, 8, 7, 1, 6, 2,
      8, 3, 6, 5, 2, 1, 4, 7, 9,
      4, 7, 3, 2, 1, 5, 8, 9, 6,
      9, 6, 5, 7, 4, 8, 3, 2, 1,
      1, 2, 8, 6, 9, 3, 7, 4, 5};

    Fr com_x_r = FrRand();
    G1 com_x = pc::ComputeCom(x, com_x_r);
    std::vector<size_t> open_positions{0, 3, 21};
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