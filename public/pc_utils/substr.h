#pragma once

#include <libff/algebra/fields/bigint.hpp>
#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/gadgetlib1/pb_variable.hpp>

#include "./details.h"
#include "./parallel_r1cs.h"
#include "./types.h"

// x: vector<Fr>, size = n
// k: Fr
// y: vector<Fr>, size = n, {0,1}
// open com(gx,x), k, com(gx,y)
// prove y[i] = IsSubStr(x[i],k)? 1:0

namespace pc_utils {

template <typename Policy>
struct Substr {
  using Sec53 = typename Policy::Sec53;
  using HyraxA = typename Policy::HyraxA;
  using R1cs = typename ParallelR1cs<Policy>;

  // if (is_any(x,0)) ret = 1; else ret = 0;
  class HasZeroGadget : public libsnark::gadget<Fr> {
   public:
    HasZeroGadget(libsnark::protoboard<Fr>& pb,
                  std::vector<libsnark::pb_linear_combination<Fr>> const& x)
        : libsnark::gadget<Fr>(pb, "HasZeroGadget"), x_(x) {
      assert(!x_.empty());
      y_.resize(x_.size() - 1);
      for (auto& i : y_) i.allocate(pb, " y");
      inv_.allocate(pb, " inv");
      ret_.allocate(pb, " ret");
    }

    libsnark::pb_variable<Fr> const& ret() const { return ret_; }

    // ret = (x0*x1...)? 0 : 1
    void generate_r1cs_constraints() {
      libsnark::pb_linear_combination<Fr> pie_of_x;
      if (x_.size() == 1) {
        pie_of_x = x_.back();
      } else {
        // x[0] * x[1] == y_[0]
        this->pb.add_r1cs_constraint(
            libsnark::r1cs_constraint<Fr>(x_[0], x_[1], y_[0]),
            "y[0] = x[0] * x[1]");
        for (size_t i = 1; i < x_.size() - 1; ++i) {
          // y[i-1] * x[i+1] == y[i]
          this->pb.add_r1cs_constraint(
              libsnark::r1cs_constraint<Fr>(y_[i - 1], x_[i + 1], y_[i]),
              "y[i] = y[i-1] * x[i+1]");
        }

        pie_of_x = y_.back();
      }

      // ret * pie_of_x == 0
      this->pb.add_r1cs_constraint(
          libsnark::r1cs_constraint<Fr>(ret_, pie_of_x, FrZero()),
          "ret_ * pie_of_x == 0");
      // inv * pie_of_x == 1 - ret
      this->pb.add_r1cs_constraint(
          libsnark::r1cs_constraint<Fr>(inv_, pie_of_x, -ret_ + FrOne()),
          "inv * pie_of_x == 1 - ret");
    }

    void generate_r1cs_witness() {
      Fr pie_of_x;
      if (x_.size() == 1) {
        // y_ is empty
        pie_of_x = this->pb.lc_val(x_[0]);
      } else {
        // y_[0] = x_[0] * x_[1];
        this->pb.val(y_[0]) = this->pb.lc_val(x_[0]) * this->pb.lc_val(x_[1]);
        for (size_t i = 1; i < x_.size() - 1; ++i) {
          // y[i] = y[i-1] * x[i+1]
          this->pb.val(y_[i]) =
              this->pb.val(y_[i - 1]) * this->pb.lc_val(x_[i + 1]);
        }

        pie_of_x = this->pb.val(y_.back());
      }
      this->pb.val(inv_) = pie_of_x == FrZero() ? FrZero() : FrInv(pie_of_x);
      this->pb.val(ret_) = pie_of_x == FrZero() ? FrOne() : FrZero();
    }

   private:
    std::vector<libsnark::pb_linear_combination<Fr>> const& x_;
    std::vector<libsnark::pb_variable<Fr>> y_;  // (x1*x2...)
    libsnark::pb_variable<Fr> inv_;             // y.back()? 1/y.back():0
    libsnark::pb_variable<Fr> ret_;             // y.back() == 0? 1 : 0
  };

  class SubstrGadget : public libsnark::gadget<Fr> {
   public:
    SubstrGadget(libsnark::protoboard<Fr>& pb, std::string const& k)
        : libsnark::gadget<Fr>(pb, "SubstrGadget"), k_(PackStrToFr(k.c_str())) {
      // NOTE: 1 Fr pack 31 char
      x_.allocate(pb, " x");
      dual_.reset(
          new libsnark::dual_variable_gadget<Fr>(this->pb, x_, 31 * 8, "dual"));

      int64_t k_len = k.size();
      sub_duals_.resize(31 - k_len + 1);
      sub_minus_k_.resize(sub_duals_.size());
      for (size_t i = 0; i < sub_duals_.size(); ++i) {
        auto begin = dual_->bits.begin() + i * 8;
        auto end = begin + k_len * 8;
        libsnark::pb_variable_array<Fr> bits(begin, end);
        sub_duals_[i].reset(
            new libsnark::dual_variable_gadget<Fr>(this->pb, bits, "sub_dual"));
        sub_minus_k_[i].assign(this->pb, sub_duals_[i]->packed - k_);
      }

      ret_.reset(new HasZeroGadget(pb, sub_minus_k_));

      generate_r1cs_constraints();
    }

    libsnark::pb_variable<Fr> const& ret() const { return ret_->ret(); }

    void Assign(Fr const& x) {
      this->pb.val(x_) = x;
      generate_r1cs_witness();
      assert(this->pb.is_satisfied());
    }

   private:
    void generate_r1cs_constraints() {
      dual_->generate_r1cs_constraints(true);

      for (auto& i : sub_duals_) {
        i->generate_r1cs_constraints(false);
      }

      ret_->generate_r1cs_constraints();
    }

    void generate_r1cs_witness() {
      dual_->generate_r1cs_witness_from_packed();

      for (auto& i : sub_duals_) {
        i->generate_r1cs_witness_from_bits();
      }

      for (auto& i : sub_minus_k_) {
        i.evaluate(this->pb);
      }

      ret_->generate_r1cs_witness();
    }

   private:
    Fr k_;

    libsnark::pb_variable<Fr> x_;

    // 248 var, 0 or 1
    std::unique_ptr<libsnark::dual_variable_gadget<Fr>> dual_;

    // 31 - k_len + 1 var
    std::vector<std::unique_ptr<libsnark::dual_variable_gadget<Fr>>> sub_duals_;

    std::vector<libsnark::pb_linear_combination<Fr>> sub_minus_k_;

    std::unique_ptr<HasZeroGadget> ret_;
  };

  struct Proof {
    typename R1cs::Proof r1cs_proof;
    std::vector<G1> com_w;
    G1 const& com_x() const { return com_w.front(); }
    G1 const& com_y() const { return com_w.back(); }
    bool operator==(Proof const& b) const {
      return r1cs_proof == b.r1cs_proof && com_w == b.com_w;
    }

    bool operator!=(Proof const& b) const { return !(*this == b); }
  };

  struct ProverInput {
    ProverInput(std::string const& k, std::vector<Fr> const& x, G1 const& com_x,
                Fr const& com_x_r, std::vector<Fr> const& y, G1 const& com_y,
                Fr const& com_y_r, int64_t g_offset)
        : k(k),
          x(x),
          com_x(com_x),
          com_x_r(com_x_r),
          y(y),
          com_y(com_y),
          com_y_r(com_y_r),
          g_offset(g_offset),
          n((int64_t)x.size()),
          gadget(SubstrGadget(pb, k)),
          s((int64_t)pb.num_variables()) {
      int64_t const primary_input_size = 0;
      pb.set_input_sizes(primary_input_size);

      w.resize(s);
      for (auto& i : w) i.resize(n);

      for (int64_t j = 0; j < n; ++j) {
        gadget.Assign(x[j]);
        assert(pb.is_satisfied());
        auto v = pb.full_variable_assignment();
        for (int64_t i = 0; i < s; ++i) {
          w[i][j] = v[i];
        }
        assert(w.front()[j] == x[j]);
        assert(w.back()[j] == y[j]);
      }
    }
    std::string const& k;
    std::vector<Fr> const& x;
    G1 const& com_x;
    Fr const& com_x_r;
    std::vector<Fr> const& y;
    G1 const& com_y;
    Fr const& com_y_r;
    int64_t const g_offset;

    int64_t const n;
    libsnark::protoboard<Fr> pb;
    SubstrGadget gadget;
    int64_t const s;
    std::vector<std::vector<Fr>> mutable w;
  };

  static void Prove(Proof& proof, h256_t seed, ProverInput const& input) {
    Tick tick(__FUNCTION__);
    std::vector<G1> com_w(input.s);
    std::vector<Fr> com_w_r(input.s);

    com_w.front() = input.com_x;
    com_w_r.front() = input.com_x_r;
    com_w.back() = input.com_y;
    com_w_r.back() = input.com_y_r;

    {
      // Tick tick2("compute com(witness)");
      auto parallel_f = [&com_w_r, &com_w, &input](int64_t i) {
        com_w_r[i] = FrRand();
        com_w[i] =
            PcComputeCommitmentG(input.g_offset, input.w[i], com_w_r[i], true);
      };
      parallel::For<int64_t>(1LL, input.s - 1, parallel_f);
    }

    ParallelR1cs<Policy>::ProverInput pr_input(input.pb, std::move(input.w), com_w,
                                        com_w_r, input.g_offset);
    ParallelR1cs<Policy>::Prove(proof.r1cs_proof, seed, std::move(pr_input));
    proof.com_w = std::move(com_w);
  }

  struct VerifierInput {
    VerifierInput(int64_t n, std::string const& k, int64_t g_offset)
        : n(n),
          k(k),
          g_offset(g_offset),
          gadget(SubstrGadget(pb, k)),
          m((int64_t)pb.num_constraints()),
          s((int64_t)pb.num_variables()) {
      int64_t const primary_input_size = 0;
      pb.set_input_sizes(primary_input_size);
    }
    int64_t n;
    std::string const& k;
    int64_t const g_offset;
    libsnark::protoboard<Fr> pb;
    SubstrGadget gadget;
    int64_t const m;
    int64_t const s;
    // since primary_input_size = 0, public_w is empty
    std::vector<std::vector<Fr>> public_w;
  };

  // NOTE: com_x and com_y can get by proof.com_x() and proof.com_y()
  static bool Verify(Proof const& proof, h256_t seed,
                     VerifierInput const& input) {
    Tick tick(__FUNCTION__);
    if ((int64_t)proof.com_w.size() != input.s) {
      assert(false);
      return false;
    }

    // if (proof.r1cs_proof.m() != (int64_t)misc::Pow2UB(input.m)) {
    //  assert(false);
    //  return false;
    //}

    ParallelR1cs<Policy>::VerifierInput pr_input(input.n, input.pb, proof.com_w,
                                          input.public_w, input.g_offset);
    return ParallelR1cs<Policy>::Verify(proof.r1cs_proof, seed, pr_input);
  }

  static bool Test();
};

// save to bin
template <typename Ar>
void serialize(Ar& ar, Substr<groth09::OrdinaryPolicy>::Proof const& t) {
  ar& YAS_OBJECT_NVP("s.p", ("hp", t.r1cs_proof), ("w", t.com_w));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Substr<groth09::OrdinaryPolicy>::Proof& t) {
  ar& YAS_OBJECT_NVP("s.p", ("hp", t.r1cs_proof), ("w", t.com_w));
}

// save to bin
template <typename Ar>
void serialize(Ar& ar, Substr<groth09::SuccinctPolicy>::Proof const& t) {
  ar& YAS_OBJECT_NVP("s.p", ("hp", t.r1cs_proof), ("w", t.com_w));
}

// load from bin
template <typename Ar>
void serialize(Ar& ar, Substr<groth09::SuccinctPolicy>::Proof& t) {
  ar& YAS_OBJECT_NVP("s.p", ("hp", t.r1cs_proof), ("w", t.com_w));
}

template <typename Policy>
bool Substr<Policy>::Test() {
  auto seed = misc::RandH256();
  int64_t g_offset = 30;

  std::vector<std::string> x_str{{"1234567890123456789012345678901", "",
                                  "abcdefg", "34356356", "234qasdfaq44uUU",
                                  "AQW34352", "ASGDFUR7erweq", "asdfq45343",
                                  "q34asdfasfa9", "4534asq34afa"}};

#ifndef _DEBUG
  while (x_str.size() < 100000) {
    uint8_t buf[33];
    int len = rand() % 32;
    for (int i = 0; i < len; ++i) {
      buf[i] = (uint8_t)rand();
    }
    buf[len] = 0;
    x_str.push_back((char*)buf);
  }
#endif

  int64_t n = (int64_t)x_str.size();
  std::vector<Fr> x(n);
  for (int64_t i = 0; i < n; ++i) {
    x[i] = PackStrToFr(x_str[i].c_str());
  }
  std::string k = "343";
  std::vector<Fr> y(x.size());
  int64_t find_count = 0;
  for (size_t i = 0; i < x_str.size(); ++i) {
    bool find = x_str[i].find(k) != std::string::npos;
    if (find) ++find_count;
    y[i] = find ? 1 : 0;
  }
  std::cout << "find_count: " << find_count << "\n";

  Tick tick(__FUNCTION__);
  Fr com_x_r = FrRand();
  G1 com_x = PcComputeCommitmentG(g_offset, x, com_x_r);
  Fr com_y_r = FrRand();
  G1 com_y = PcComputeCommitmentG(g_offset, y, com_y_r, true);  // y is {0,1}

  ProverInput prover_input(k, x, com_x, com_x_r, y, com_y, com_y_r, g_offset);
  Proof proof;
  Prove(proof, seed, prover_input);

  if (proof.com_x() != com_x || proof.com_y() != com_y) {
    assert(false);
    return false;
  }

#ifndef DISABLE_SERIALIZE_CHECK
  // serialize to buffer
  yas::mem_ostream os;
  yas::binary_oarchive<yas::mem_ostream, YasBinF()> oa(os);
  oa.serialize(proof);
  std::cout << "proof size: " << os.get_shared_buffer().size << "\n";
  // serialize from buffer
  yas::mem_istream is(os.get_intrusive_buffer());
  yas::binary_iarchive<yas::mem_istream, YasBinF()> ia(is);
  Proof proof2;
  ia.serialize(proof2);
  if (proof != proof2) {
    assert(false);
    std::cout << "oops, serialize check failed\n";
    return false;
  }
#endif

  VerifierInput verifier_input(n, k, g_offset);
  bool success = Verify(proof, seed, verifier_input);
  std::cout << __FILE__ << " " << __FUNCTION__ << ": " << success << "\n\n\n\n\n\n";
  return success;
}
}  // namespace pc_utils