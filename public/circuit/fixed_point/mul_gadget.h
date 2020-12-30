#pragma once

#include "./abs_gadget.h"

namespace circuit::fixed_point {

// ret = a[0]*a[1]...*a[i]*c[0]*c[1]...

template <size_t D, size_t N>
class MulGadget : public libsnark::gadget<Fr> {
 public:
  MulGadget(libsnark::protoboard<Fr>& pb,
            libsnark::pb_linear_combination_array<Fr> const& a,
            std::vector<Fr> const& c, const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix),
        a_(a),
        c_(c),
        final_c_(FrOne()) {
    CheckMaxNumOfBits();

    auto constances = RationalConst<D, N>();

    for (auto& i : c_) {
      final_c_ *= i;
    }

    Fr kFrDxN = constances.GetFrDxN(a_.size() + c_.size());
    if (a_.size() > 1) {
      products_.allocate(pb, a_.size() - 1,
                         FMT(this->annotation_prefix, " products"));
      libsnark::linear_combination<Fr> last_product =
          libsnark::linear_combination<Fr>(products_[products_.size() - 1] *
                                           final_c_);
      product_off_.assign(this->pb, last_product + kFrDxN);
    } else {
      libsnark::linear_combination<Fr> last_product(a_[0] * final_c_);
      product_off_.assign(this->pb, last_product + kFrDxN);
    }

    bits_.allocate(this->pb, D + (a_.size() + c_.size()) * N + 1,
                   FMT(this->annotation_prefix, " bits"));
    p1_gadget_.reset(new libsnark::packing_gadget<Fr>(
        this->pb, bits_, product_off_, FMT(this->annotation_prefix, " p1")));

    packed_.allocate(pb, FMT(this->annotation_prefix, " packed"));
    p2_gadget_.reset(new libsnark::packing_gadget<Fr>(
        this->pb,
        libsnark::pb_variable_array<Fr>(
            bits_.begin() + (a_.size() + c_.size() - 1) * N, bits_.end()),
        packed_, FMT(this->annotation_prefix, " p2")));

    ret_.assign(this->pb,
                libsnark::linear_combination<Fr>(packed_) - constances.kFrDN);
    sign_.assign(this->pb, p2_gadget_->bits[D + N]);
  }

  void generate_r1cs_constraints() {
    if (a_.size() > 1) {
      this->pb.add_r1cs_constraint(
          libsnark::r1cs_constraint<Fr>(a_[0], a_[1], products_[0]),
          FMT(this->annotation_prefix, " a[0]*a[1] = products[0]"));
      for (size_t i = 2; i < a_.size(); ++i) {
        this->pb.add_r1cs_constraint(
            libsnark::r1cs_constraint<Fr>(a_[i], products_[i - 2],
                                          products_[i - 1]),
            FMT(this->annotation_prefix, " a[%uz]*a[%zu]*c = products[%zu]", i,
                i - 2, i - 1));
      }
    }

    p1_gadget_->generate_r1cs_constraints(true);
    p2_gadget_->generate_r1cs_constraints(false);
  }

  void generate_r1cs_witness() {
    for (auto& i : a_) {
      i.evaluate(this->pb);
    }

    if (!products_.empty()) {
      this->pb.val(products_[0]) =
          this->pb.lc_val(a_[0]) * this->pb.lc_val(a_[1]);

      for (size_t i = 2; i < a_.size(); ++i) {
        this->pb.val(products_[i - 1]) =
            this->pb.lc_val(a_[i]) * this->pb.val(products_[i - 2]);
      }
    }

    product_off_.evaluate(this->pb);

    p1_gadget_->generate_r1cs_witness_from_packed();
    p2_gadget_->generate_r1cs_witness_from_bits();
    ret_.evaluate(this->pb);
    sign_.evaluate(this->pb);
  }

  libsnark::pb_linear_combination<Fr> ret() const { return ret_; }

  // 1: >=0; 0: <=0
  libsnark::pb_linear_combination<Fr> sign() const { return sign_; };

  static bool Test(std::vector<double> const& double_a,
                   std::vector<double> const& double_c);

 private:
  void CheckMaxNumOfBits() {
    CHECK(!a_.empty(), "");

    size_t total_bits_of_c = 0;
    for (auto& i : c_) {
      Fr abs_i = i.isNegative() ? -i : i;
      size_t bits_of_i =
          libff::bigint<Fr::num_limbs>(abs_i.getMpz().get_mpz_t()).num_bits();
      total_bits_of_c += bits_of_i;
    }
    size_t total_bits = a_.size() * (D + N) + total_bits_of_c;

    CHECK(total_bits < 253, "");
  }

 private:
  libsnark::pb_linear_combination_array<Fr> a_;
  std::vector<Fr> c_;
  Fr final_c_;
  libsnark::pb_linear_combination<Fr> ret_;
  libsnark::pb_linear_combination<Fr> sign_;
  libsnark::pb_variable_array<Fr> products_;
  libsnark::pb_linear_combination<Fr> product_off_;  // = product_ + 2^(D+2N)
  std::unique_ptr<libsnark::packing_gadget<Fr>> p1_gadget_;
  std::unique_ptr<libsnark::packing_gadget<Fr>> p2_gadget_;
  libsnark::pb_variable_array<Fr> bits_;
  libsnark::pb_variable<Fr> packed_;
};

template <size_t D, size_t N>
bool MulGadget<D, N>::Test(std::vector<double> const& double_a,
                           std::vector<double> const& double_c) {
  Tick tick(__FN__);
  double double_ret = 1.0;
  for (auto& i : double_a) {
    double_ret = double_ret * i;
  }
  for (auto& i : double_c) {
    double_ret = double_ret * i;
  }

  std::cout << Tick::GetIndentString() << "double_ret: " << double_ret << "\n";

  std::vector<Fr> fr_a(double_a.size());
  for (size_t i = 0; i < double_a.size(); ++i) {
    fr_a[i] = DoubleToRational<D, N>(double_a[i]);
  }
  std::vector<Fr> fr_c(double_c.size());
  for (size_t i = 0; i < double_c.size(); ++i) {
    fr_c[i] = DoubleToRational<D, N>(double_c[i]);
  }

  libsnark::protoboard<Fr> pb;
  libsnark::pb_variable_array<Fr> pb_a;
  pb_a.allocate(pb, double_a.size(), "MulGadget::Test a");

  MulGadget<D, N> gadget(pb, pb_a, fr_c, "MulGadget::Test gadget");

  gadget.generate_r1cs_constraints();

  for (size_t i = 0; i < double_a.size(); ++i) {
    pb.val(pb_a[i]) = fr_a[i];
  }

  gadget.generate_r1cs_witness();
  CHECK(pb.is_satisfied(), "");

  Fr fr_ret = pb.lc_val(gadget.ret());
  std::cout << Tick::GetIndentString() << "fr_ret: " << fr_ret << "\t"
            << RationalToDouble<D, N>(fr_ret) << "\n";
  Fr fr_sign = pb.lc_val(gadget.sign());
  std::cout << Tick::GetIndentString() << "sign: " << fr_sign << "\n";
  CHECK(fr_sign == (double_ret >= 0 ? 1 : 0), "");

  std::cout << Tick::GetIndentString()
            << "num_constraints: " << pb.num_constraints() << "\n";
  std::cout << Tick::GetIndentString()
            << "num_variables: " << pb.num_variables() << "\n";
  return true;
}

inline bool TestMul() {
  Tick tick(__FN__);
  bool ret;
  std::vector<bool> rets;
  std::vector<double> double_a;
  std::vector<double> double_c;

  double_a = std::vector<double>{-7.3, 21.1};
  double_c = std::vector<double>{-2.4, -1.4};
  ret = MulGadget<32, 32>::Test(double_a, double_c);
  rets.push_back(ret);

  double_a = std::vector<double>{-7.3, 21.1};
  double_c = std::vector<double>{};
  ret = MulGadget<32, 32>::Test(double_a, double_c);
  rets.push_back(ret);

  double_a = std::vector<double>{-7.3};
  double_c = std::vector<double>{-2.4};
  ret = MulGadget<32, 32>::Test(double_a, double_c);
  rets.push_back(ret);

  return std::all_of(rets.begin(), rets.end(), [](auto i) { return i; });
}

template <size_t D, size_t N>
class MulVCGadget : public MulGadget<D, N> {
 public:
  MulVCGadget(libsnark::protoboard<Fr>& pb,
              libsnark::pb_linear_combination<Fr> const& a, double const& c,
              const std::string& annotation_prefix = "")
      : MulGadget<D, N>(pb, vec_a(a), vec_c(c), annotation_prefix) {}

 private:
  std::vector<Fr> vec_c(double const& c) {
    return std::vector<Fr>{DoubleToRational<D, N>(c)};
  }
  libsnark::pb_linear_combination_array<Fr> vec_a(
      libsnark::pb_linear_combination<Fr> const& a) {
    libsnark::pb_linear_combination_array<Fr> ret(1);
    ret[0] = a;
    return ret;
  }
};

template <size_t D, size_t N>
class MulVVGadget : public MulGadget<D, N> {
 public:
  MulVVGadget(libsnark::protoboard<Fr>& pb,
              libsnark::pb_linear_combination<Fr> const& a1,
              libsnark::pb_linear_combination<Fr> const& a2,
              const std::string& annotation_prefix = "")
      : MulGadget<D, N>(pb, vec_a(a1, a2), vec_c(), annotation_prefix) {}

 private:
  std::vector<Fr> vec_c() { return std::vector<Fr>{}; }
  libsnark::pb_linear_combination_array<Fr> vec_a(
      libsnark::pb_linear_combination<Fr> const& a1,
      libsnark::pb_linear_combination<Fr> const& a2) {
    libsnark::pb_linear_combination_array<Fr> ret(2);
    ret[0] = a1;
    ret[1] = a2;
    return ret;
  }
};

template <size_t D, size_t N>
class MulVVCGadget : public MulGadget<D, N> {
 public:
  MulVVCGadget(libsnark::protoboard<Fr>& pb,
               libsnark::pb_linear_combination<Fr> const& a1,
               libsnark::pb_linear_combination<Fr> const& a2, double c,
               const std::string& annotation_prefix = "")
      : MulGadget<D, N>(pb, vec_a(a1, a2), vec_c(c), annotation_prefix) {}

 private:
  std::vector<Fr> vec_c(double c) {
    return std::vector<Fr>{DoubleToRational<D, N>(c)};
  }

  libsnark::pb_linear_combination_array<Fr> vec_a(
      libsnark::pb_linear_combination<Fr> const& a1,
      libsnark::pb_linear_combination<Fr> const& a2) {
    libsnark::pb_linear_combination_array<Fr> ret(2);
    ret[0] = a1;
    ret[1] = a2;
    return ret;
  }
};

}  // namespace circuit::fixed_point