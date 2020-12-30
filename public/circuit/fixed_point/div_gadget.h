#pragma once

#include "./abs_gadget.h"

namespace circuit::fixed_point {

// ret = a / b, c = a / b ... d
// C * B + D * 2^N == A * 2^N && D < B && C is valid
// a, b must >=0
template <size_t D, size_t N>
class DivBaseGadget : public libsnark::gadget<Fr> {
  static_assert(2 * D + 2 * N < 253, "invalid D,N");

 public:
  DivBaseGadget(libsnark::protoboard<Fr>& pb,
                libsnark::pb_linear_combination<Fr> const& a,
                libsnark::pb_linear_combination<Fr> const& b,
                const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix), a_(a), b_(b) {
    auto constances = RationalConst<D, N>();
    c_.allocate(pb, FMT(this->annotation_prefix, " c"));
    dn_.allocate(pb, FMT(this->annotation_prefix, " dn"));

    libsnark::pb_linear_combination<Fr> bn;
    bn.assign(this->pb, b_ * constances.kFrN);

    less_.allocate(pb, FMT(this->annotation_prefix, " less"));
    less_or_eq_.allocate(pb, FMT(this->annotation_prefix, " less_or_eq"));
    comparison_gadget_.reset(new libsnark::comparison_gadget<Fr>(
        pb, D + N + N, dn_, bn, less_, less_or_eq_,
        FMT(this->annotation_prefix, " comparison")));

    // even we known c >=0, we need check if c is valid
    sign_gadget_.reset(new SignGadget<D, N>(
        pb, c_, FMT(this->annotation_prefix, " type_gadget")));
  }

  void generate_r1cs_constraints() {
    auto constances = RationalConst<D, N>();
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(c_, b_, a_ * constances.kFrN - dn_),
        FMT(this->annotation_prefix, " C * B + D * 2^N == A * 2^N"));
    comparison_gadget_->generate_r1cs_constraints();
    this->pb.add_r1cs_constraint(
        libsnark::r1cs_constraint<Fr>(less_, libsnark::pb_variable<Fr>(0),
                                      libsnark::pb_variable<Fr>(0)),
        FMT(this->annotation_prefix, " less == 1"));
    sign_gadget_->generate_r1cs_constraints();
  }

  void generate_r1cs_witness() {
    auto constances = RationalConst<D, N>();
    a_.evaluate(this->pb);
    b_.evaluate(this->pb);

    Fr fr_a = this->pb.lc_val(a_);
    mpz_class mpz_a = FrToSignedMpz(fr_a);
    Fr fr_b = this->pb.lc_val(b_);
    mpz_class mpz_b = FrToSignedMpz(fr_b);
    assert(!fr_a.isNegative() && !fr_b.isNegative());

    mpz_class mpz_c = mpz_a * constances.kMpzN / mpz_b;
    mpz_class mpz_dn = mpz_a * constances.kMpzN - mpz_b * mpz_c;
    Fr fr_c = SignedMpzToFr(mpz_c);
    this->pb.val(c_) = fr_c;
    Fr fr_dn = SignedMpzToFr(mpz_dn);
    this->pb.val(dn_) = fr_dn;

    // std::cout << "fr_a: " << fr_a << "\n";
    // std::cout << "fr_b: " << fr_b << "\n";
    // std::cout << "fr_c: " << fr_c << "\n";
    auto num = libff::bigint<Fr::num_limbs>(mpz_c.get_mpz_t()).num_bits();
    (void)num;

    comparison_gadget_->generate_r1cs_witness();

    sign_gadget_->generate_r1cs_witness();
  }

  libsnark::pb_variable<Fr> ret() const { return c_; }

  libsnark::pb_variable<Fr> sign() const { return sign_gadget_->ret(); }

  static bool Test(double double_a, double double_b);

 private:
  libsnark::pb_linear_combination<Fr> a_;
  libsnark::pb_linear_combination<Fr> b_;
  libsnark::pb_variable<Fr> c_;
  libsnark::pb_variable<Fr> dn_;
  libsnark::pb_variable<Fr> less_;
  libsnark::pb_variable<Fr> less_or_eq_;
  std::unique_ptr<libsnark::comparison_gadget<Fr>> comparison_gadget_;
  std::unique_ptr<SignGadget<D, N>> sign_gadget_;
};

template <size_t D, size_t N>
bool DivBaseGadget<D, N>::Test(double double_a, double double_b) {
  Tick tick(__FN__);
  auto double_ret = double_a / double_b;
  std::cout << Tick::GetIndentString() << "double_ret: " << double_ret << "\n";

  Fr a = DoubleToRational<D, N>(double_a);
  Fr b = DoubleToRational<D, N>(double_b);

  libsnark::protoboard<Fr> pb;
  libsnark::pb_variable<Fr> pb_a;
  libsnark::pb_variable<Fr> pb_b;
  pb_a.allocate(pb, "TestDiv a");
  pb_b.allocate(pb, "TestDiv b");
  DivBaseGadget<D, N> gadget(pb, pb_a, pb_b, "TestDiv");
  gadget.generate_r1cs_constraints();
  pb.val(pb_a) = a;
  pb.val(pb_b) = b;
  gadget.generate_r1cs_witness();
  CHECK(pb.is_satisfied(), "");
  if (!pb.is_satisfied()) return false;

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

inline bool TestDivBase() {
  Tick tick(__FN__);
  double double_a = 7.3;
  double double_b = 32.1;
  return DivBaseGadget<32, 32>::Test(double_a, double_b);
}

/////
template <size_t D, size_t N>
class DivGadget : public libsnark::gadget<Fr> {
  static_assert(2 * D + 2 * N < 253, "invalid D,N");

 public:
  DivGadget(libsnark::protoboard<Fr>& pb,
            libsnark::pb_linear_combination<Fr> const& a,
            libsnark::pb_linear_combination<Fr> const* sign_a,
            libsnark::pb_linear_combination<Fr> const& b,
            libsnark::pb_linear_combination<Fr> const* sign_b,
            const std::string& annotation_prefix = "")
      : libsnark::gadget<Fr>(pb, annotation_prefix), a_(a), b_(b) {
    if (sign_a) {
      sign_a_ = *sign_a;
    } else {
      sign_a_gadget_.reset(new SignGadget<D, N>(
          this->pb, a_, FMT(this->annotation_prefix, " sign_a_gadget")));
      sign_a_ = sign_a_gadget_->ret();
    }

    if (sign_b) {
      sign_b_ = *sign_b;
    } else {
      sign_b_gadget_.reset(new SignGadget<D, N>(
          this->pb, b_, FMT(this->annotation_prefix, " sign_b_gadget")));
      sign_b_ = sign_b_gadget_->ret();
    }

    if (sign_a_.is_constant()) {
      libsnark::linear_combination<Fr> abs_a;
      sign_a_.evaluate(pb);
      if (pb.lc_val(sign_a_) == 1) {
        abs_a = a_;
      } else {
        abs_a = -a_;
      }
      lc_abs_a_.assign(pb, abs_a);
    } else {
      abs_a_.allocate(pb, FMT(this->annotation_prefix, " abs_a"));
      lc_abs_a_.assign(pb, abs_a_);
    }

    if (sign_b_.is_constant()) {
      libsnark::linear_combination<Fr> abs_b;
      sign_b_.evaluate(pb);
      if (pb.lc_val(sign_b_) == 1) {
        abs_b = b_;
      } else {
        abs_b = -b_;
      }
      lc_abs_b_.assign(pb, abs_b);
    } else {
      abs_b_.allocate(pb, FMT(this->annotation_prefix, " abs_b"));
      lc_abs_b_.assign(pb, abs_b_);
    }

    divbase_gadget_.reset(
        new DivBaseGadget<D, N>(this->pb, lc_abs_a_, lc_abs_b_,
                                FMT(this->annotation_prefix, " div_gadget")));

    if (sign_a_.is_constant() && sign_b_.is_constant()) {
      bool sa = pb.lc_val(sign_a_) == FrOne() ? true : false;
      bool sb = pb.lc_val(sign_b_) == FrOne() ? true : false;
      bool s = !(sa ^ sb);
      lc_sign_ret_.assign(pb, s ? FrOne() : FrZero());
      lc_ret_.assign(pb, divbase_gadget_->ret() * (s ? 1 : -1));
    } else {
      sign_ret_.allocate(pb, FMT(this->annotation_prefix, " sign_ret"));
      lc_sign_ret_.assign(pb, sign_ret_);
      ret_.allocate(pb, FMT(this->annotation_prefix, " ret"));
      lc_ret_.assign(pb, ret_);
    }
  }

  void generate_r1cs_constraints() {
    if (sign_a_gadget_) {
      sign_a_gadget_->generate_r1cs_constraints();
    }

    if (sign_b_gadget_) {
      sign_b_gadget_->generate_r1cs_constraints();
    }

    if (!sign_a_.is_constant()) {
      // abs_a = sign_a? a:-a
      this->pb.add_r1cs_constraint(
          libsnark::r1cs_constraint<Fr>(sign_a_, a_ * 2, a_ + abs_a_),
          FMT(this->annotation_prefix, " abs_a = sign_a? a:-a"));
    }

    if (!sign_b_.is_constant()) {
      // abs_b = sign_b? b:-b
      this->pb.add_r1cs_constraint(
          libsnark::r1cs_constraint<Fr>(sign_b_, b_ * 2, b_ + abs_b_),
          FMT(this->annotation_prefix, " abs_b = sign_b? b:-b"));
    }

    if (!sign_a_.is_constant() || !sign_b_.is_constant()) {
      // sign_ret = 2 * sign_a * sign_b + 1 - sign_a - sign_b
      this->pb.add_r1cs_constraint(
          libsnark::r1cs_constraint<Fr>(sign_a_, sign_b_ * 2,
                                        sign_a_ + sign_b_ - 1 + sign_ret_),
          FMT(this->annotation_prefix,
              " sign_ret = 2 * sign_a * sign_b + 1 - sign_a - sign_b"));

      // ret = sign_ret? positive_ret: -positive_ret
      auto positive_ret = divbase_gadget_->ret();
      this->pb.add_r1cs_constraint(
          libsnark::r1cs_constraint<Fr>(sign_ret_, positive_ret * 2,
                                        positive_ret + ret_),
          FMT(this->annotation_prefix,
              " ret = sign_ret? positive_ret: -positive_ret"));
    }

    divbase_gadget_->generate_r1cs_constraints();
  }

  void generate_r1cs_witness() {
    a_.evaluate(this->pb);
    b_.evaluate(this->pb);

    if (sign_a_gadget_) {
      sign_a_gadget_->generate_r1cs_witness();
    }
    if (sign_b_gadget_) {
      sign_b_gadget_->generate_r1cs_witness();
    }

    sign_a_.evaluate(this->pb);
    sign_b_.evaluate(this->pb);

    Fr a = this->pb.lc_val(a_);
    Fr sign_a = this->pb.lc_val(sign_a_);
    Fr check_sign_a = a.isNegative() ? 0 : 1;
    CHECK(sign_a == check_sign_a, "");
    if (!sign_a_.is_constant()) {
      this->pb.val(abs_a_) = a.isNegative() ? -a : a;
    }
    lc_abs_a_.evaluate(this->pb);
    // std::cout << "lc_abs_a: " << this->pb.lc_val(lc_abs_a_) << "\n";

    Fr b = this->pb.lc_val(b_);
    Fr sign_b = this->pb.lc_val(sign_b_);
    // std::cout << "sign_b: " << sign_b << "\n";
    Fr check_sign_b = b.isNegative() ? 0 : 1;
    CHECK(sign_b == check_sign_b, "");
    if (!sign_b_.is_constant()) {
      this->pb.val(abs_b_) = b.isNegative() ? -b : b;
    }
    lc_abs_b_.evaluate(this->pb);
    // std::cout << "lc_abs_b: " << this->pb.lc_val(lc_abs_b_) << "\n";

    if (!sign_a_.is_constant() || !sign_b_.is_constant()) {
      this->pb.val(sign_ret_) = sign_a * sign_b * 2 + 1 - sign_a - sign_b;
    }
    lc_sign_ret_.evaluate(this->pb);

    divbase_gadget_->generate_r1cs_witness();

    if (!sign_a_.is_constant() || !sign_b_.is_constant()) {
      Fr positive_ret = this->pb.val(divbase_gadget_->ret());
      // std::cout << "positive_ret: " << positive_ret << "\n";
      this->pb.val(ret_) =
          this->pb.val(sign_ret_) == 1 ? positive_ret : -positive_ret;
    }
    lc_ret_.evaluate(this->pb);

    // std::cout << "ret: " << this->pb.lc_val(lc_ret_) << "\n";
  }

  libsnark::pb_linear_combination<Fr> ret() const { return lc_ret_; }

  libsnark::pb_linear_combination<Fr> sign() const { return lc_sign_ret_; }

  enum { kNUllSign, kConstSign, kVarSign };
  static bool Test(double double_a, double double_b, int sign_a_flag,
                   int sign_b_flag);

 private:
  libsnark::pb_linear_combination<Fr> a_;
  libsnark::pb_linear_combination<Fr> sign_a_;
  std::unique_ptr<SignGadget<D, N>> sign_a_gadget_;
  libsnark::pb_linear_combination<Fr> b_;
  libsnark::pb_linear_combination<Fr> sign_b_;
  std::unique_ptr<SignGadget<D, N>> sign_b_gadget_;
  libsnark::pb_variable<Fr> abs_a_;
  libsnark::pb_linear_combination<Fr> lc_abs_a_;
  libsnark::pb_variable<Fr> abs_b_;
  libsnark::pb_linear_combination<Fr> lc_abs_b_;
  libsnark::pb_variable<Fr> sign_ret_;
  libsnark::pb_linear_combination<Fr> lc_sign_ret_;
  libsnark::pb_variable<Fr> ret_;
  libsnark::pb_linear_combination<Fr> lc_ret_;
  std::unique_ptr<DivBaseGadget<D, N>> divbase_gadget_;
};

template <size_t D, size_t N>
bool DivGadget<D, N>::Test(double double_a, double double_b, int sign_a_flag,
                           int sign_b_flag) {
  Tick tick(__FN__);
  auto double_ret = double_a / double_b;
  std::cout << Tick::GetIndentString() << "double_ret: " << double_ret << "\n";

  Fr a = DoubleToRational<D, N>(double_a);
  Fr b = DoubleToRational<D, N>(double_b);

  libsnark::protoboard<Fr> pb;
  libsnark::pb_variable<Fr> pb_a;
  libsnark::pb_variable<Fr> pb_b;
  libsnark::pb_variable<Fr> pb_sign_a;
  libsnark::pb_variable<Fr> pb_sign_b;
  std::shared_ptr<libsnark::pb_linear_combination<Fr>> pb_lc_sign_a;
  std::shared_ptr<libsnark::pb_linear_combination<Fr>> pb_lc_sign_b;

  pb_a.allocate(pb, "a");
  pb_b.allocate(pb, "b");

  if (sign_a_flag == kVarSign) {
    pb_sign_a.allocate(pb, "sign_a");
    pb_lc_sign_a.reset(new libsnark::pb_linear_combination<Fr>(pb_sign_a));
  } else if (sign_a_flag == kConstSign) {
    pb_lc_sign_a.reset(new libsnark::pb_linear_combination<Fr>());
    pb_lc_sign_a->assign(pb, double_a >= 0 ? Fr(1) : Fr(0));
  }

  if (sign_b_flag == kVarSign) {
    pb_sign_b.allocate(pb, "sign_b");
    pb_lc_sign_b.reset(new libsnark::pb_linear_combination<Fr>(pb_sign_b));
  } else if (sign_b_flag == kConstSign) {
    pb_lc_sign_b.reset(new libsnark::pb_linear_combination<Fr>());
    pb_lc_sign_b->assign(pb, double_b >= 0 ? Fr(1) : Fr(0));
  }

  DivGadget<D, N> gadget(pb, pb_a, pb_lc_sign_a.get(), pb_b, pb_lc_sign_b.get(),
                         "TestDiv2");
  gadget.generate_r1cs_constraints();

  pb.val(pb_a) = a;
  pb.val(pb_b) = b;
  if (sign_a_flag == kVarSign) {
    pb.val(pb_sign_a) = Fr(double_a >= 0 ? 1 : 0);
  }

  if (sign_b_flag == kVarSign) {
    pb.val(pb_sign_b) = Fr(double_b >= 0 ? 1 : 0);
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

inline bool TestDiv() {
  Tick tick(__FN__);
  enum { kNUllSign, kConstSign, kVarSign };
  bool ret;
  std::vector<bool> rets;
  double double_a;
  double double_b;

  double_a = -7.3;
  double_b = -32.1;
  ret = DivGadget<32, 32>::Test(double_a, double_b, kNUllSign, kNUllSign);
  rets.push_back(ret);
  ret = DivGadget<32, 32>::Test(double_a, double_b, kVarSign, kNUllSign);
  rets.push_back(ret);
  ret = DivGadget<32, 32>::Test(double_a, double_b, kNUllSign, kConstSign);
  rets.push_back(ret);
  ret = DivGadget<32, 32>::Test(double_a, double_b, kConstSign, kVarSign);
  rets.push_back(ret);
  ret = DivGadget<32, 32>::Test(double_a, double_b, kVarSign, kConstSign);
  rets.push_back(ret);

  double_a = 7.3;
  double_b = -32.1;
  ret = DivGadget<32, 32>::Test(double_a, double_b, kNUllSign, kNUllSign);
  rets.push_back(ret);
  ret = DivGadget<32, 32>::Test(double_a, double_b, kVarSign, kNUllSign);
  rets.push_back(ret);
  ret = DivGadget<32, 32>::Test(double_a, double_b, kNUllSign, kConstSign);
  rets.push_back(ret);
  ret = DivGadget<32, 32>::Test(double_a, double_b, kConstSign, kVarSign);
  rets.push_back(ret);
  ret = DivGadget<32, 32>::Test(double_a, double_b, kVarSign, kConstSign);
  rets.push_back(ret);

  double_a = -7.3;
  double_b = 32.1;
  ret = DivGadget<32, 32>::Test(double_a, double_b, kNUllSign, kNUllSign);
  rets.push_back(ret);
  ret = DivGadget<32, 32>::Test(double_a, double_b, kVarSign, kNUllSign);
  rets.push_back(ret);
  ret = DivGadget<32, 32>::Test(double_a, double_b, kNUllSign, kConstSign);
  rets.push_back(ret);
  ret = DivGadget<32, 32>::Test(double_a, double_b, kConstSign, kVarSign);
  rets.push_back(ret);
  ret = DivGadget<32, 32>::Test(double_a, double_b, kVarSign, kConstSign);
  rets.push_back(ret);

  double_a = 7.3;
  double_b = 32.1;
  ret = DivGadget<32, 32>::Test(double_a, double_b, kNUllSign, kNUllSign);
  rets.push_back(ret);
  ret = DivGadget<32, 32>::Test(double_a, double_b, kVarSign, kNUllSign);
  rets.push_back(ret);
  ret = DivGadget<32, 32>::Test(double_a, double_b, kNUllSign, kConstSign);
  rets.push_back(ret);
  ret = DivGadget<32, 32>::Test(double_a, double_b, kConstSign, kVarSign);
  rets.push_back(ret);
  ret = DivGadget<32, 32>::Test(double_a, double_b, kVarSign, kConstSign);
  rets.push_back(ret);

  return std::all_of(rets.begin(), rets.end(), [](auto i) { return i; });
}

///
template <size_t D, size_t N>
class InvGadget : public DivGadget<D, N> {
 public:
  InvGadget(libsnark::protoboard<Fr>& pb,
            libsnark::pb_linear_combination<Fr> const& b,
            libsnark::pb_linear_combination<Fr> const* sign_b,
            const std::string& annotation_prefix = "")
      : DivGadget<D, N>(pb, GetA(pb), GetSignA(pb), b, sign_b,
                        annotation_prefix) {}

  enum { kNUllSign, kConstSign, kVarSign };
  static bool Test(double double_b, int sign_b_flag);

 private:
  libsnark::pb_linear_combination<Fr> const& GetA(
      libsnark::protoboard<Fr>& pb) {
    a_.assign(pb, RationalConst<D, N>().kFrN);
    return a_;
  }

  libsnark::pb_linear_combination<Fr>* GetSignA(libsnark::protoboard<Fr>& pb) {
    sign_a_.assign(pb, FrOne());
    return &sign_a_;
  }
  libsnark::pb_linear_combination<Fr> a_;
  libsnark::pb_linear_combination<Fr> sign_a_;
};

template <size_t D, size_t N>
bool InvGadget<D, N>::Test(double double_b, int sign_b_flag) {
  Tick tick(__FN__);
  auto double_ret = 1.0 / double_b;
  std::cout << Tick::GetIndentString() << "double_ret: " << double_ret << "\n";

  Fr b = DoubleToRational<D, N>(double_b);

  libsnark::protoboard<Fr> pb;
  libsnark::pb_variable<Fr> pb_b;
  libsnark::pb_variable<Fr> pb_sign_b;
  std::shared_ptr<libsnark::pb_linear_combination<Fr>> pb_lc_sign_b;

  pb_b.allocate(pb, "b");

  if (sign_b_flag == kVarSign) {
    pb_sign_b.allocate(pb, "sign_b");
    pb_lc_sign_b.reset(new libsnark::pb_linear_combination<Fr>(pb_sign_b));
  } else if (sign_b_flag == kConstSign) {
    pb_lc_sign_b.reset(new libsnark::pb_linear_combination<Fr>());
    pb_lc_sign_b->assign(pb, double_b >= 0 ? Fr(1) : Fr(0));
  }

  InvGadget<D, N> gadget(pb, pb_b, pb_lc_sign_b.get(), "TestDiv");
  gadget.generate_r1cs_constraints();

  pb.val(pb_b) = b;
  if (sign_b_flag == kVarSign) {
    pb.val(pb_sign_b) = Fr(double_b >= 0 ? 1 : 0);
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

inline bool TestInv() {
  Tick tick(__FN__);
  bool ret;
  std::vector<bool> rets;

  enum { kNUllSign, kConstSign, kVarSign };
  double double_b = -1.34;
  ret = InvGadget<32, 32>::Test(double_b, kNUllSign);
  rets.push_back(ret);
  ret = InvGadget<32, 32>::Test(double_b, kVarSign);
  rets.push_back(ret);
  ret = InvGadget<32, 32>::Test(double_b, kConstSign);
  rets.push_back(ret);

  double_b = 1.34;
  ret = InvGadget<32, 32>::Test(double_b, kNUllSign);
  rets.push_back(ret);
  ret = InvGadget<32, 32>::Test(double_b, kVarSign);
  rets.push_back(ret);
  ret = InvGadget<32, 32>::Test(double_b, kConstSign);
  rets.push_back(ret);

  return std::all_of(rets.begin(), rets.end(), [](auto i) { return i; });
}

}  // namespace circuit::fixed_point