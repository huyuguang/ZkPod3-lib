#pragma once

#include <stdlib.h>

#include <iostream>
#include <algorithm>
#include <random>
#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/gadgetlib1/gadgets/routing/as_waksman_routing_gadget.hpp>
#include <libsnark/gadgetlib1/pb_variable.hpp>

namespace circuit {

class PermutationGadget : public libsnark::gadget<Fr> {
 private:
  std::vector<std::vector<libsnark::pb_variable<Fr>>> routed_packets;

  libsnark::as_waksman_topology neighbors;

 private:
  static void MarkVarInvalid(libsnark::pb_variable<Fr> &var) {
    var.index = (size_t)(-1);
  }
  bool IsVarValid(libsnark::pb_variable<Fr> const &var) {
    return var.index != (size_t)(-1);
  }
 public:
  const size_t num_packets;
  const size_t num_columns;

  PermutationGadget(
      libsnark::protoboard<Fr> &pb,
      const std::vector<libsnark::pb_variable<Fr>> &input_packets,
      const std::vector<libsnark::pb_variable<Fr>> &output_packets,
      const std::string &annotation_prefix = "")
      : gadget<Fr>(pb, annotation_prefix),
        num_packets(input_packets.size()),
        num_columns(libsnark::as_waksman_num_columns(num_packets)) {
    neighbors = libsnark::generate_as_waksman_topology(num_packets);
    routed_packets.resize(num_columns + 1);

    for (size_t column_idx = 0; column_idx < num_columns + 1; ++column_idx) {
      routed_packets[column_idx].resize(num_packets);
      for (size_t packet_idx = 0; packet_idx < num_packets; ++packet_idx) {
        MarkVarInvalid(routed_packets[column_idx][packet_idx]);
      }
    }

    /* Two pass allocation. First allocate LHS packets, then for every
       switch either copy over the variables from previously allocated
       to allocate target packets */
    for (size_t packet_idx = 0; packet_idx < num_packets; ++packet_idx) {
      routed_packets[0][packet_idx] = input_packets[packet_idx];
    }

    for (size_t packet_idx = 0; packet_idx < num_packets; ++packet_idx) {
      routed_packets[num_columns][packet_idx] = output_packets[packet_idx];
    }

    size_t saved_count = 0;
    for (size_t column_idx = num_columns - 1; column_idx > 0; --column_idx) {
      for (size_t row_idx = 0; row_idx < num_packets; ++row_idx) {
        if (neighbors[column_idx][row_idx].first ==
            neighbors[column_idx][row_idx].second) {
          auto const &next =
              routed_packets[column_idx + 1]
                            [neighbors[column_idx][row_idx].first];
          if (IsVarValid(next)) {
            /* This is a straight edge, so just copy over the previously
             * allocated subpackets */
            routed_packets[column_idx][row_idx] = next;
            ++saved_count;
          }
        }
      }
    }

    for (size_t column_idx = 0; column_idx < num_columns - 1; ++column_idx) {
      routed_packets[column_idx + 1].resize(num_packets);

      for (size_t row_idx = 0; row_idx < num_packets; ++row_idx) {
        if (neighbors[column_idx][row_idx].first ==
            neighbors[column_idx][row_idx].second) {
          /* This is a straight edge, so just copy over the previously allocated
           * subpackets */
          auto &next = routed_packets[column_idx + 1]
                                     [neighbors[column_idx][row_idx].first];
          assert(!IsVarValid(next) ||
                 next.index == routed_packets[column_idx][row_idx].index);
          if (!IsVarValid(next)) {
            next = routed_packets[column_idx][row_idx];
          } else {
            ++saved_count;
            assert(next.index == routed_packets[column_idx][row_idx].index);
          }
        } else {
          const size_t straight_edge = neighbors[column_idx][row_idx].first;
          const size_t cross_edge = neighbors[column_idx][row_idx].second;

          auto &next_straight = routed_packets[column_idx + 1][straight_edge];
          if (!IsVarValid(next_straight)) {
            next_straight.allocate(
                pb, FMT(annotation_prefix, " routed_packets_%zu_%zu",
                        column_idx + 1, straight_edge));
          } else {
            ++saved_count;
            //printf("check1\n");
          }

          auto &next_cross = routed_packets[column_idx + 1][cross_edge];
          if (!IsVarValid(next_cross)) {
            next_cross.allocate(
                pb, FMT(annotation_prefix, " routed_packets_%zu_%zu",
                        column_idx + 1, cross_edge));
          } else {
            ++saved_count;
            //printf("check2\n");
          }

          ++row_idx; /* skip the next idx, as it to refers to the same packets
                      */
        }
      }
    }

    std::cout << "saved count: " << saved_count << "\n";
  }

  void generate_r1cs_constraints() {

    /* actual routing constraints */
    for (size_t column_idx = 0; column_idx < num_columns; ++column_idx) {
      for (size_t row_idx = 0; row_idx < num_packets; ++row_idx) {
        if (neighbors[column_idx][row_idx].first ==
            neighbors[column_idx][row_idx].second) {
          /* if there is no switch at this position, then just continue with
           * next row_idx */
          continue;
        }

        /* easy case: require that
           (cur-straight_edge)*(cur-cross_edge) = 0 for both
           switch inputs */
        for (size_t switch_input : {row_idx, row_idx + 1}) {
          const size_t straight_edge =
              neighbors[column_idx][switch_input].first;
          const size_t cross_edge = neighbors[column_idx][switch_input].second;

          this->pb.add_r1cs_constraint(
              libsnark::r1cs_constraint<Fr>(
                  routed_packets[column_idx][switch_input] -
                      routed_packets[column_idx + 1][straight_edge],
                  routed_packets[column_idx][switch_input] -
                      routed_packets[column_idx + 1][cross_edge],
                  0),
              FMT(this->annotation_prefix, " easy_route_%zu_%zu", column_idx,
                  switch_input));
        }

        /* we processed both switch inputs at once, so skip the next iteration
         */
        ++row_idx;
      }
    }
  }

  libsnark::integer_permutation compute_permutation(
      std::vector<Fr> const &input, std::vector<Fr> const &output) {
    assert(input.size() == output.size() && !input.empty());
    libsnark::integer_permutation permutation(input.size());
    std::set<size_t> exist_pos;
    for (size_t packet_idx = 0; packet_idx < input.size(); ++packet_idx) {
      auto const &v = input[packet_idx];
      size_t offset = 0;
      for (;;) {
        auto it = std::find(output.begin() + offset, output.end(), v);
        assert(it != output.end());
        offset = std::distance(output.begin(), it);
        if (exist_pos.find(offset) == exist_pos.end()) {
          permutation.set(packet_idx, offset);
          exist_pos.insert(offset);
          break;
        }
        ++offset;
      }
    }

    for (size_t packet_idx = 0; packet_idx < num_packets; ++packet_idx) {
      assert(output[permutation.get(packet_idx)] == input[packet_idx]);
    }
   
    std::cout << "input:\n";
    for (size_t i = 0; i < num_packets; ++i) {
      std::cout << input[i] << "\t";
    }
    std::cout << "\n\n";

    std::cout << "output:\n";
    for (size_t i = 0; i < num_packets; ++i) {
      std::cout << output[i] << "\t";
    }
    std::cout << "\n\n";

    std::cout << "permutation:\n";
    for (size_t i = 0; i < num_packets; ++i) {
      std::cout << permutation.get(i) << "\t";
    }
    std::cout << "\n\n";

    return permutation;
  }

  void generate_r1cs_witness(std::vector<Fr> const &input,
                             std::vector<Fr> const &output) {
    auto permutation = compute_permutation(input, output);

    for (size_t packet_idx = 0; packet_idx < num_packets; ++packet_idx) {
      pb.val(routed_packets[0][packet_idx]) = input[packet_idx];
      pb.val(routed_packets[num_columns][packet_idx]) = output[packet_idx];
    }

    /* do the routing */
    libsnark::as_waksman_routing routing = libsnark::get_as_waksman_routing(permutation);

    for (size_t column_idx = 0; column_idx < num_columns-1; ++column_idx) {
      for (size_t row_idx = 0; row_idx < num_packets; ++row_idx) {
        if (neighbors[column_idx][row_idx].first ==
            neighbors[column_idx][row_idx].second) {
          /* this is a straight edge, so just pass the values forward */
          const size_t next = neighbors[column_idx][row_idx].first;

          this->pb.val(routed_packets[column_idx + 1][next]) =
              this->pb.val(routed_packets[column_idx][row_idx]);
        } else {
          /* route according to the switch bit */
          const bool switch_val = routing[column_idx][row_idx];

          for (size_t switch_input : {row_idx, row_idx + 1}) {
            const size_t straight_edge =
                neighbors[column_idx][switch_input].first;
            const size_t cross_edge =
                neighbors[column_idx][switch_input].second;

            const size_t switched_edge =
                (switch_val ? cross_edge : straight_edge);

            this->pb.val(
                routed_packets[column_idx + 1][switched_edge]) =
                this->pb.val(
                    routed_packets[column_idx][switch_input]);
          }

          /* we processed both switch inputs at once, so skip the next iteration
           */
          ++row_idx;
        }
      }
    }

    printf("routed packets\n");
    for (auto const &i : routed_packets) {
      for (auto const &j : i) {
        std::cout << pb.val(j) << "\t";
      }
      printf("\n");
    }
    printf("\n\n");
  }
};


inline void test_as_waksman_routing_gadget(const size_t num_packets) {
  libsnark::protoboard<Fr> pb;

  std::vector<libsnark::pb_variable<Fr>> input_packets(num_packets);
  std::vector<libsnark::pb_variable<Fr>> output_packets(num_packets);
  for (size_t packet_idx = 0; packet_idx < num_packets; ++packet_idx) {
    input_packets[packet_idx].allocate(pb, FMT("", "input_%zu", packet_idx));
  }
  for (size_t packet_idx = 0; packet_idx < num_packets; ++packet_idx) {
    output_packets[packet_idx].allocate(pb, FMT("", "output_%zu", packet_idx));
  }

  PermutationGadget r(pb, input_packets, output_packets,
                      "main_routing_gadget");
  r.generate_r1cs_constraints();

  std::vector<Fr> input(num_packets);
  std::vector<Fr> output(num_packets);
  for (size_t packet_idx = 0; packet_idx < num_packets; ++packet_idx) {
    input[packet_idx] = packet_idx + 1;
  }
  srand((int)time(nullptr));
  input[rand()%num_packets] = input[rand()%num_packets];
  output = input;
  std::random_device rng;
  std::mt19937 urng(rng());
  std::shuffle(input.begin(), input.end(), urng);

  r.generate_r1cs_witness(input, output);

  printf("positive test\n");
  assert(pb.is_satisfied());

  //printf("negative test\n");
  //pb.val(libsnark::pb_variable<Fr>(10)) = Fr(12345);
  //assert(!pb.is_satisfied());

  printf("num_constraints = %zu, num_variables = %zu\n", pb.num_constraints(),
         pb.get_constraint_system().num_variables());
}
}  // namespace circuit