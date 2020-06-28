#pragma once

namespace clink::vgg16 {
using Policy = groth09::OrdinaryPolicy; // SuccinctPolicy
using R1cs = typename clink::ParallelR1cs<Policy>;
using HyraxA = typename Policy::HyraxA;
}  // namespace clink::vgg16