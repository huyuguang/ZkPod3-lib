#pragma once

namespace clink::vgg16 {
// using Policy = groth09::OrdinaryPolicy;
using Policy = groth09::SuccinctPolicy;
using R1cs = typename clink::ParallelR1cs<Policy>;
using HyraxA = typename Policy::HyraxA;
using Sec51 = typename Policy::Sec51;
using Sec53 = typename Policy::Sec53;
using Sec43 = typename Policy::Sec43;

}  // namespace clink::vgg16