#pragma once

#include "./policy.h"

namespace clink::vgg16 {

struct ConvLayerInfo {
  ConvLayerInfo(size_t order, size_t dimension, size_t channel_count,
                size_t kernel_count)
      : order(order),
        dimension(dimension),
        channel_count(channel_count),
        kernel_count(kernel_count) {}
  size_t const order;
  size_t const dimension;
  size_t const channel_count;
  size_t const kernel_count;
  size_t input_count() const { return dimension * dimension * channel_count; }
  size_t output_count() const { return dimension * dimension * kernel_count; }
};

struct BnLayerInfo {
  BnLayerInfo(size_t order, size_t dimension, size_t channel_count)
      : order(order), dimension(dimension), channel_count(channel_count) {}
  size_t const order;
  size_t const dimension;
  size_t const channel_count;
};

struct DenseLayerInfo {
  DenseLayerInfo(size_t order, size_t input_count, size_t type_count)
      : order(order), input_count(input_count), type_count(type_count) {}
  size_t const order;
  size_t const input_count;
  size_t const type_count;
};

struct ImageInfo {
  ImageInfo(size_t order, size_t dimension, size_t channel_count)
      : order(order), dimension(dimension), channel_count(channel_count) {}
  size_t const order;
  size_t const dimension;
  size_t const channel_count;
};

inline const std::array<ConvLayerInfo, 13> kConvLayerInfos{
    ConvLayerInfo(0, 32, 3, 64),    ConvLayerInfo(1, 32, 64, 64),
    ConvLayerInfo(2, 16, 64, 128),  ConvLayerInfo(3, 16, 128, 128),
    ConvLayerInfo(4, 8, 128, 256),  ConvLayerInfo(5, 8, 256, 256),
    ConvLayerInfo(6, 8, 256, 256),  ConvLayerInfo(7, 4, 256, 512),
    ConvLayerInfo(8, 4, 512, 512),  ConvLayerInfo(9, 4, 512, 512),
    ConvLayerInfo(10, 2, 512, 512), ConvLayerInfo(11, 2, 512, 512),
    ConvLayerInfo(12, 2, 512, 512)};

inline const std::array<BnLayerInfo, 14> kBnLayerInfos{
    BnLayerInfo(0, 32, 64),  BnLayerInfo(1, 32, 64),  BnLayerInfo(2, 16, 128),
    BnLayerInfo(3, 16, 128), BnLayerInfo(4, 8, 256),  BnLayerInfo(5, 8, 256),
    BnLayerInfo(6, 8, 256),  BnLayerInfo(7, 4, 512),  BnLayerInfo(8, 4, 512),
    BnLayerInfo(9, 4, 512),  BnLayerInfo(10, 2, 512), BnLayerInfo(11, 2, 512),
    BnLayerInfo(12, 2, 512), BnLayerInfo(13, 1, 512)};

inline const std::array<DenseLayerInfo, 2> kDenseLayerInfos{
    DenseLayerInfo(0, 512, 512), DenseLayerInfo(1, 512, 10)};

inline const std::array<ImageInfo, 35> kImageInfos{
    ImageInfo(0, 32, 3),   ImageInfo(1, 32, 64),  ImageInfo(2, 32, 64),
    ImageInfo(3, 32, 64),  ImageInfo(4, 32, 64),  ImageInfo(5, 16, 64),
    ImageInfo(6, 16, 128), ImageInfo(7, 16, 128), ImageInfo(8, 16, 128),
    ImageInfo(9, 16, 128), ImageInfo(10, 8, 128), ImageInfo(11, 8, 256),
    ImageInfo(12, 8, 256), ImageInfo(13, 8, 256), ImageInfo(14, 8, 256),
    ImageInfo(15, 8, 256), ImageInfo(16, 8, 256), ImageInfo(17, 4, 256),
    ImageInfo(18, 4, 512), ImageInfo(19, 4, 512), ImageInfo(20, 4, 512),
    ImageInfo(21, 4, 512), ImageInfo(22, 4, 512), ImageInfo(23, 4, 512),
    ImageInfo(24, 2, 512), ImageInfo(25, 2, 512), ImageInfo(26, 2, 512),
    ImageInfo(27, 2, 512), ImageInfo(28, 2, 512), ImageInfo(29, 2, 512),
    ImageInfo(30, 2, 512), ImageInfo(31, 1, 512), ImageInfo(32, 1, 512),
    ImageInfo(33, 1, 512), ImageInfo(34, 1, 10)};

enum LayerType { kConv, kReluBn, kPooling, kDense };

inline const std::array<std::pair<LayerType, size_t>, 34> kLayerTypeOrders{
    {{kConv, 0},    {kReluBn, 0},  {kConv, 1},    {kReluBn, 1},  {kPooling, 0},
     {kConv, 2},    {kReluBn, 2},  {kConv, 3},    {kReluBn, 3},  {kPooling, 1},
     {kConv, 4},    {kReluBn, 4},  {kConv, 5},    {kReluBn, 5},  {kConv, 6},
     {kReluBn, 6},  {kPooling, 2}, {kConv, 7},    {kReluBn, 7},  {kConv, 8},
     {kReluBn, 8},  {kConv, 9},    {kReluBn, 9},  {kPooling, 3}, {kConv, 10},
     {kReluBn, 10}, {kConv, 11},   {kReluBn, 11}, {kConv, 12},   {kReluBn, 12},
     {kPooling, 4}, {kDense, 0},   {kReluBn, 13}, {kDense, 1}}};

inline const std::array<size_t, 13> kConvLayers{0,  2,  5,  7,  10, 12, 14,
                                                17, 19, 21, 24, 26, 28};

inline const std::array<size_t, 2> kDenseLayers{31, 33};

inline const std::array<size_t, 14> kReluBnLayers{1,  3,  6,  8,  11, 13, 15,
                                                  18, 20, 22, 25, 27, 29, 32};

inline const std::array<size_t, 5> kPoolingLayers{4, 9, 16, 23, 30};
}  // namespace clink::vgg16