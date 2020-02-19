# Hyrax

https://eprint.iacr.org/2017/1132.pdf

### 概要
本文仅包含论文附录中的A1、A2以及A3算法。本文不关心原论文中和sumcheck和GKR以及电路相关的部分。  
[基础概念和符号约定](../ecc/readme.md)

### A1
- x、y、z是F<sub>r</sub>中3个秘密的标量
- Prover公开com(x)、com(y)、com(z)并证明z=<x,y>
  - com(x)=h<sup>r<sub>x</sub></sup>g<sub>x</sub><sup>x</sup>
  - com(y)=h<sup>r<sub>y</sub></sup>g<sub>y</sub><sup>y</sup>
  - com(z)=h<sup>r<sub>z</sub></sup>g<sub>z</sub><sup>z</sup>
  - g<sub>x</sub>和g<sub>y</sub>可以相同

### A2
- a是公开的F<sub>r</sub>向量，向量长度为n
- x是秘密的F<sub>r</sub>向量，向量长度为n
- y是秘密的F<sub>r</sub>标量，且y=<x,a>
- Prover公开com(x)、com(y)并证明y=<x,a>
  - com(x)=h<sup>r<sub>x</sub></sup>g<sub>x</sub><sup>x</sup>
  - com(y)=h<sup>r<sub>y</sub></sup>g<sub>y</sub><sup>y</sup>
  - g<sub>x</sub>（向量）和g<sub>y</sub>（标量）可以重叠

### A3
A3是A2的递归折叠版本，有更大的证明开销和验证开销，同时有对数级别的证明长度。

### 性能
- A1
  - Prove cost：6EccExp
  - Verify cost：6EccExp
  - Proof size：3G1+5F<sub>r</sub>
- A2
  - Prove cost：MulExp(n)
  - Verify cost：MulExp(n)
  - Proof size：2G1+(n+2)F<sub>r</sub>
- A3
  - Prove cost： 2MulExp(n)~4MulExp(n)
  - Verify cost： 2nEccExp
  - Proof size：(2+2log(n))G1+2F<sub>r</sub>
