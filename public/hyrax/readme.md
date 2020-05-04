# Hyrax

https://eprint.iacr.org/2017/1132.pdf

### 概要
本文仅包含论文附录中的A1、A2以及A3算法。本文不关心原论文中和sumcheck和GKR以及电路相关的部分。  
[基础概念和符号约定](../ecc/readme.md)

-------

### A1

##### 语义

- $x,y,z$是$F_r$中的3个秘密标量
- Prover公开$com(x),com(y),com(z)$并证明$z=xy$，其中有：
  - $X=com(x)=h^{r_x}{g_x}^x$
  - $Y=com(y)=h^{r_y}g_y^y$
  - $Z=com(z)=h^{r_z}{g_x}^z$
  - $g_x,g_y$是$G1$中预定义的两个独立标量。
  - 原论文中$com(x),com(y),com(z)$使用了相同的$g$，但实际上$g_x,g_y$可以不相同。注意$com(z)$所使用的$g$必须和$com(x)$所使用的$g$相同。

##### 协议

协议完全遵循原论文。不同之处在于用$FS$变换把协议转为了非交互的。

1. Prover选取5个$F_r$随机数$b_1,b_2,b_3,b_4,b_5$，并且计算：
   1. $\alpha=h^{b_2}g_x^{b_1}$
   2. $\beta=h^{b_4}g_y^{b_3}$
   3. $\delta=h^{b_5}X^{b3}$
   4. 以上三个值相当于锁定了$b_1,b_2,b_3,b_4,b_5 $
2. Prover基于$FS$变换产生挑战随机数$c=FS(X,Y,Z,\alpha,\beta,\delta)$
3. Prover计算：
   1. $z_1=b_1+cx$
   2. $z_2=b_2+cr_x$
   3. $z_3=b_3+cy$
   4. $z_4=b_4+cr_y$
   5. $z_5=b_5+c(r_z-r_xy)$
4. Prover发送$Proof={\alpha,\beta,\delta,z_1,z_2,z_3,z_4,z_5}$
5. Verifier验证如下等式是否成立：
   1. ${\alpha}X^c \overset ? = h^{z_2}g_x{z_1}$
   2. $\beta Y^c \overset ? = h^{z_4}g_y^{z_3}$
   3. $\delta Z^c \overset ? = X^{z_3}h^{z_5}$

##### 性能

- 证明开销：$6EccExp$
- 验证开销：$6EccExp$
- 证据长度：$3G1+5F_r$

##### 说明

Prover选取5个随机数，并且通过$\alpha,\beta,\delta$锁定。Verifier做的前2个验证，是在验证$z1,z2,z3,z4$的定义。Verifiier做的第3个判定验证$z_5$的定义，同时验证$com(zc)=com(xyc)$也即$zc=xyc $。只要第三个验证成立，那么$z5$只能用这个方法产生，且必然有$zc=xyc$，除非Prover能破解$dlp$问题。

##### 备注

据说这个方法来源于Bayer's PhD Thesis，改自于Chaum&Pedersen[CP92]。

---------

### A2
##### 语义

- $\vec a$是公开的$F_r$向量，向量长度为$n$
- $\vec x$是秘密的$F_r$向量，向量长度为$n$
- $y$是秘密的$F_r$标量，且$y=\vec x \cdot \vec a$
- Prover公开$com(x)、com(y)$并证明$y=\vec x \cdot \vec a$
  - $com(\vec x)=h^{r_x}\vec{g_x}^{\vec x}$
  - $com(y)=h^{r_y}g_y^y$
  - $\vec{g_x}, g_y$可以重叠

##### 协议

协议完全遵循原论文。不同之处在于用$FS$变换把协议转为了非交互的。

1. TODO

##### 性能

- 证明开销：$MulExp(n)$
- 验证开销：$MulExp(n)$
- 证据长度：$2G1+(n+2)F_r$

-------

### A3
##### 语义

A3是A2的递归折叠版本，有更大的证明开销和验证开销，同时有对数级别的证明长度。

##### 协议

1. TODO

##### 性能

- 证明开销： $\sim 3MulExp(n)$

- 验证开销： $2nEccExp$
- 证据长度：$(2+2\lg(n))G1+2F_r$