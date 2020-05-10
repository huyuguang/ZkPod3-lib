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

1. Prover选取长度为$n$的随机向量$\vec d$，随机标量$r_\beta,r_\delta$，并且计算：
   1. $\delta=com(d)=h^{r_\delta}\vec{g_x}^{d}$
   2. $\beta=com(\vec a\cdot \vec d)=h^{r_\beta}g_y^{\vec a \cdot \vec d}$
2. Prover基于$FS$变换得到挑战数$c$。
3. Prover计算并发送：
   1. $\vec z=\vec x c+\vec d$
   2. $z_\delta=c r_{\xi}+r_\delta$
   3. $z_\beta=cr_{\tau}+r_\beta$
4. Verifier检查：
   1. $\xi^c \delta \overset{?}=com(\vec z)=h^{z_\delta}\vec {g_x}^{\vec z}$
   2. $\tau^c \beta \overset{?}=com(\vec z \cdot \vec a)=h^{z_\beta}g_y^{\vec z \cdot \vec a}$

##### 性能

- 证明开销：$MulExp(n)$
- 验证开销：$MulExp(n)$
- 证据长度：$2G1+(n+2)F_r$

-------

### A3
##### 语义

A3是A2的递归折叠版本，有更大的证明开销和验证开销，同时有对数级别的证明长度。有趣的是，同一个协议，在另外一篇稍晚的论文`Halo (Recursive Proof Composition without a Trusted
Setup)` 的`Modified Inner Product`一节中又被重新发明了一次。两者的区别在于`Halo`中有一个实现上的优化，从而使得验证速度更快，和A2相当，因此值得特意提一下。

##### 协议

协议完全遵循原论文。不同之处在于用$FS$变换把协议转为了非交互的，原论文中是一个递归过程，实现时改成了循环。具体过程在此略，可以直接参考原始论文。

##### Halo的优化

在A3中，Verifier使用和Prover一样的方法（原论文图8第4步）计算每一步的$\vec g，\vec a$。但实际上，只有Prover才需要每一次迭代的$\vec g, \vec a$，对于Verifier来说，只需要计算得到最后一步的$g,a$即可。因此在Halo中，这个过程被合并为两个内积。也即，$g=\vec s\cdot \vec g, a=\vec s\cdot \vec a$。具体实现可以参考代码，在这里简单的描述一下$\vec s$的构造，也即`Halo`论文中`Modified Inner Product`一节中的式子(1)。

1. 在每一轮中基于$FS$变换得到挑战标量$c$，合并这些挑战标量就构成一个向量$\vec c$。

2. 对$\vec c$中的每一个标量求倒数，得到向量$\vec d$。注意求一个向量的倒数有快速算法，并不需要对每个标量逐一求倒数。

3. 以$n=16$为例，一共迭代4次，每次递归过程中求$\vec g,\vec a$的过程如下：

   1. 一开始$\vec {g_0}$长度为$n$第一次递归：

      ```g1[0] = g0[0]* d0 + g0[8] * c0
      g1[1] = g0[1]*d[0]+g0[9]*c[0]
      g1[2] = g0[2]*d[0]+g0[10]*c[0]
      g1[3] = g0[3]*d[0]+g0[11]*c[0]
      g1[4] = g0[4]*d[0]+g0[12]*c[0]
      g1[5] = g0[5]*d[0]+g0[13]*c[0]
      g1[6] = g0[6]*d[0]+g0[14]*c[0]
      g1[7] = g0[7]*d[0]+g0[15]*c[0]
      ```

   2. 第二次迭代：

      ```
      g2[0] = g1[0]*d[1]+g1[4]*c[1]
      g2[1] = g1[1]*d[1]+g1[5]*c[1]
      g2[2] = g1[2]*d[1]+g1[6]*c[1]
      g2[3] = g1[3]*d[1]+g1[7]*c[1]
      ```

   3. 第三次迭代：

      ```
      g3[0] = g2[0]*d[2]+g2[2]*c[2]
      g3[1] = g2[1]*d[2]+g2[3]*c[2]
      ```

   4. 第四次迭代：

      ```
      g4 = g3[0]*d[3]+g3[1]*c[3]
      ```

   5. `g4`就是Verifier最后要求的$g$，展开并合并同类项有：

   ```
   g4 =
   g0[0]*d[0]*d[1]*d[2]*d[3]
   g0[1]*d[0]*d[1]*d[2]*c[3]
   g0[2]*d[0]*d[1]*c[2]*d[3]
   g0[3]*d[0]*d[1]*c[2]*c[3]
   g0[4]*d[0]*c[1]*d[2]*d[3]
   g0[5]*d[0]*c[1]*d[2]*c[3]
   g0[6]*d[0]*c[1]*c[2]*d[3]
   g0[7]*d[0]*c[1]*c[2]*c[3]
   g0[8]*c[0]*d[1]*d[2]*d[3]
   g0[9]*c[0]*d[1]*d[2]*c[3]
   g0[10]*c[0]*d[1]*c[2]*d[3]
   g0[11]*c[0]*d[1]*c[2]*c[3]
   g0[12]*c[0]*c[1]*d[2]*d[3]
   g0[13]*c[0]*c[1]*d[2]*c[3]
   g0[14]*c[0]*c[1]*c[2]*d[3]
   g0[15]*c[0]*c[1]*c[2]*c[3]
   ```

4. 可以注意到$s$d的每个元素正好和其所在位置的2进制位有关。于是求`s[k]`有一个简单的做法如下：

   ```
   Fr ret = FrOne();
   std::bitset<64> bits(k);
   for (size_t i = 0; i < round; ++i) {
     if (bits[round - i - 1]) {
       ret *= c[i];
     } else {
       ret *= d[i];
     }
   }
   return ret;
   ```

   

5. 以上方法所需要的有限域计算太多，因为每一项都能找到一项和它仅相差一个bit，因此可以利用此前计算的结果。

   ```
       for (size_t k = 0; k < s.size(); ++k) {
         if (k == 0) {
           s[k] = FrOne();
           for (size_t i = 0; i < round; ++i) {
             s[k] *= d[i];
           }
         } else {
           std::bitset<64> bits(k);
           for (size_t i = 0; i < round; ++i) {
             if (bits[round - i - 1]) {
               bits[round - i - 1] = 0;
               auto kk = bits.to_ulong();
               s[k] = s[kk] * cc[i];
               break;
             }
           }
         }
       }
   
   ```

   

6. 这个优化使得Verifier仅需要做一次$MulExp(n)$，而不需要像Prover一样，反复折叠，以至于无法充分利用$MulExp$所带来的优势。此外，有限域的计算量相比A3原版也大大降低了。

##### 性能

- 证明开销： $MulExp(n) \sim 3MulExp(n)$

- 验证开销： $MulExp(n )$
- 证据长度：$(2+2\lg(n))G1+2F_r$

