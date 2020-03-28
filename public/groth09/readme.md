# Groth09

http://www0.cs.ucl.ac.uk/staff/J.Groth/MatrixZK.pdf

### 概要
本文仅包含论文中的4.3、5.1、5.2、5.3中涉及到的协议或算法。本文不关心原论文中的其他部分。在实现论文中涉及到的算法时，进行了一定的扩展和优化。论文中的算法次序发展次序是5.1->5.2->5.3->4.3，因此在本文中也按照这个次序。

论文中阐述的协议都是交互式的，很显然这是不方便使用的。因此在实现时，都通过Fiat-Shamir-Transform转为了非交互式证明，以下简称为$FS$变换。

我非常希望除了能详细的介绍协议细节以及实现细节，还能介绍安全性证明以及动机，也即为什么论文作者要这么设计。但是目前似乎无法做到，一方面要讲清楚动机是个很困难的事情，需要对问题有充分的理解，而我现在还达不到这个水平；另一方面也缺乏时间和精力。因此先用这样一个偏代码实现的角度介绍一下源码中涉及到的协议。希望以后能有时间和精力再来补齐这个留白。

[基础概念和符号约定](../ecc/readme.md)



### 5.1b (5.1a)
##### 语义

- Prover公开$com(\vec{x})、com(\vec{y})、com(z)$并证明$z=\vec{x}*\vec{y}$。

- $\vec t$是$F_r$中的1个公开向量，长度为$n$。
- $\vec x、\vec y$是$F_r$中2个秘密向量，长度均为$n$。
- $z$是$F_r$中的1个秘密标量。
- $\vec{x}*\vec{y}$表示的是一个双线性映射(bilinear map)，例如$\vec{x}\cdot\vec{y}$或者$\vec x\cdot (\vec y\circ \vec t)$，其中$\vec x\cdot \vec y$表示两个向量的内积，而$\vec y\circ \vec t$表示两个向量的Hadamard Product。
- 实现时出于性能和代码清晰的考虑，把$\vec{x}*\vec{y}$ 代表的$\vec x\cdot \vec y$和$\vec x\cdot (\vec y\circ \vec t)$两种语义分开成了5.1a和5.1b。5.1a只是5.1b的一个特例，很明显当5.1b中的$\vec t$每个元素都为1时，5.1b就退化成了5.1a。
- $com(\vec x)=h^{r}\vec{g_{x}}^\vec{x}$。
- $com(\vec y)=h^{s}\vec{g_y}^\vec{y}$。
- $com(z)=h^{t}g_z^z$。
- $r、s、t$都是$F_r$中的随机秘密标量。
- $\vec{g_x}、\vec{g_y}$是$G$中公开的独立向量，长度为$n$。
- $g_z$是G中公开的独立标量。
- $\vec{g_x}^{\vec x}$指$\prod_{i=0}^{n-1} g_{xi}^{x_i}$，这是计算开销最大的部分，通常也称之为$multiexp (n)$。
- 原论文中$com(\vec x)、com(\vec y)$使用了相同的$\vec g$，但实际上$\vec{g_x}、\vec{g_y}$可以互不相同。

##### 证明

?	证明方法完全遵循原论文5.1节。不同之处在于用$FS$变换把协议转为了非交互的。

1. Prover选择3个$F_r$随机数$r,s,t$。
2. 计算$a=com(\vec{x}),b=com(\vec{y}),c=com(z)$。
3. Prover 选择长度为$n$的随机$F_r$向量$\vec {d_x}, \vec {d_y}$，随机$F_r$标量$d_z$，随机$F_r$标量$r_d,s_d,t_1,t_0$。以上随机数都是Prover的秘密。
4. Prover锁定这些随机数，也即计算$a_d,b_d,c_1,c_0$：
   - $a_d=h^{r_d}\vec{g_x}^\vec{d_x}$
   - $b_d=h^{s_d}\vec{g_y}^\vec{d_y}$
   - $c1=h^{t_1}g_z^{\vec{x}\cdot({\vec{d_y}\circ{\vec{t}})}+\vec{d_x}\cdot(\vec{y}\circ{\vec{t}})}$，注意在这里论文中的$\vec{x}*\vec{d_y}$和$\vec{d_x}*\vec{y}$实现成了$\vec{x}\cdot (\vec{d_y}\circ{\vec{t}})$以及$\vec{d_x}\cdot (\vec{y}\circ{\vec{t}})$。在5.1a中，没有这个$\vec{t}$，因此就退化成了$\vec{x}\cdot{\vec{d_y}}$。
   - $c_0=h^{t_0}{g_z}^{\vec{d_x}\cdot{(\vec{d_y}\circ{\vec{t}})}}$。注意在这里论文中的$\vec{d_x}*\vec{d_y}$实现成了$\vec{d_x}\cdot{(\vec{d_y}\circ{\vec{t}})}$，在5.1a中，退化成了$\vec{d_x}\cdot{\vec{d_y}}$。
   - 计算$a_d, b_d$时有两次$multiexp(n)$。
5. 原论文中Prover发送$a_d,b_d,c_1,c_0$给Verifier，然后Verifier挑战一个随机数$e$。在这里我们用$FS$变换把协议转为非交互的。$e=sha256(seed, a,b,c,a_d,b_d,c_1,c_0)$。其中$seed$是Prover和Verifier交互时双方约定的一个随机数，类似于session cookie这样的角色。

5. Prover计算并发送$Proof$。
   - $\vec{f_x}=e\vec{x}+\vec{d_x}$
   - $\vec{f_y}=e\vec{y}+\vec{d_y}$
   - $r_x=er+{r_d}$ 注意这个$r$是$com(\vec{x})$中的$r$
   - $s_y=es+s_d$ 注意这个$s$是$com(\vec{y})$中的$s$
   - $t_z=e^2t+e{t_1}+t_0$ 注意这个$t$是$com(z)$中的$t$
   - $Proof={a,b,c,a_d,b_d,c_1,c_0,\vec{f_x},\vec{f_y},r_x,s_y,t_z}$

6. Verifier验证Proof。
   - 使用同样的方法计算挑战随机数$e$
   - 验证如下三个等式是否成立，如果都成立则验证通过。
     - $a^e{a_d}?=h^{r_x}\vec{g_x}^\vec{f_x}$
     - $b^e{b_d}?=h^{s_y}\vec{g_y}^\vec{f_y}$
     - $e^{e^2}{c_1}^e{c_0}?=h^{t_z}{g_z}^{f_z}, f_z=\vec{f_x}\cdot(\vec{f_y}\circ{\vec{t}})$，在5.1a中没有$\vec{t}$，因此退化为$\vec{f_x}\cdot{\vec{f_y}}$
     - 第一个等式和第二个等式有两次$multiexp(n )$

##### 性能

- 证明开销和证明长度不包括$a=com(\vec{x}), b=com(\vec{y}),c=com(z)$部分。
- 证明开销：2次$multiexp(n)$。
- 验证开销：2次$multiexp(n)$。
- 证明长度：$4G1+(2N+3)F_r$。
- 可以看出，当退化为5.1a时，区别只是$O(n)$级别的$F_r$中乘法计算（几个$F_r$向量内积），因此可以认为5.1a和5.1b的性能是一样的。

### 5.1c
##### 语义

- 和5.1b几乎完全相同。
- $g_z\notin \vec{g_x}$ 。这个要求是由证明中所用到的Bullet Proof P3引入的。通常计算$com(z)$时可以随意的选择$g_z $，但如果$com(z)$是从上一级证明中引入的（例如5.1c是多个串联证明中的一环，也即5.1c的输入是另一个证明的输出），则可能不满足这个条件。那么也很容易用一个小变换将$g_z$转换为$g_{z'}$同时将$com(z)$转换为$com(z')$并证明$z=z'$：
  - 令$z'=z$，重新选择一个$g_{z'}$以及随机$F_r $标量$t'$并且计算得到$com(z)' =h^{t'}g_{z'}^{z'}$ 。
  - 基于$com(z), com(z')$使用$FS$变换得到挑战数$c$。
  - 使用两次Hyrax A2证明$zc==z'c$，注意Hyrax A2的语义是证明向量内积，此时$z,c$都是标量，也即只有一个元素的向量。详细细节参看clink目录下的readme.md。
- 原始论文（groth09）中并不存在该算法，该算法是本文概要中所提到的扩展的一部分。注意到原论文中5.1所述的算法的证明长度都是$O(N)$，我们希望能牺牲部分性能换取证明长度降为$O(logN)$。
- 该算法使用HyraxA3和Bullet Proof P3组合而来，和原始论文中的5.1毫无关系。

##### 证明

可以注意到5.1a的语义非常类似于[Bullet Proof P3](../bp/readme.md)，但存在一定的差别，主要体现在commitment时所用的$\vec g$。和5.1b的区别则是那个额外的$\vec t$，则可以通过Hyrax A3来解决。Bullet Proof发明了一种递归折叠的技术（Hyrax A3也是借鉴了Bullet Proof的折叠方法），用更大的证明开销和更大的验证开销换来了证明长度从$O(N)$降到$O(logN)$，在某些情况下这是有利的。大致思路如下：

  - 公开$com(\vec{\underline{yt}})$，并且通过Hyrax A3证明$\vec{\underline{yt}}=\vec y\circ \vec t$。
  - 用Bullet Proof Protocol3证明$z=\vec x\cdot \vec{\underline{yt}}$。

以下是详细步骤：

1. Prover选择3个随机数$r,s,t\in F_r$。
2. 计算$a=com(\vec{x}),b=com(\vec{y}),c=com(z)$。
3. 定义并计算$\underline{\vec{yt}}=\vec{y}\circ{\vec{t}}$，选择随机数$r_{yt}\in{F_r}$，计算$com(\underline{\vec{yt}})=h^{r_{\underline{yt} }}\vec{g_{\underline{yt} }}^{\underline{\vec{yt}}}$，其中$\vec{g_{yt}}$是$G$中长度为n的公开向量。
4. 为了满足Bullet Proof Protocol3 的要求，上式的$\vec{g_{\underline{yt}}}$不能和$\vec{g_x}$重叠。
5. 用HyraxA3证明对于被$com(\underline{\vec{yt}})$锁定的$\underline{\vec{yt}}$，满足$\underline{\vec{yt}}=\vec y\circ \vec t$：
   1. 用$FS$变换基于刚计算出来的$com(\underline{\vec{yt}})$产生一个长度为$n$的$F_r$向量$\vec{e}$作为挑战向量。
   2. 证明$\underline{\vec{yt}}\cdot\vec{e} == \vec{y}\cdot({\vec{e}\circ\vec{t}})$。
6. 现在有$com(\underline{\vec{yt}}), com(\vec{x}), com(z)$，需要证明$z==\vec{x}\cdot{\underline{\vec{yt}}}$：
   1. 注意到这非常接近Bullet Proof P3的语义。
   2. Bullet Proof要求内积的两个向量由不同的$\vec{g}$锁定，而$\vec{g_{\underline{yt}}}$和$\vec{g_x}$不重叠正好满足此条件。
   3. $com(\underline{\vec{yt}})* com(\vec{x})=h^{r_{\underline{ yt} }}\vec{g_{\underline{yt} }}^{\underline{\vec{yt}}}h^{r}\vec{g_{x}}^\vec{x}=h^{r+r_{\underline{ yt} }}\vec{g_{\underline{yt} }}^{\underline{\vec{yt}}}\vec{g_{x}}^\vec{x}$。
   4. 回顾一下Bullet Proof P3的语义
      1. 公开$P=com(\vec{a},\vec{b})=h^\alpha\vec{g_a}^\vec{a}\vec{g_b}^{\vec{b}}$
      2. 公开$Q=com(c)=h^{\beta}g_c^c$
      3. 证明$c=\vec{a}\cdot\vec{b}$
   5. 转为Bullet Proof P3:
      1. $\alpha=r+r_{\underline{yt}}$
      2. $\beta=t$
      3. $P=com(\underline{\vec{yt}})* com(\vec{x})$
      4. $Q=com(z)$
7. Prover构造Proof并公开，Proof由三部分组成：
   1. $com(\underline{\vec{yt}})$
   2. 步骤5产生的proof1
   3. 步骤6产生的proof2
8. Verifier验证Proof
   1. 使用相同的方法通过$FS$变换得到挑战随机数。
   2. 验证proof1，也即确认$\underline{\vec{yt}}=\vec y\circ \vec t$。
   3. 验证proof2，也即确认$z==\vec{x}\cdot{\underline{\vec{yt}}}$。

##### 性能

- 证明开销和证明长度不包括$a=com(\vec{x}), b=com(\vec{y}),c=com(z)$部分。
- 可以看出证明开销为两个子证明的开销之和，其中HyraxA3部分的开销是$3Multiexp(n)$
- 可以看出验证开销为两个子证明的开销之和，其中HyraxA3部分的开销是$2nEccexp$
- 证明长度为$4G1+2log(n)F_r$。

### 5.2b(5.2a)

##### 语义

- Prover公布$com(\vec{x_i}), com(\vec{y_i}), com(z)$并证明$z=\sum_{i=0}^{m-1}\vec{x_i}*\vec{y_i}$。
- $x、y$是$F_r$中2个秘密矩阵，大小均为$m$行$n$列。
- $\vec{x_i}$表示矩阵$x$第$i$行，$\vec{y_i}$表示矩阵$y$第$i$行。
- $z$是$F_r$中的秘密标量。
- $a=com(\vec {x_i})=h^{r_i}\vec{g_x}^{\vec x_i}, i\in[0,m-1]$ 。
- $b=com(\vec {y_i})=h^{s_i}\vec{g_y}^{\vec y_i}, i\in[0,m-1]$ 。
- $c=com(z)=h^{t }g_z^z$。
- $\vec{x}*\vec{y}$表示的是一个双线性映射(bilinear map)，例如$\vec{x}\cdot\vec{y}$或者$\vec x\cdot (\vec y\circ \vec t)$，其中$\vec x\cdot \vec y$表示两个向量的内积，而$\vec y\circ \vec t$表示两个向量的Hadamard Product。
- 实现时出于性能和代码清晰的考虑，把$\vec{x}*\vec{y}$ 代表的$\vec x\cdot \vec y$和$\vec x\cdot (\vec y\circ \vec t)$两种语义分开成了5.2a和5.2b。5.2a只是5.2b的一个特例，很明显当5.2b中的$\vec t$每个元素都为1时，5.2b就退化成了5.2a。
- $\vec{r}、\vec{s}$都是$F_r$中长度为$m$的随机秘密$F_r$向量，$r_i, s_i$指第$i$个$F_r$标量。
- $t$是秘密$F_r$标量。
- $\vec{g_x}、\vec{g_y}$是$G$中公开的独立向量，长度为$n$。
- $g_z$是G中公开的独立标量。
- $\vec{g_x}^{\vec{x_i}}$指$\prod_{j=0}^{n-1} g_x^{x_{ij}}$，这是计算开销最大的部分，通常也称之为$multiexp (n)$。
- 原论文中$com(\vec x)、com(\vec y)$使用了相同的$\vec g$，但实际上$\vec{g_x}、\vec{g_y}$可以互不相同。

##### 证明

证明方法完全遵循原论文5.2节。不同之处在于用$FS$变换把协议转为了非交互的。

详细步骤略。

##### 性能

- 算法特点是简洁但有O($m^2n$)级别的$F_r$运算，因此性能不佳，主要目的是为了引出随后的5.3方案。
- 最终将问题转换为5.1a或5.1b。

### 5.3b(5.3a)

##### 语义

- 语义和5.2b或者5.2a完全相同。
- 5.2b(5.2a)的递归版本，通过增加交互次数，降低了5.2a中$F_r$中的计算量。
- 最终将问题转换为5.1c或5.1b或5.1a。也即，根据模板参数的不同，5.3b具有不同的证明、验证开销以及证明长度。

##### 证明

证明方法完全遵循原论文5.2节。不同之处在于用$FS$变换把协议转为了非交互的。注意5.3是个递归算法，因此有多次$FS$变换，同时需要对齐数据。

以下是详细步骤：

1. Prover选择一个$Fr_r$随机数$t$，选择两个长度为$m$的$F_r$随机向量$\vec r, \vec s$。

2. 计算$a_i=com(\vec{x_i}),b_i=com(\vec{y_i}),c=com(z), i\in[0,m-1]$。$\vec{a},\vec{b}$都是长度为$m$的$G1$向量，$c$是$G1$标量。

3. 因为5.3是一个递归折叠算法，因此在这里需要对齐数据，使得矩阵$x,y$的行数从$m$变成2的幂，不够的填充垃圾数据。令$m'\geq m$且$m'$是2的幂。如果$m'\geq{m}$，则需要补齐：
   - $x,y$用$F_r(0)$补齐为$m'$行$n$列的矩阵。
   - $\vec{a},\vec{b}$用$G1(1)$补齐为长度为$m'$的向量。
   - $\vec{r},\vec{s}$用$F_r(0)$补齐为长度为$m'$的向量。
   - 毫无疑问，这种补齐方法并不会破坏$com()$的计算方法，也即，既然$\vec{x_m}$里面所有的元素都是$F_r(0)$显然有$G1(1)=h^{F_r(0)}\vec{g_x}^{\vec{x_{m}}}$。
   - 注意代码中使用了椭圆曲线的加法群表示法而不是本文中使用的乘法群表示法，因此代码中的零元不是$G1(1)$而是$G1(0) $。
   
4. 递归折叠
   1. $\underline{xy1}=\sum_{i=0}^{m/2-1}{\vec{x_{2i+1}}*\vec{y_{2i}}}$。
   2. $\underline{xy2}=\sum_{i=0}^{m/2-1}{\vec{x_{2i}}*\vec{y_{2i+1}}}$。
   3. 选择两个$F_r$随机标量$t_l, t_u$。
   4. 计算$c_l=h^{t_l}g_z^\underline{xy1}, c_u=h^{t_u}g_z^\underline{xy2}$。
   5. 基于$\vec{a},\vec{b},c,c_l,c_u$进行$FS$变换，得到挑战值e。
   6. 更新$x,y,z$：
      1. $\vec{x_i}=\vec{x_{2i+1}}{e}+\vec{x_{2i}}, \vec{y_i}=\vec{y_{2i}}{e}+\vec{y_{2i+1}}, i\in[0,m/2-1]$。
      2. $z=e^2\underline{xy1}+ez+\underline{xy2}$。
      3. 注意此时$x,y$成了$m/2$行$n$列。
   7. 更新$\vec a,\vec b,c,\vec r,\vec s,t$，也即$com(x),com(y),com(z)$：
      1. $\vec{a_i}=\vec{a_{2i}}+\vec{a_{2i+1}}e,\vec{b_i}=\vec{b_{2i}}e+\vec{b_{2i+1}},r_i=r_{2i}+r_{2i+1}e,s_i=s_{2i}e+s_{2i+1},i\in[0,m/2-1]$。
      2. $c=c_le^2+ce+c_u, t=t_le^2+te+t_u。$
      3. 注意$\vec{a},\vec{b}$此时成了长度为$m/2$的向量。
   8. 新的$x,y,z,\vec{a},\vec{b},c,\vec{r},\vec{s},t$是自洽的。
   9. 回到步骤1，直到$x,y$变成了1行$n$的矩阵，也即，一个长度为$n$的向量。同样，此时$\vec{a},\vec{b},\vec{r},\vec{s}$成了长度为1的向量，也即，成了标量。
   10. 重新定义$\vec{x}=x_0, \vec{y}=y_0,z=z,a=a_0,b=b_0,c=c,r=r_0,s=s_0,t=t$。
   
5. 基于当前的$\vec{x},\vec{y},z,a,b,c,r,s,t$，使用5.1b或者5.1c证明。

6. Prover构造Proof并公开，Proof由两部分组成:

   - 递归折叠中每次计算得到的$c_l,c_u$

   - 5.1b或者5.1c的证明。

7. Verifier验证Proof：

   1. 用同样的折叠方法计算得到最后的$a,b,c$
   2. 用5.1b或者5.1c验证。

##### 性能

递归折叠之后使用5.1证明时可以选择5.1b和5.1c，两种方案的性能特点不同。前者证明和验证开销小而证明长度大，后者证明和验证开销大而证明长度小。

- 5.1b方案 TODO
- 5.1c方案 TODO



### 4.3b

##### 语义

- Prover公布$com(\vec{x_i}), com(\vec{y_i}), com(\vec{z_i})$并且证明$z=x \circ y$
- $x、y、z$是$F_r$中3个秘密矩阵，大小均为$m$行$n$列。
- $com(\vec{x_i})=h^{r_i}\vec{g_x}^\vec{x_i}, i\in[0,m-1] $。
- $com(\vec{y_i})=h^{s_i}\vec{g_y}^\vec{y_i}, i\in[0,m-1]$。
- $com(\vec{z_i})=h^{t_i}\vec{g_z}^\vec{z_i}, i\in[0,m-1]$。
- $x \circ y$指Hadamard Product，也即$z_{i,j}=x_{i,j}y_{i,j}, i\in[0,m-1],j\in[0,n-1]$。
- $\vec{r},\vec{s},\vec{t}$ 是三个长度为$n$的随机$F_r$向量，$r_i,s_i,t_i$表示向量中第$i$个元素。
- 原论文中$com(\vec{x_i}), com(\vec{y_i}), com(\vec{z_i})$使用了相同的$\vec g$，但实际上$\vec{g_x}, \vec{g_y}, \vec{g_z}$可以互不相同。

##### 证明

证明方法大体遵循原论文，也即将问题转换为5.3，但在具体实现时略微做了调整，最终将问题转换为5.3和Hyrax A2/3。也即，根据模板参数的不同，4.3b具有不同的证明、验证开销以及证明长度。原始论文中并未用到Hyrax A2/3，此处是对原始论文的一点小改进。当然，同样用$FS$变换把协议转为了非交互的。

以下是详细步骤：

1. Prover选择三个长度为$m$的$F_r$随机向量$\vec r, \vec s, \vec t$。

2. 计算$a_i=com(\vec{x_i}),b_i=com(\vec{y_i}),c_i=com(\vec{z_i}), i\in[0,m-1]$。$\vec{a},\vec{b,\vec{c}}$都是长度为$m$的$G1$向量。

3. 因为最终会将问题转换为5.3的语义，而5.3是一个递归折叠算法，因此在这里需要对齐数据，使得矩阵$x,y,z$的行数从$m$变成2的幂，不够的填充垃圾数据。令$m'\geq m$且$m'$是2的幂。如果$m'\geq{m}$，则需要补齐：

   - $x,y,z$用$F_r(0)$补齐为$m'$行$n$列的矩阵。
   - $\vec{a},\vec{b},\vec{c}$用$G1(1)$补齐为长度为$m'$的向量。
   - $\vec{r},\vec{s},\vec{t}$用$F_r(0)$补齐为长度为$m'$的向量。
   - 毫无疑问，这种补齐方法并不会破坏$com()$的计算方法，也即，既然$\vec{x_m}$里面所有的元素都是$F_r(0)$显然有$G1(1)=h^{F_r(0)}\vec{g_x}^{\vec{x_{m}}}$。
   - 注意代码中使用了椭圆曲线的加法群表示法而不是本文中使用的乘法群表示法，因此代码中的零元不是$G1(1)$而是$G1(0) $。

4. 基于$\vec a,\vec b,\vec c,m,n$进行$FS$变换得到两个$F_r$挑战向量：

   1. 长度为$m$的$\vec k$
   2. 长度为$n$的$\vec t $

5. 用$\vec k,\vec t$对$x,y$进行双向扰动：$Z=\sum_{i=0}^{m-1}(\vec{x_i}{k_i})\cdot (\vec{y_i}\circ \vec t))$，我们将从两个方向去证明它。  

6. 用5.3证明$Z=\sum_{i=0}^{m-1}(\vec{x_i}{k_i})\cdot (\vec{y_i}\circ \vec t))$  :

   1. 令$\underline{a53}_i=com(\vec{x_i}k_i)$，可以直接通过已有的$\vec a,\vec r$计算得到$\underline{a53}_i, \underline{r53}_i$：
      1. $\underline{a53}_i={a_i}^{k_i}$
      2. $\underline{r53}_i={r_i}{k_i}$
      3. 很显然有$\underline{a53}_i={a_i}^{k_i}={(h^{r_i}\vec{g_x}^{\vec x_i})}^{k_i}=h^{\underline{r53}_i}\vec{g_x}^{\vec{x_i}{k_i}}$，注意这个计算Verifier也可以进行，并不需要知道秘密。

   2. 选择$F_r$随机标量$\underline{t53}$并计算$\underline{c53}=com(Z)=h^{\underline{t53}}{g_z}^Z$。
   3. 使用5.3证明，得到$proof_{53}$。

7. 用HyraxA2或HyraxA3证明$Z=\vec{\underline{zk}}\cdot \vec t$ ，其中$\underline{\vec{zk}}$是一个长度为$n$的向量，且有 $\underline{zk}_j=\sum_{i=0}^{m-1}z_{ij}k_i, j\in[0,n-1]$：

   1. 可以直接通过已有的$\vec{c}$计算$\xi=com(\underline{zk})$:
      1. $\xi=\vec{c}^k=\prod_{i=0}^{m-1}c_i^{k_i}，$注意这个计算Verifier也可以进行，并不需要知道秘密。
      2. $r_{\xi}=\vec{t}\cdot{\vec k}$

   2. 复用前面所计算的$com(Z)$:
      1. $r_{\tau}=\underline{t53}$
      2. $\tau=\underline{c53}$
   3. 使用Hyrax A2或Hyrax A3证明，得到$proof_a$。

8. Prover构造Proof并公开，Proof由三部分组成:

   1. $\underline{c53}$，因为Verifier无法自行计算出这个值。
   2. $proof_{53}$
   3. $proof_a$

9. Verifier验证Proof：

   1. 基于$\vec a,\vec b,\vec c,m,n$进行$FS$变换得到两个$F_r$挑战向量$\vec{k}, \vec{t}$。
2. 用和Prover相同的方法计算得到$\vec{\underline{a53}}$。
   
   3. 验证5.3。
4. 用和Prover相同的方法计算得到$\xi$。
   5. 验证Hyrax A2/A3。

   

##### 性能

4.3b的实现会根据模板参数的不同选择不同的子证明。最主要的有如下两种：

- 5.1b+HyraxA2 TODO
- 5.1c+HyraxA3 TODO



