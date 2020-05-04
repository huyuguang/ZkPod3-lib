# Bulletproofs

https://eprint.iacr.org/2017/1066.pdf


### 概要

本文不关心和电路有关的部分。本文仅包含论文中的protocol 1和protocol 2，以及基于protocol 1证明内积（在文中被称之为protocol 3.1，3.2，3.3）。其中Protocol 3.1并未直接出现在论文中，Protocol 3.2是基于3.1的一个小变化，Protocol 3.3是对多个内积证明的聚合优化。

[基础概念和符号约定](../ecc/readme.md)

-------

### Protocol 2

先介绍Protocol 2而不是Protocol 1是因为Protocol 1基于Protocol 2。

##### 语义

Prover有两个$F_r$中的秘密向量$\vec a,\vec b $，公开$P_2=\vec{g_1}^{\vec a}\vec{g_2}^{\vec b}u^{\vec a \cdot \vec b}$，证明自己知道被$P_2$锁定的$\vec a,\vec b $。

##### 协议

算法细节请参考论文16页。

-------

### Protocol 1

##### 语义

Prover有两个$F_r$中的秘密向量$\vec a,\vec b $，公开$P_1=\vec {g_1}^{\vec a}\vec{g_2}^{\vec b}, c=\vec a \cdot \vec b$，证明自己知道被$P_1$锁定的$\vec a,\vec b$，且$c=\vec{a} \cdot \vec b$。

##### 协议

1. Verifier产生一个随机挑战值$x$，令$U=u^x$，或者直接使用$FS$变换产生$x=FS(P1,C,n)$。

2. Prover和Verifier同时计算$P_2=\vec{g_1}^{\vec a} \vec{g_2}^{\vec b} u^{xc}=P_1 u^{xc}=P_1 U^c$。
3. Prover运行Protocol 2证明。

算法细节请参考论文15页。注意Protocol 1和Protocol 2是soundness的，但不是zero knownledge的。

------

### Protocol 3.1

Protocol 3.1并未直接出现在论文中，以下内容根据论文第4章整理得来。

##### 语义

- $\vec a, \vec b$是两个长度为$n$的$F_r $秘密向量。
- Prover公开$com(\vec a,\vec b), com(c)$，证明$c=\vec a \cdot \vec b$
- $P=com(\vec a,\vec b)=h^\alpha {\vec {g_1}}^{\vec a} \vec{g_2}^\vec{b}$
- $Q=com(c)=h^\beta u^c$
- $\alpha, \beta$是两个随机$F_r$标量
- 注意$\vec{g_1},\vec{g_2}$不能重叠

##### 协议

1. Prover选择长度为$n$的$F_r$随机向量$\vec {d_a},\vec {d_b}$，$F_r$随机标量$\rho,\tau_1,\tau_2$。
2. 计算$r=h^\rho \vec {g_1}^{\vec {d_a}} \vec{g_2}^{\vec {d_b}}$。
3. 计算$t_1=h^{\tau_1}u^{\vec a \cdot \vec {d_b} + \vec b \cdot \vec{d_a}}, t_2=h^{\tau_2}u^{\vec {d_a} \cdot \vec{d_b}}$。
4. 使用$FS$变换得到挑战$x=FS(P,Q,r,t_1,t_2,n)$
5. 计算$\vec {a_2}=\vec a+\vec {d_a} x, \vec{b_2}=\vec b+ \vec{d_b} x, c_2=\vec {a_2} \cdot \vec {b_2}$。
6. 计算$\mu=\alpha+x \rho, \tau=\beta + x {\tau_1} + x^2{\tau_2} $。
7. 计算$P_1=Pr^xh^{-\mu}$。
8. 用$P_1,c_2,\vec{a_2},\vec{b_2}$做Protocol1证明，得到$proof1$。
9. 公开$Proof=c_2,\mu,\tau,r,t_1,t_2,c_2,proof1$。
10. Verifier验证$Q t_1 ^x t_2 ^ {x^2} \overset ? = h^\tau u^{c_2}$。
11. Verifier计算$P_1=P r^xh^{-\mu}$。
12. Verifier验证$proof1$。

##### 性能

- 证明开销：
- 验证开销：
- 证据长度：$(4+2\lg(n))G1+6F_r$

----------

### Protocol 3.2

Protocol 3.2是一个更接近实际使用的Protocol 3.1的变种。其语义为Prover有两个$F_r$中的秘密向量$\vec a,\vec b$，Prover公开$com(\vec{a}),com(\vec b),com(c)$并且证明$c=\vec a \cdot \vec b $。

- $h,u$是公开的独立的$G1$生成元，$\vec{g_1},\vec{g_2}$是公开的独立的$G1$生成元向量
- $A=com(\vec a)=h^{r_1}{\vec{g_1}}^\vec a$
- $B=com(\vec b)=h^{r_2}{\vec {g_2}}^\vec b$
- $C=com(c)=h^{r_3}u^c$

很显然很容易将$com(\vec a),com(\vec b)$转为Protocol 3所需要的输入，也即$P=h^\alpha \vec{g_1}^\vec a \vec{g_2}^\vec b$ 以及 $Q=h^\beta u^c$。令$\alpha=r_1+r_2$，那么有$P=com(\vec a)com(\vec b)=h^\alpha \vec {g_1}^\vec a \vec {g_2}^\vec b$，同样令$\beta=r_3$，那么有$Q=com(c)$。

但是有时候我们的$com(a), com(b)$会使用相同的向量$\vec {g_1}$（或者$\vec {g_1},\vec{g_2}$存在部分项相同甚至完全相同），也即：

- $com(a)=h^{r_1}\vec{g_1}^\vec a$
- $com(b)=h^{r_2}\vec{g_1}^\vec b$
- $com(c)=h^{r_3}u^c$

这种情况下，首先需要把$com(\vec b)$换成以$\vec{g_2}$为底，也即，公开$B'=com(b')=h^{r'}\vec{g_2}^\vec{b'}$，并且证明对于公开的$com(\vec{b'}),com(\vec b)$，有$\vec{b'}=\vec b$。具体做法是由Verifier挑战一个随机值$x$，或者用$FS$变换直接生成$x=FS(A,B,C,B')$。然后基于$x$得到一个随机向量$\vec v$，Prover证明$\vec b \cdot \vec v = \vec{b'} \cdot \vec v$（例如，使用Hyrax A2或者A3，参见hyrax目录下的readme.md）。如果$com(c)$所使用的$u$和$com(\vec a)$所使用的$\vec{g_1}$也存在重叠，例如$u$和$\vec{g_1}$中的某个元素相等，那么也可以用类似的方法进行换底。

--------

### Protocol 3.3

有些情况下需要证明一个向量和多个向量的内积，且不希望重复多次Protocol 3.2。

##### 语义

- $\vec{a}$是$F_r$中的秘密向量。
- $\vec{b}_i, i\in [1,n]$是$F_r$中的n个秘密向量。
- $c_i, i\in[1,n]$是n个$F_r$中的秘密标量。
- Prover公开$com(\vec{a}),com(\vec{b}_i), com(\vec c),i\in[1,n]$，证明$c_i=\vec{a}_i \cdot \vec{b}_i,i\in[1,n]$。
- $A=com(\vec{a})=h^{r_{a}}\vec{g}^{\vec{a}}$。
- $B_i=com(\vec{b}_i)=h^{r_{bi}}\vec{g}^{\vec{b}_i}, i\in[1,n]$。
- $C=com(\vec c)=h^{r_c }\vec g^{\vec c}$。

##### 协议

1. 使用$FS$变换得到挑战$x=FS(A,\vec B,C)$。
2. 令$\vec x=\{x^0,x^1,x^2...\}$。
3. 令$\vec{d}_i=\vec{b}_ix_i, D_i=com(r_{di},\vec d_i)$，显然有$D_i=com(r_{bi},\vec b_i)^{x_i}, r_{di}=r_{bi}^{x_i} $。
5. 令$ \vec{e}=\Sigma_{i=1}^{n}\vec{d}_i, E=com(\vec{e})$，显然有$$。
5. 显然如下$n$个等式成立：$c_ix^i=\vec{a}\cdot (\vec{b}_ix^i)$，于是可以把它们合并为$\vec c \cdot \vec x=\Sigma_{i=1}^{n}c_ix^i=\Sigma_{i=1}^{n}\vec{a}\cdot(\vec{b}_ix^i)=\vec{a}\cdot (\Sigma_{i=1}^{n}\vec{b}_ix^i)=\vec{a}\cdot \vec{e}$，令等式左边为$L$，等式右边为$R$。
6. 通过HyraxA3证明$com(L)$。
7. 证明com(R)：
   - 通过线性组合得到$E$。
   - 通过Protocol 3.2证明$R=\vec{a}\cdot \vec{e}$。
8. Verifier验证4,5的两个证明，并且验证$com(L)\overset{?}=com(R)$。



