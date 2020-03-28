# Bulletproofs

https://eprint.iacr.org/2017/1066.pdf


### 概要
本文仅包含论文中的protocol 1和protocol 2，以及基于protocol 1证明内积（在文中被称之为protocol 3）。本文不关心和电路有关的部分。简而言之，本文最终要解决的问题如下：  
Prover有两个Fr中长度为n的秘密向量a、b，prover公开com(a)、com(b)、com(c)并且证明c=<a,b>。  
- h、u是公开的独立的G生成元，g<sub>1</sub>、g<sub>2</sub>是公开的独立的G生成元向量
- com(a)=h<sup>r<sub>1</sub></sup>g<sub>1</sub><sup>a</sup>
- com(b)=h<sup>r<sub>2</sub></sup>g<sub>2</sub><sup>b</sup>
- com(c)=h<sup>r<sub>3</sub></sup>u<sup>c</sup>  

[基础概念和符号约定](../ecc/readme.md)

### Protocol 2
先介绍Protocol 2而不是Protocol 1是因为Protocol 1基于Protocol 2。  
Prover有两个F<sub>r</sub>中的秘密向量a、b，公开P<sub>2</sub>=g<sub>1</sub><sup>a</sup>g<sub>2</sub><sup>b</sup>u<sup><a,b></sup>，证明知道a、b。  
算法细节请参考论文16页。  

### Protocol 1
- Prover公开P<sub>1</sub>=g<sub>1</sub><sup>a</sup>g<sub>2</sub><sup>b</sup>
- Prover公开c=<a,b>
- Verifier产生一个随即挑战x
- Prover和Verifier同时计算P<sub>2</sub>=g<sub>1</sub><sup>a</sup>g<sub>2</sub><sup>b</sup>u<sup>xc</sup>=P<sub>1</sub>u<sup>xc</sup>
- Prover运行Protocol 2证明
- 最终Prover证明了被P<sub>1</sub>锁定的向量a、b的内积是公开值c

算法细节请参考论文15页。   
注意Protocol 1和Protocol 2是soundness的，但不是zero knownledge的。

### Protocol 3
Protocol 3并未直接出现在论文中，以下内容根据论文第4章整理得来。  

Prover  
- 公开P=h<sup>*alpha*</sup>g<sub>1</sub><sup>a</sup>g<sub>2</sub><sup>b</sup> 以及 Q=h<sup>*beta*</sup>u<sup>c</sup>，其中*alpha*,*beta*是两个F<sub>r</sub>随机标量。  
- 选择F<sub>r</sub>随机向量d<sub>a</sub>、d<sub>b</sub>，F<sub>r</sub>随机标量*rho*、*tau<sub>1</sub>*、*tau<sub>2</sub>*
- 计算R=h<sup>*rho*</sup>g<sub>1</sub><sup>d<sub>a</sub></sup>g<sub>2</sub><sup>d<sub>b</sub></sup>
- 定义t<sub>1</sub>=<a,d<sub>b</sub>>+<b,d<sub>a</sub>>, t<sub>2</sub>=<d<sub>a</sub>,d<sub>b</sub>>
- 计算T<sub>1</sub>=h<sup>*tau<sub>1</sub>*</sup>u<sup>t<sub>1</sub></sup>，T<sub>2</sub>=h<sup>*tau<sub>2</sub>*</sup>u<sup>t<sub>2</sub></sup>
- 发布R，T<sub>1</sub>，T<sub>2</sub>

Verifier
- 挑战x

Prover
- 令a'=a+d<sub>a</sub>x，b'=b+d<sub>b</sub>x
- 计算c'=<a',b'>
- 计算*mu*=*alpha*+x*rho*，*tau*=*beta*+x*tau<sub>1</sub>*+x<sup>2</sup>*tau<sub>2</sub>*
- 发送c'，*mu*，*tau*

Verifier
- 检查是否满足QT<sub>1</sub><sup>x</sup>T<sub>2</sub><sup>x<sup>2</sup></sup> == h<sup>*tau*</sup>u<sup>c'</sup>

Prover
- 令**P<sub>1</sub>**=PR<sup>x</sup>h<sup>*-mu*</sup>，**c**=c'
- 令**a**=a'，**b**=b'
- 发送**P<sub>1</sub>**和**c**

随后Prover和Verify运行Protocol 1。  

### 最终问题
Prover有两个Fr中的秘密向量a、b，prover公开com(a)、com(b)、com(c)并且证明c=<a,b>。  
- h、u是公开的独立的G生成元，g<sub>1</sub>、g<sub>2</sub>是公开的独立的G生成元向量
- com(a)=h<sup>r<sub>1</sub></sup>g<sub>1</sub><sup>a</sup>
- com(b)=h<sup>r<sub>2</sub></sup>g<sub>2</sub><sup>b</sup>
- com(c)=h<sup>r<sub>3</sub></sup>u<sup>c</sup>  

很显然很容易将com(a)和com(b)转为Protocol 3所需要的输入，也即P=h<sup>*alpha*</sup>g<sub>1</sub><sup>a</sup>g<sub>2</sub><sup>b</sup> 以及 Q=h<sup>*beta*</sup>u<sup>c</sup>。令*alpha*=r<sub>1</sub>+r<sub>2</sub>，那么有P=com(a)com(b)=h<sup>*alpha*</sup>g<sub>1</sub><sup>a</sup>g<sub>2</sub><sup>b</sup>，同样令*beta*=r<sub>3</sub>，那么有Q=com(c)。  

但是有时候我们的com(a), com(b)会使用相同的向量g，也即：
- com(a)=h<sup>r<sub>1</sub></sup>g<sub>1</sub><sup>a</sup>
- com(b)=h<sup>r<sub>2</sub></sup>g<sub>1</sub><sup>b</sup>
- com(c)=h<sup>r<sub>3</sub></sup>u<sup>c</sup>  

这种情况下，首先需要把com(b)换成以g<sub>2</sub>为底，也即，公开com(b')=h<sup>r'</sup>g<sub>2</sub><sup>b'</sup>，并且证明对于公开的com(b')和com(b)，有b' == b。具体做法是由Verifier挑战一个随机值，以该随机值为种子得到一个随机向量v，Prover证明<b,v> == <b',v>（例如，使用Hyrax A2或者A3，参见hyrax目录下的readme.md）。  
如果com(c)所使用的u和com(a)所使用的g<sub>1</sub>也存在重叠，例如u == g<sub>1</sub>[0]，那么也可以用类似的方法进行换底。

### 性能
- Protocol 3
  - Prove cost：TODO
  - Verify cost：TODO
  - Proof size：(4+2log(n))G1+6F<sub>r</sub>
