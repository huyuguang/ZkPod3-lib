# Bulletproofs

https://eprint.iacr.org/2017/1066.pdf

### 概要
本文仅包含论文中的protocol 1和protocol 2，以及基于protocol 1证明内积（在文中被称之为protocol 3）。本文不关心和电路有关的部分。简而言之，本文最终要解决的问题如下：  
Prover有两个Fr中的秘密向量a、b，prover公开com(a)、com(b)、com(c)并且证明c=<a,b>。  
- h、u是公开的独立的G生成元，g<sub>1</sub>、g<sub>2</sub>是公开的独立的G生成元矢量
- com(a)=h<sup>r<sub>1</sub></sup>g<sub>1</sub><sup>a</sup>
- com(b)=h<sup>r<sub>2</sub></sup>g<sub>2</sub><sup>b</sup>
- com(c)=h<sup>r<sub>3</sub></sup>u<sup>c</sup>  

### 术语约定
G是一条椭圆曲线（y<sup>2</sup>=ax<sup>3</sup>+bx+c）上的点构建的群。F<sub>p</sub>是x和y的定义域，F<sub>r</sub>是该群的阶构造的有限域。  

g、h、u是独立的G生成元，独立生成元的含义是任何人都不知道它们之间的指数关系。例如无人知道方程g<sub>1</sub>=g<sub>2</sub><sup>x</sup>中的x是多少。某些时候为简化起见，会用g表示一个生成元向量，也会用g<sub>1</sub>表示一个生成元向量，用g<sub>2</sub>表示另一个生成元向量。当g表示一个向量时，g<sup>a</sup>表示g[0]<sup>a[0]</sup>g[1]<sup>a[1]</sup>......。  

通过如下方式可以产生满足这种条件的g和h：基于某个hash函数，例如sha256反复hash得到一组随机字节。然后将这些字节映射到G上的某个点，由此得到的点就是g和h。  

可以这样将字节映射到G上的点：首先将随机字节映射到一个F<sub>r</sub>，然后以该F<sub>r</sub>作为x，求解椭圆曲线方程y<sup>2</sup>=ax<sup>3</sup>+bx+c。在这里需要进行一次有限域上的开方运算，最后得到y。有较小概率ax<sup>3</sup>+bx+c得到的数不是一个平方数，这种情况下反复++x，直到是平方数为止。

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
- h、u是公开的独立的G生成元，g<sub>1</sub>、g<sub>2</sub>是公开的独立的G生成元矢量
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