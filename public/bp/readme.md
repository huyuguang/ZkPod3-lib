# Bulletproofs

https://eprint.iacr.org/2017/1066.pdf

#### 概要
本文仅包含论文中的protocol 1和protocol 2，以及基于protocol 1证明内积，在本文中被称之为protocol 3。本文不关心电路的部分。简而言之，本文要解决的问题如下：  
Prover有两个Fr中的秘密向量a、b，prover公开com(a)、com(b)、com(c)并且证明c=<a,b>。  
- com(a)=h<sup>r<sub>1</sub></sup>g<sub>1</sub><sup>a</sup>
- com(b)=h<sup>r<sub>2</sub></sup>g<sub>2</sub><sup>b</sup>
- com(c)=h<sup>r<sub>3</sub></sup>u<sup>c</sup>  

#### 术语约定
G是一条椭圆曲线（y<sup>2</sup>=ax<sup>3</sup>+bx+c）上的点构建的群。Fp是x和y的定义域，Fr是该群的阶构造的有限域。  

g、h、u是独立的G生成元，独立生成元的含义是任何人都不知道它们之间的指数关系。例如无人知道方程g<sub>1</sub>=g<sub>2</sub><sup>x</sup>中的x是多少。某些时候为简化起见，会用g表示一个生成元向量，也会用g<sub>1</sub>表示一个生成元向量，用g<sub>2</sub>表示另一个生成元向量。当g表示一个向量时，g<sup>a</sup>表示g[0]<sup>a[0]</sup>g[1]<sup>a[1]</sup>......。  

通过如下方式可以产生满足这种条件的g和h：基于某个hash函数，例如sha256反复hash得到一组随机字节。然后将这些字节映射到G上的某个点，由此得到的点就是g和h。  

可以这样将字节映射到G上的点：首先将随机字节映射到一个Fr，然后以该Fr作为x，求解椭圆曲线方程y<sup>2</sup>=ax<sup>3</sup>+bx+c。在这里需要进行一次有限域上的开方运算，最后得到y。有较小概率ax<sup>3</sup>+bx+c得到的数不是一个平方数，这种情况下反复++x，直到是平方数为止。

### Protocol 2
先介绍Protocol 2而不是Protocol 1是因为Protocol 1基于Protocol 2。  
Prover有两个Fr中的秘密向量a、b，公开P<sub>2</sub>=g<sub>1</sub><sup>a</sup>g<sub>2</sub><sup>b</sup>u<sup><a,b></sup>，证明知道a、b。  
证明细节请参考论文16页。  

### Protocol 1
Prover公开P<sub>1</sub>=g<sub>1</sub><sup>a</sup>g<sub>2</sub><sup>b</sup>，Prover同时公开c，Verifier产生一个随即挑战x，Prover和Verifier同时计算P<sub>2</sub>=g<sub>1</sub><sup>a</sup>g<sub>2</sub><sup>b</sup>u<sup>xc</sup>=P<sub>1</sub>u<sup>xc</sup>。然后Prover运行Protocol 1证明。最终Prover证明了被P<sub>1</sub>锁定的向量a、b的内积是公开值c。 
证明细节请参考论文15页。   
注意首先Protocol 1和Protocol 2是soundness的，但不是zero knownledge的。其次，我们并不希望公开c。

### Protocol 3



