# Bulletproofs

https://eprint.iacr.org/2017/1066.pdf

#### ��Ҫ
���Ľ����������е�protocol 1��protocol 2���Լ�����protocol 1֤���ڻ����ڱ����б���֮Ϊprotocol 3�����Ĳ����ĵ�·�Ĳ��֡������֮������Ҫ������������£�  
Prover������Fr�е���������a��b��prover����com(a)��com(b)��com(c)����֤��c=<a,b>��  
- com(a)=h<sup>r<sub>1</sub></sup>g<sub>1</sub><sup>a</sup>
- com(b)=h<sup>r<sub>2</sub></sup>g<sub>2</sub><sup>b</sup>
- com(c)=h<sup>r<sub>3</sub></sup>u<sup>c</sup>  

#### ����Լ��
G��һ����Բ���ߣ�y<sup>2</sup>=ax<sup>3</sup>+bx+c���ϵĵ㹹����Ⱥ��Fp��x��y�Ķ�����Fr�Ǹ�Ⱥ�Ľ׹����������  

g��h��u�Ƕ�����G����Ԫ����������Ԫ�ĺ������κ��˶���֪������֮���ָ����ϵ����������֪������g<sub>1</sub>=g<sub>2</sub><sup>x</sup>�е�x�Ƕ��١�ĳЩʱ��Ϊ�����������g��ʾһ������Ԫ������Ҳ����g<sub>1</sub>��ʾһ������Ԫ��������g<sub>2</sub>��ʾ��һ������Ԫ��������g��ʾһ������ʱ��g<sup>a</sup>��ʾg[0]<sup>a[0]</sup>g[1]<sup>a[1]</sup>......��  

ͨ�����·�ʽ���Բ�����������������g��h������ĳ��hash����������sha256����hash�õ�һ������ֽڡ�Ȼ����Щ�ֽ�ӳ�䵽G�ϵ�ĳ���㣬�ɴ˵õ��ĵ����g��h��  

�����������ֽ�ӳ�䵽G�ϵĵ㣺���Ƚ�����ֽ�ӳ�䵽һ��Fr��Ȼ���Ը�Fr��Ϊx�������Բ���߷���y<sup>2</sup>=ax<sup>3</sup>+bx+c����������Ҫ����һ���������ϵĿ������㣬���õ�y���н�С����ax<sup>3</sup>+bx+c�õ���������һ��ƽ��������������·���++x��ֱ����ƽ����Ϊֹ��

### Protocol 2
�Ƚ���Protocol 2������Protocol 1����ΪProtocol 1����Protocol 2��  
Prover������Fr�е���������a��b������P<sub>2</sub>=g<sub>1</sub><sup>a</sup>g<sub>2</sub><sup>b</sup>u<sup><a,b></sup>��֤��֪��a��b��  
֤��ϸ����ο�����16ҳ��  

### Protocol 1
Prover����P<sub>1</sub>=g<sub>1</sub><sup>a</sup>g<sub>2</sub><sup>b</sup>��Proverͬʱ����c��Verifier����һ���漴��սx��Prover��Verifierͬʱ����P<sub>2</sub>=g<sub>1</sub><sup>a</sup>g<sub>2</sub><sup>b</sup>u<sup>xc</sup>=P<sub>1</sub>u<sup>xc</sup>��Ȼ��Prover����Protocol 1֤��������Prover֤���˱�P<sub>1</sub>����������a��b���ڻ��ǹ���ֵc�� 
֤��ϸ����ο�����15ҳ��   
ע������Protocol 1��Protocol 2��soundness�ģ�������zero knownledge�ġ���Σ����ǲ���ϣ������c��

### Protocol 3



