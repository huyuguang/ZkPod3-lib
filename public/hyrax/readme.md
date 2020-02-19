# Hyrax

https://eprint.iacr.org/2017/1132.pdf

### ��Ҫ
���Ľ��������ĸ�¼�е�A1��A2�Լ�A3�㷨�����Ĳ�����ԭ�����к�sumcheck��GKR�Լ���·��صĲ��֡�  
[��������ͷ���Լ��](../ecc/readme.md)

### A1
- x��y��z��F<sub>r</sub>��3�����ܵı���
- Prover����com(x)��com(y)��com(z)��֤��z=<x,y>
  - com(x)=h<sup>r<sub>x</sub></sup>g<sub>x</sub><sup>x</sup>
  - com(y)=h<sup>r<sub>y</sub></sup>g<sub>y</sub><sup>y</sup>
  - com(z)=h<sup>r<sub>z</sub></sup>g<sub>z</sub><sup>z</sup>
  - g<sub>x</sub>��g<sub>y</sub>������ͬ

### A2
- a�ǹ�����F<sub>r</sub>��������������Ϊn
- x�����ܵ�F<sub>r</sub>��������������Ϊn
- y�����ܵ�F<sub>r</sub>��������y=<x,a>
- Prover����com(x)��com(y)��֤��y=<x,a>
  - com(x)=h<sup>r<sub>x</sub></sup>g<sub>x</sub><sup>x</sup>
  - com(y)=h<sup>r<sub>y</sub></sup>g<sub>y</sub><sup>y</sup>
  - g<sub>x</sub>����������g<sub>y</sub>�������������ص�

### A3
A3��A2�ĵݹ��۵��汾���и����֤����������֤������ͬʱ�ж��������֤�����ȡ�

### ����
- A1
  - Prove cost��6EccExp
  - Verify cost��6EccExp
  - Proof size��3G1+5F<sub>r</sub>
- A2
  - Prove cost��MulExp(n)
  - Verify cost��MulExp(n)
  - Proof size��2G1+(n+2)F<sub>r</sub>
- A3
  - Prove cost�� 2MulExp(n)~4MulExp(n)
  - Verify cost�� 2nEccExp
  - Proof size��(2+2log(n))G1+2F<sub>r</sub>
