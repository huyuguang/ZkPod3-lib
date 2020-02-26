# zkpod-clink

##### 安装依赖库  
- Linux(Ubuntu18 or above)
  - sudo apt-get update   
  - sudo apt-get install libcrypto++-dev
  - sudo apt-get install libboost-all-dev   
  - sudo apt-get install libgmp-dev
  - sudo apt install libtbb-dev  

- Mac
  - brew install libcryptopp
  - brew install boost
  - brew install gmp
  - brew install tbb

- Windows
  1. 安装vcpkg，最好安装在c:/目录下，否则将需要修改./win/props/目录下vcpkg.64.debug.props以及vcpkg.64.release.props里相应的路径。
  2. 用vcpkg安装crypto++、boost、mpir、tbb库：  
      - vcpkg install cryptopp:x64-windows
      - vcpkg install boost:x64-windows
      - vcpkg install mpir:x64-windows
      - vcpkg install tbb:x64-windows 

##### 拉libsnark代码并拉取依赖，注意不要递归拉取  
1. git submodule init && git submodule update
2. cd depends/libsnark
3. git submodule init && git submodule update

##### 编译
- Linux or Mac
  - 编译mcl: cd depends/libsnark/depends/mcl; make clean && make lib/libmcl.a MCL_USE_OPENSSL=0
  - 需要c++17支持
  - cd pod_dummy && make -j8  
  - 编译好的可执行程序在linux/bin目录下
- Windows
  - 不需要单独编译mcl
  - 需要msvc2017或msvc2019
  - 直接用msvc2019打开all.sln，编译即可
  - 编译好的可执行程序在win/bin目录下

##### 运行pod_dummy  
./pod_dummy -h  
具体每个功能的详细说明请参看相关代码所在目录的readme.md。
