# pod_all

Core code for https://github.com/sec-bit/zkPoD-lib  

##### 安装依赖库  
ubuntu18.04  
需要gcc 7.x以上  
sudo apt-get update   
sudo apt-get install libcrypto++-dev libcrypto++-doc libcrypto++-utils  
sudo apt-get install libboost-all-dev   

安装tbb
sudo apt install libtbb-dev  
或者可以参考如下URL：  
https://askubuntu.com/questions/1170054/install-newest-tbb-thread-building-blocks-on-ubuntu-18-04
https://stackoverflow.com/questions/3181468/how-do-you-install-intel-tbb-on-os-x

##### 拉libsnark代码并拉取依赖，注意不要递归拉取  
git submodule init && git submodule update  
cd depends/libsnark  
git submodule init && git submodule update

##### 编译pod_dummy，vrs_cache（linux or osx）  
cd pod_dummy  
make -j8  

cd vrs_cache  
make -j8  


编译好的代码在linux/bin目录下。  

##### 运行pod_dummy  
cd linux/bin  
./pod_dummy -h  

##### vrs_cache
vrs_cache的作用是预先生成一些gro09的commitment。比如说如果希望获取10行，也即demand_ranges x-10，且该文件每行是1024（bulletin里的s），那么 -c 10240。  
这会在vrs_cache目录下产生一个cache文件。  
pod_core运行时会自动寻找这个文件，如果存在，那么会把文件后缀改成using，在交易完成后会改成used以表明不能再用。如果交易失败则会删除后缀以重复使用。  
通常并不能预知用户要获取的demand_ranges，所以vrs_cache的策略是建立一组cache文件，pod_core会自动匹配最接近地cache文件，并且基于该cache文件进行修补，同样能起到很好的加速效果。  
注意交易成功之后cache文件后缀会被改成used，从而不可再用。因此测试前需要把used后缀去掉。

##### Windows平台下的编译 (msvc2019)：  
1，安装vcpkg，最好安装在c:/目录下，否则将需要修改./win/props/下的两个props文件中的路径。  
2，用vcpkg安装crypto++和boost  
cd c:/vcpkg  
vcpkg install boost boost:x64-windows boost:x86-windows-static boost:x64-windows-static  
vcpkg install cryptopp cryptopp:x64-windows cryptopp:x86-windows-static cryptopp:x64-windows-static  
vcpkg install tbb tbb:x64-windows tbb:x86-windows-static tbb:x64-windows-static  

3，直接用msvc2019打开all.sln，编译。  
