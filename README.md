# pod_all

Core code for https://github.com/sec-bit/zkPoD-lib  

##### 安装依赖库  
ubuntu18.04  
需要gcc 7.x以上  
sudo apt-get update   
sudo apt-get install libcrypto++-dev libcrypto++-doc libcrypto++-utils  
sudo apt-get install libboost-all-dev   

如果希望用TBB处理并行，那么还需要安装TBB库（可选，使用TBB性能会略微高5%~10%）：  
sudo apt install libtbb-dev  
或者可以参考如下URL：  
https://askubuntu.com/questions/1170054/install-newest-tbb-thread-building-blocks-on-ubuntu-18-04
https://stackoverflow.com/questions/3181468/how-do-you-install-intel-tbb-on-os-x

##### 拉libsnark代码并拉取依赖，注意不要递归拉取  
git submodule init && git submodule update  
cd depends/libsnark  
git submodule init && git submodule update

##### 编译pod_core，pod_setup，pod_publish（linux or osx）  
cd pod_core  
make -j8 或者 make -j8 USE_TBB=1  

cd ../vrs_cache  
make -j8 或者 make -j8 USE_TBB=1  

cd ../pod_publish  
make -j8 或者 make -j8 USE_TBB=1  

编译好的代码在linux/bin目录下。  

##### 运行pod_publish发布一个文件  
cd linux/bin  
./pod_publish -m plain -f test.txt -o plain_data -c 1023  
./pod_publish -m table -f test100000.csv -o table_data -t csv -k 0 1  
有两种模式，一种是plain，用于无结构的数据文件，另一种是数据表，目前仅支持csv格式。  
运行完后会在plain_data或者table_data目录下产生发布文件。cat bulletin可以查看相关信息。  
例如：  
{  
    "mode": "plain",  
    "size": "140614293",  
    "s": "1024",  
    "n": "4434",  
    "sigma_mkl_root":   "28accd27712939fb46e6313d18d67af36a523f44659ced9c27e2d3fd0e9cd3a5"  
}  

##### 运行pod_core  
cd linux/bin  
./pod_core -m plain -a atomic_swap_pod_vc -p plain_data -o plain_output --demand_ranges 0-10  
./pod_core -m table -a atomic_swap_pod_vc -p table_data -o table_output --demand_ranges 1-10  

demand_ranges支持多个区间。1-10表示从第2行开始取10行。  
如果希望获取整个文件，那么是0-n，n可以从bulletin文件中获得。

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
