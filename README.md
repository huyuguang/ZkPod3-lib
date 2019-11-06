# pod_all

Core code for https://github.com/sec-bit/zkPoD-lib  

### 获取代码
0，安装依赖库  
ubuntu18.04  
需要gcc 7.x以上  
sudo apt-get update   
sudo apt-get install libcrypto++-dev libcrypto++-doc libcrypto++-utils  
sudo apt-get install libboost-all-dev   

1，拉libsnark代码并拉取依赖，注意不要递归拉取  
git submodule init && git submodule update  
cd depends/libsnark  
git submodule init && git submodule update

2，编译pod_core，pod_setup，pod_publish（linux or osx）  
cd pod_core  
make  

cd ../vrs_cache  
make  

cd ../pod_publish  
make  

编译好的代码在linux/bin目录下。  

3，运行pod_publish发布一个文件  
cd linux/bin  
./pod_publish -m table -f test100000.csv -o table_data -t csv -k 0 1  

4，运行pod_core  
cd linux/bin  
./pod_core -m table -a atomic_swap_pod_vc -p table_data -o table_output --demand_ranges 1-10  


### 编译 (windows + msvc2019)：  
直接用msvc2019打开all.sln。  
