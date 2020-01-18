# pod_core

-m plain -a complaint_pod -p plain_data -o plain_output --demand_ranges 0-10  
-m table -a complaint_pod -p table_data -o table_output --demand_ranges 1-20  

-m plain -a atomic_swap_pod -p plain_data -o plain_output --demand_ranges 1-10  
-m table -a atomic_swap_pod -p table_data -o table_output --demand_ranges 1-10  

-m plain -a atomic_swap_pod_vc -p plain_data -o plain_output --demand_ranges 1-10  
-m table -a atomic_swap_pod_vc -p table_data -o table_output --demand_ranges 1-10  

-m plain -a ot_complaint_pod -p plain_data -o plain_output --demand_ranges 1-2 --phantom_ranges 0-3  
-m table -a ot_complaint_pod -p table_data -o table_output --demand_ranges 1-2 --phantom_ranges 0-3  

-m table -a vrf_query -p table_data -o table_output -k first_name -v Kathy  
-m table -a vrf_query -p table_data -o table_output -k "Emp ID" -v 614227  
-m table -a ot_vrf_query -p table_data -o table_output -k "Emp ID" -v 313736 964888 abc -n 350922 aaa eee bbb  

--dump_ecc_pub  



https://askubuntu.com/questions/1170054/install-newest-tbb-thread-building-blocks-on-ubuntu-18-04  
To upgrade to the latest version, please do the following:

Add the Ubuntu repository that contains the latest version 2019~U8-1, run the following command in terminal:

echo "deb http://cz.archive.ubuntu.com/ubuntu eoan main universe" | sudo tee -a  /etc/apt/sources.list

Update the repositories, run the following command in terminal:

sudo apt update

Upgrade to the latest version, run the following command in terminal:

sudo apt install libtbb-dev

After this you should have the latest libtbb-dev installed.

