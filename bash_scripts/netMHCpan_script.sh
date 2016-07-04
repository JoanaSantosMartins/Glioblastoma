#!/bin/bash

tput setaf 1
#

echo "Hello!"
echo """Please make sure you have the file base.txt in this directory, and the input file (example idh1.fasta)
If you lost base.txt file you can recover it from this script ( at the end of netMHCpan_script.sh)
If your input file has a name differnt from idh1.fasta you have to replace it in base.txt"""
echo "Starting..."
echo "Getting human HLAs for netMHCpan"

# get human HLAs for netmhcIIpan
/home/manager/netMHC/netMHCpan-2.8/netMHCpan -listMHC | grep "HLA" | grep -v "#" > human_HLA.txt

echo "Done, file human_hla_netmhcpan.txt created!"
echo "Creating run file!"

# create run file
cat human_HLA.txt | while read line
do
   sed "s/old/'$line'/g" base.txt >> run_file.sh
done

sed -i "s/'//g" run_file.sh
echo 'mkdir out' | cat - run_file.sh > temp && mv temp run_file.sh
echo '#!/bin/bash' | cat - run_file.sh > temp && mv temp run_file.sh
chmod 777 run_file.sh

echo "Run file created as run_file.sh"
echo "To run the job pls input ./run_file.sh"
echo "My job here is Done"


<<COMMENT1
base.txt -> copy the following:
../netMHCpan -a old -f 8peptide_IDH1.fasta -l 8 -s -th 50.0 -lt 500.0 -rth 0.50 -rlt 2.0 -ic50 > ./out/old_8.txt
../netMHCpan -a old -f 9peptide_IDH1.fasta -l 9 -s -th 50.0 -lt 500.0 -rth 0.50 -rlt 2.0 -ic50 > ./out/old_9.txt
../netMHCpan -a old -f 10peptide_IDH1.fasta -l 10 -s -th 50.0 -lt 500.0 -rth 0.50 -rlt 2.0 -ic50 > ./out/old_10.txt
../netMHCpan -a old -f 11peptide_IDH1.fasta -l 11 -s -th 50.0 -lt 500.0 -rth 0.50 -rlt 2.0 -ic50 > ./out/old_11.txt
../netMHCpan -a old -f 12peptide_IDH1.fasta -l 12 -s -th 50.0 -lt 500.0 -rth 0.50 -rlt 2.0 -ic50 > ./out/old_12.txt
../netMHCpan -a old -f 13peptide_IDH1.fasta -l 13 -s -th 50.0 -lt 500.0 -rth 0.50 -rlt 2.0 -ic50 > ./out/old_13.txt
../netMHCpan -a old -f 14peptide_IDH1.fasta -l 14 -s -th 50.0 -lt 500.0 -rth 0.50 -rlt 2.0 -ic50 > ./out/old_14.txt
COMMENT1

tput sgr0


