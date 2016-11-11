#!/bin/bash

tput setaf 1
echo "Hello!"
echo """Please make sure you have the file base.txt in this directory, and the input file (example idh1.fasta)
If you lost base.txt file you can recover it from this script ( at the end of netMHCIIpan_script.sh)
If your input file has a name differnt from idh1.fasta you have to replace it in base.txt"""
echo "Starting..."
echo "Getting human HLAs for netMHCpan"

# get human HLAs for netmhcpan
../netMHCIIpan -list | grep "DRB" | grep -v "#" > human_HLA.txt 


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
../netMHCIIpan -a old -f 9peptide_idh1.fasta -length 9 -affS 50.0 -affW 500.0 -rankS 0.50 -rankW 2.0 > ./out/old_9.txt
../netMHCIIpan -a old -f 10peptide_idh1.fasta -length 10 -affS 50.0 -affW 500.0 -rankS 0.50 -rankW 2.0 > ./out/old_10.txt
../netMHCIIpan -a old -f 11peptide_idh1.fasta -length 11 -affS 50.0 -affW 500.0 -rankS 0.50 -rankW 2.0 > ./out/old_11.txt
../netMHCIIpan -a old -f 12peptide_idh1.fasta -length 12 -affS 50.0 -affW 500.0 -rankS 0.50 -rankW 2.0 > ./out/old_12.txt
../netMHCIIpan -a old -f 13peptide_idh1.fasta -length 13 -affS 50.0 -affW 500.0 -rankS 0.50 -rankW 2.0 > ./out/old_13.txt
../netMHCIIpan -a old -f 14peptide_idh1.fasta -length 14 -affS 50.0 -affW 500.0 -rankS 0.50 -rankW 2.0 > ./out/old_14.txt
../netMHCIIpan -a old -f 15peptide_idh1.fasta -length 15 -affS 50.0 -affW 500.0 -rankS 0.50 -rankW 2.0 > ./out/old_15.txt
../netMHCIIpan -a old -f 16peptide_idh1.fasta -length 16 -affS 50.0 -affW 500.0 -rankS 0.50 -rankW 2.0 > ./out/old_16.txt
../netMHCIIpan -a old -f 17peptide_idh1.fasta -length 17 -affS 50.0 -affW 500.0 -rankS 0.50 -rankW 2.0 > ./out/old_17.txt
../netMHCIIpan -a old -f 18peptide_idh1.fasta -length 18 -affS 50.0 -affW 500.0 -rankS 0.50 -rankW 2.0 > ./out/old_18.txt
../netMHCIIpan -a old -f 19peptide_idh1.fasta -length 19 -affS 50.0 -affW 500.0 -rankS 0.50 -rankW 2.0 > ./out/old_19.txt
../netMHCIIpan -a old -f 20peptide_idh1.fasta -length 20 -affS 50.0 -affW 500.0 -rankS 0.50 -rankW 2.0 > ./out/old_20.txt

COMMENT1

tput sgr0


