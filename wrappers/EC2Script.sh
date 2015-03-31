#!/usr/bin/env bash

sudo yum -y install mysql-server gcc tmux
sudo chkconfig mysqld on
sudo service mysqld start
curl -O ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.2.29+-x64-linux.tar.gz
tar xzf ncbi-*.tar.gz
sudo cp ncbi-blast-2.2.29+/bin/* /usr/local/bin

curl -O http://www.micans.org/mcl/src/mcl-latest.tar.gz
tar xzf mcl*.tar.gz
cd mcl-12-068
./configure
make
sudo make install
cd ..

curl -O http://orthomcl.org/common/downloads/software/v2.0/orthomclSoftware-v2.0.9.tar.gz
tar xzf orthomcl*.tar.gz

mkdir my_othromcl_dir
cd my_othromcl_dir

# this is my orthoMCL config file
curl -O http://chrisprobert.com.s3.amazonaws.com/chlamy/orthomcl.config.template
mv orthomcl.config.template orthomcl.config

# this is my "all proteins" input fasta
curl -O http://chrisprobert.com.s3.amazonaws.com/chlamy/allProts.fa

mkdir compliantFasta
mkdir inputFasta
mv allProts.fa inputFasta

../orthomclSoftware-v2.0.9/bin/orthomclInstallSchema my_othromcl_dir/orthomcl.config
../orthomclSoftware-v2.0.9/bin/orthomclAdjustFasta ecu inputFasta/allProts.fa 1
mv *.fasta compliantFasta
../orthomclSoftware-v2.0.9/bin/orthomclFilterFasta compliantFasta 10 20

makeblastdb -in goodProteins.fasta -dbtype prot -out goodProteinsDB
blastp -query goodProteins.fasta -db goodProteinsDB -outfmt 6 -num_threads 8 -out blastOut

../orthomclSoftware-v2.0.9/bin/orthomclBlastParser blastOut compliantFasta/ >> similarSequences.txt
../orthomclSoftware-v2.0.9/bin/orthomclLoadBlast orthomcl.config similarSequences.txt
../orthomclSoftware-v2.0.9/bin/orthomclPairs orthomcl.config orthoLogFile cleanup=yes
../orthomclSoftware-v2.0.9/bin/orthomclDumpPairsFiles orthomcl.config

mcl mclInput --abc -I 1.5 -o mclOutput

../orthomclSoftware-v2.0.9/bin/orthomclMclToGroups groupID 1 < mclOutput > groups.txt

cd ..