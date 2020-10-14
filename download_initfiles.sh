# Need to figure out a way to install CheckM via conda and place the 
pip3 install checkm-genome
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvzf https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz -C databases/checkm_db

checkm data setRoot databases/checkm_db