/home/siyang/Bin/software_pip/YH_script/AsmVar/lastdb -c -m1111110 -v b36_1_1000000_3000000.lastdb b36_1_1000000_3000000.fa 
/home/siyang/Bin/MyPipe/Assembly/develop/1.1.0/third_party/variantstudio/bin/lastal -e25 -v -q3 -j4 b36_1_1000000_3000000.lastdb query.fa | /home/siyang/Bin/MyPipe/Assembly/develop/1.1.0/third_party/variantstudio/bin/last-split -s35 -v > align.maf 2> align.log
rm b36_1_1000000_3000000.lastdb.*
