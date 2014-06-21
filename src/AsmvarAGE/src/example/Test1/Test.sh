/home/siyang/Bin/MyPipe/Assembly/develop/1.1.0/third_party/variantstudio/bin/AGE/age_align  -indel -both -match=1 -mismatch=-1 -go=-10 -ge=-1 -coor1=1-1000 -coor2=1-534 target.1.fa query.1.fa > t1.age

../../age_align --indel --both --match=1 --mismatch=-1 --go=-10 --ge=-1 --target target.1.fa --query query.1.fa var.txt
../../age_align -i -b -m 1 -M -1 -g -10 -e -1 -t target.1.fa -q query.1.fa var.txt > t2.age
