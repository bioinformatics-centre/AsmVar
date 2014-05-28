
g++ TestBamtoolsAPI.cpp utility.cpp -I ~/Bin/MyPipe/VariantDetect/bin/Genotyping/third_party/pezmaster31-bamtools-7f8b301/include -L ~/Bin/MyPipe/VariantDetect/bin/Genotyping/third_party/pezmaster31-bamtools-7f8b301/lib -o TestBamtoolsAPI -lz -lbamtools
echo "Complie Complete"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/Bin/MyPipe/VariantDetect/bin/Genotyping/third_party/pezmaster31-bamtools-7f8b301/lib && ./TestBamtoolsAPI test.bam
