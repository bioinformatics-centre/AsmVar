#INCLUDE_DIR = $(CURDIR)/../include
#UTILS_DIR   = $(INCLUDE_DIR)/utils

CC  = g++ -lz
obj = main.o utility.o Variant.o VarUnit.o gzstream.o AGEaligner.o vcf.o Region.o

ASV_VariantDetector : $(obj)

	$(CC) -o ASV_VariantDetector $(obj)

main.o : gzstream.h utility.h Variant.h
VarUnit.o  : utility.h VarUnit.h AGEaligner.h Sequence.h Scorer.h Region.h
AGEaligner.o : AGEaligner.h Sequence.h Scorer.h
Variant.o  : Variant.h VarUnit.h vcf.h Region.h
Region.o   : Region.h
vcf.o      : vcf.h utility.h
utility.o  : utility.h
gzstream.o : gzstream.h

.PHONY : clean

clean:
	rm $(obj)

