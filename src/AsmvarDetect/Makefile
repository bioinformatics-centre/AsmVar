#INCLUDE_DIR = $(CURDIR)/../include
#UTILS_DIR   = $(INCLUDE_DIR)/utils

CC  = g++ -lz
obj = main.o utility.o Variant.o VarUnit.o gzstream.o AGEaligner.o

ASV_VariantDetector : $(obj)

	$(CC) -o out3 $(obj)

main.o : gzstream.h utility.h Variant.h
VarUnit.o  : VarUnit.h
AGEaligner.o : AGEaligner.h
utility.o  : utility.h
gzstream.o : gzstream.h
Variant.o  : Variant.h

.PHONY : clean

clean:
	rm $(obj)

