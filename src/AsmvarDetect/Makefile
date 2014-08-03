#INCLUDE_DIR = $(CURDIR)/../include
#UTILS_DIR   = $(INCLUDE_DIR)/utils

CC  = g++ -lz
obj = main.o utility.o Variant.o gzstream.o

out : $(obj)

	$(CC) -o out $(obj)

main.o : gzstream.h utility.h Variant.h
utility.o  : utility.h
gzstream.o : gzstream.h
Variant.o  : Variant.h

.PHONY : clean

clean:
	rm $(obj)

