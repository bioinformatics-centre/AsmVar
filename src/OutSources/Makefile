BAM_LIB:=$(LIB_DIR)/libbamtools.a
SAM_LIB:=$(LIB_DIR)/libbam.a
SSW:=$(OBJ_DIR)/ssw.o $(OBJ_DIR)/ssw_cpp.o
FASTA:=$(OBJ_DIR)/Fasta.o
UTIL:=$(OBJ_DIR)/md5.o

libs: $(SSW) $(FASTA) $(BAM_LIB) $(SAM_LIB) $(UTIL)

$(SSW): stripedSW/*
	@$(MAKE) all --no-print-directory --directory=stripedSW;

$(FASTA): fasta/*
	@$(MAKE) all --no-print-directory --directory=fasta;

$(UTIL): util/*
	@$(MAKE) all --no-print-directory --directory=util;

$(BAM_LIB):
	@echo "  * generating libbamtools.a"
	@cd bamtools/ && mkdir -p build && cd build && cmake .. && $(MAKE)
	@echo "  * copying libbamtools.a to the library directory"
	@cp bamtools/lib/libbamtools.a $(LIB_DIR)
	@echo "  * copying include files to the library directory"
	@cp -r bamtools/include/* $(INCLUDE_DIR)
	@echo " "

$(SAM_LIB):
	@echo "  * generating libbam.a"
	@$(MAKE) lib --no-print-directory -C samtools
	@echo "  * copying libbam.a to the library directory"
	@cp samtools/libbam.a $(LIB_DIR)
	@echo "  * copying include files to the library directory"
	@cp samtools/bgzf.h samtools/bam.h samtools/khash.h samtools/kseq.h $(INCLUDE_DIR)

clean:
	@rm -rf bamtools/bin bamtools/build bamtools/lib/ bamtools/include
	@$(MAKE) clean --no-print-directory -C samtools

.PHONY: all libs clean

