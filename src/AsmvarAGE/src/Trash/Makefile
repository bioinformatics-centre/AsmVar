VERSION = v0.4
DEFAULT_FLAGS = -DAGE_VERSION=\"$(VERSION)\" -DAGE_TIME
CXX	= g++ -fopenmp -DOMP $(DEFAULT_FLAGS)

MESS = "Compiling with parallel support."
ifeq ($(OMP),no)
        CXX  = g++ $(DEFAULT_FLAGS)
	MESS = "Compiling with NO parallel support."
endif

OBJDIR = obj
OBJS   = $(OBJDIR)/AGEaligner.o \
	 $(OBJDIR)/age_align.o

DISTRIBUTION    = $(PWD)/AGE_$(VERSION).zip
TESTFILE        = TestDecription_$(VERSION).txt
TMPDIR  	= /tmp
AGEDIR   	= age_$(VERSION)
MAINDIR		= $(TMPDIR)/$(AGEDIR)
SRCDIR   	= $(MAINDIR)/src
EXAMPLEDIR 	= $(MAINDIR)/examples


all: mess age_align

mess:
	@echo ""
	@echo $(MESS)
	@echo ""

age_align: $(OBJS)
	$(CXX) -o $@ $(OBJS)

$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(OBJDIR)
	$(CXX) -c $< -o $@

clean:
	rm -f $(OBJS)

distribution: clean all
	@echo Creating directory ...
	@rm -rf $(MAINDIR)
	@mkdir  $(MAINDIR)
	@mkdir  $(SRCDIR)
	@mkdir  $(EXAMPLEDIR)
	@echo Making tests ...
	@./test.sh > $(TESTFILE)
	@echo Copying files ...
	@cp *.cpp *.hh    $(SRCDIR)
	@cp Makefile      $(SRCDIR)
	@cp README        $(MAINDIR)
	@cp license.rtf   $(MAINDIR)
	@cp test*.fa $(TESTFILE) $(EXAMPLEDIR)
	@echo Zipping ...
	@ln -s $(MAINDIR)
	@zip -qr $(DISTRIBUTION) $(AGEDIR)
	@rm $(AGEDIR)
