SOURCES:=$(shell ls *.cpp)

OBJS:=$(addprefix $(OBJ_DIR)/,$(SOURCES:.cpp=.o))

all: $(OBJS)

$(OBJS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(CXXFLAGS) $(INCLUDES)
