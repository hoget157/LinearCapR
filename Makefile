PROG := LinCapR
OBJDIR := ./temp
SRCS := $(wildcard *.cpp)
OBJS := $(addprefix $(OBJDIR)/, $(SRCS:%.cpp=%.o))
DEPS := $(addprefix $(OBJDIR)/, $(SRCS:%.cpp=%.d))

ifeq ($(OS),Windows_NT)
OBJDIR := .\temp
OBJS = $(addprefix $(OBJDIR)\, $(SRCS:%.cpp=%.o))
DEPS = $(addprefix $(OBJDIR)\, $(SRCS:%.cpp=%.d))
endif

CXX ?= $(if $(filter $(OS),Windows_NT),,$(shell if command -v clang++ >/dev/null 2>&1; then echo clang++; elif command -v g++ >/dev/null 2>&1; then echo g++; fi))
ifeq ($(strip $(CXX)),)
CXX := g++
endif

CXXFLAGS := -O3 -std=c++17 -Wall #-pg -g
# INCLUDEPATH := -I/usr/local/include
# LIBPATH := -L/usr/local/lib
# LIBS := -framework Cocoa -framework OpenGL -lz -ljpeg -lpng

all: $(DEPENDS) $(PROG)

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBPATH) $(LIBS)

ifeq ($(OS),Windows_NT)
$(OBJDIR)\\%.o: %.cpp
	if not exist temp mkdir temp
	$(CXX) $(CXXFLAGS) $(INCLUDEPATH) -MMD -MP -MF $(<:%.cpp=temp\\%.d) -c $< -o $(<:%.cpp=temp\\%.o)
else
$(OBJDIR)/%.o: %.cpp
	mkdir -p temp
	$(CXX) $(CXXFLAGS) $(INCLUDEPATH) -MMD -MP -MF $(<:%.cpp=temp/%.d) -c $< -o $(<:%.cpp=temp/%.o)
endif

.PHONY: clean
clean:
ifeq ($(OS),Windows_NT)
	del $(PROG).exe $(OBJS) $(DEPS)
else
	rm $(PROG) $(OBJS) $(DEPS)
endif

-include $(DEPS)
