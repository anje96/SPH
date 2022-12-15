# Set compiler
CC := g++

# global include
INCLDIR := include
HDF5DIR :=
HDF5LIB := 

# resource on flags: https://stackoverflow.com/questions/2855121/what-is-the-purpose-of-using-pedantic-in-the-gcc-g-compiler
# Set compiler flags
CFLAGS := -std=c++11 -pedantic -Wall -O3 -I$(INCLDIR) -I$(HDF5DIR)
# Set linker flags
LDFLAGS := -std=c++11 -lncurses -lm -L$(HDF5LIB) -lhdf5

# Define extensions of files
SRCEXT := cpp
OBJEXT := o
HEADEXT := h
DEPEXT := d

# Define directories to look for certain types of files
SRCDIR := src
OBJDIR := obj
DEPDIR := dep

# Autogenerate lists for all files in the project 
SOURCES := $(wildcard $(SRCDIR)/*.$(SRCEXT))
DEPENDS := $(SOURCES:$(SRCDIR)/%.$(SRCEXT)=$(DEPDIR)/%.$(DEPEXT))
OBJECTS := $(SOURCES:$(SRCDIR)/%.$(SRCEXT)=$(OBJDIR)/%.$(OBJEXT))

# get kernel name to be able to run sed correctly on Darwin (MacOS) or Linux kernels
KERNEL := $(shell uname -s)
ifeq ($(KERNEL), Darwin) 
	SED := sed -i "~"
else
	SED := sed -i
endif

# Define executable directory and executable file
EXEDIR := bin
EXE := exe

# MAIN TARGET
all: $(OBJECTS) $(EXEDIR)/$(EXE)

# build all with debug flags
debug: CFLAGS := -g -std=c99 -pedantic -Wall -O0 -I$(INCLDIR)
# show linker invocation when building debug target
debug: LDFLAGS += -v
debug: all

# define phony target for cleaning project
# use this e.g. to rebuild for debugging
clean:
	@rm $(EXEDIR)/$(EXE) $(OBJECTS) $(DEPENDS)

# use this to delete all autogenerated files and directories
superclean: 
	@rm -r $(DEPDIR) $(OBJDIR) $(EXEDIR)

# autogenerate directories if they don't exist
$(EXEDIR):
	@mkdir -p $@

$(DEPDIR):
	@mkdir -p $@

$(OBJDIR):
	@mkdir -p $@

# generating dependency files for all sources
# sed changes '%.o: ...' to '%.o %.d: ...' in dependency file
$(DEPDIR)/%.$(DEPEXT): $(SRCDIR)/%.$(SRCEXT) | $(DEPDIR)
	@echo "Generating dependency file '$@' ..."
	@$(CC) -MM $(CFLAGS) $< -MF $@
	@$(SED) 's,$(*F)\.$(OBJEXT),$*\.$(OBJEXT) $@,' $@
	@rm -f $@~
	@echo "... done."

# include targets from generated dependency files
-include $(DEPENDS)

# build main target
# check if target directory 'bin' already exist via prerequisite
$(EXEDIR)/$(EXE): $(OBJECTS) | $(EXEDIR)
	@echo "Linking binary '$@'..."
	@$(CC) $(CFLAGS) $(OBJECTS) -o $(EXEDIR)/$(EXE) $(LDFLAGS)
	@echo "... done."

# pattern rule to build object from source file
$(OBJDIR)/%.$(OBJEXT): $(DEPDIR)/%.$(DEPEXT) | $(OBJDIR)
	@echo "Compiling '$@' ..."
	@$(CC) -c $(CFLAGS) $(@:$(OBJDIR)/%.$(OBJEXT)=$(SRCDIR)/%.$(SRCEXT)) -o $@
	@echo "... done."

.PHONY: all debug clean cleaner 
