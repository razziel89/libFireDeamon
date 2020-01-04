# makefile for libFireDeamon

#import build parameters
include make.vars
#-----------------------------------------------------
#         SOURCE FOR MAIN C++ LIBRARY
#-----------------------------------------------------
MAINSRC	:= \
	$(SRCDIR)/parallel_generic.cpp $(SRCDIR)/irregular_grid_interpolation.cpp \
	$(SRCDIR)/arbitrary_grid_local_minima.cpp $(SRCDIR)/set_procname.cpp
ifeq ($(MINIMUMDEPS), false) 
MAINSRC += \
	$(SRCDIR)/skin_surface_deamon.cpp $(SRCDIR)/electrostatic_potential_charges.cpp $(SRCDIR)/electron_density.cpp \
	$(SRCDIR)/isosurface.cpp $(SRCDIR)/halfnum/angular_integral.cpp $(SRCDIR)/halfnum/radial_integral.cpp \
	$(SRCDIR)/electrostatic_potential_orbitals.cpp $(SRCDIR)/constants.cpp $(SRCDIR)/orbital_overlap.cpp
endif
MAINOBJ  := $(MAINSRC:$(SRCDIR)%.cpp=$(OBJDIR)%.o)
MAINHDR  := $(MAINSRC:$(SRCDIR)%.cpp=$(PREFIX)/include/FireDeamon%.h)
#-----------------------------------------------------
#         SOURCE FOR TEST FILES
#-----------------------------------------------------
#First Test
TEST1SRC  := $(TESTDIR)/test1.cpp 
TEST1OBJ  := $(TEST1SRC:cpp=o)
#Second Test
TEST2SRC  := $(TESTDIR)/test2.cpp
TEST2OBJ  := $(TEST2SRC:cpp=o)
#Third Test
TEST3SRC  := $(TESTDIR)/test3.cpp
TEST3OBJ  := $(TEST3SRC:cpp=o)
#Fourth Test
TEST4SRC  := $(TESTDIR)/test4.cpp
TEST4OBJ  := $(TEST4SRC:cpp=o)
#------------------------------------------------------
# 	DEFAULT VARIABLES WITHOUT BINDINGS
#------------------------------------------------------
RUNSWIG   :=
RUNBIND   :=
INSTALL   := main_install
UNINSTALL := main_uninstall
LDFLAGS   := $(LDFLAGSMIN)
CPPFLAGS  := $(CPPFLAGSMIN)
TEST      :=
ifeq ($(MINIMUMDEPS), false) 
TEST      += test
LDFLAGS   += $(LDFLAGSMAX)
CPPFLAGS  += $(CPPFLAGSMAX)
endif
#-----------------------------------------------------
#    VARIABLE DEFINITIONS IF PYTHON IS REQUESTED
#-----------------------------------------------------
ifeq ($(USEPYTHON), true) 
PYTHONBINDCXX     := $(PYTHONDIR)/FireDeamon_wrap.cxx
PYTHONBINDCXXCP   := $(PYTHONDIR)/FireDeamon_wrap.cxx.min
PYTHONBINDPY      := $(PYTHONDIR)/cpp.py
PYTHONBINDPYCP    := $(PYTHONDIR)/cpp.py.min
PYTHONBINDOBJ     := $(PYTHONDIR)/FireDeamon_wrap.o
PYTHONCAT         := $(PYTHONDIR)/FireDeamonMin.py
PYTHONCATDEST     := $(PYTHONDIR)/FireDeamon.py
PYTHONBINDCAT     := $(PYTHONDIR)/min.i
PYTHONBINDCATDEST := $(PYTHONDIR)/cpp.i
INSTALL           += python_install
UNINSTALL         += python_uninstall
RUNBIND           += python_bindings
PYTHON_DEST_DIR   := $(PREFIX)/lib/python$(shell $(PYTHON) -c 'import sys; sys.stdout.write(sys.version[:3]+"\n")')/site-packages
ifeq ($(MINIMUMDEPS), false) 
TEST              += python_test
PYTHONCAT         += $(PYTHONDIR)/FireDeamonMax.py
PYTHONBINDCAT     += $(PYTHONDIR)/max.i
PYTHONBINDCXXCP   := $(PYTHONDIR)/FireDeamon_wrap.cxx.max
PYTHONBINDPYCP    := $(PYTHONDIR)/cpp.py.max
endif
ifeq ($(SWIGRUN), true) 
RUNSWIG           += python_swig
else
RUNSWIG           += python_pseudoswig
endif
endif
#-----------------------------------------------------
#              DEFAULT BUILD RULES
#-----------------------------------------------------
.PHONY : default
default : mainlib bindings
#-----------------------------------------------------
#       INSTALLATION AND UNINSTALLATION RULES 
#-----------------------------------------------------
.PHONY : install 
install : $(INSTALL)

.PHONY : uninstall 
uninstall : $(UNINSTALL)

.PHONY : main_install 
main_install : $(MAINHDR)
	mkdir -p $(PREFIX)/lib
	cp $(LIBDIR)/libFireDeamon.so $(PREFIX)/lib/libFireDeamon.so
	cp $(LIBDIR)/libFireDeamon.a $(PREFIX)/lib/libFireDeamon.a

.PHONY : main_uninstall 
main_uninstall : 
	rm -f $(PREFIX)/lib/libFireDeamon.so $(PREFIX)/lib/libFireDeamon.a
	rm -rf $(PREFIX)/include/FireDeamon

.PHONY : python_install 
python_install :
	@echo "Will install in $(PYTHON_DEST_DIR)"
	mkdir -p $(PYTHON_DEST_DIR)/FireDeamon
	cp $(PYTHONDIR)/FireDeamon.py $(PYTHON_DEST_DIR)/FireDeamon/__init__.py
	cp $(PYTHONDIR)/cpp.py $(PYTHON_DEST_DIR)/FireDeamon/cpp.py
	cp $(PYTHONDIR)/_cpp.so $(PYTHON_DEST_DIR)/FireDeamon/_cpp.so

.PHONY : python_uninstall 
python_uninstall : 
	rm -rf $(PYTHON_DEST_DIR)/FireDeamon
#-----------------------------------------------------
#              BUILD RULES FOR LANGUAGE BINDINGS
#-----------------------------------------------------
.PHONY : bindings
bindings : $(MAINOBJ) $(RUNSWIG)
	@if [ "$(RUNBIND)" ]; then $(MAKE) $(RUNBIND); fi
#-----------------------------------------------------
#              BUILD RULES FOR PYTHON MODULE
#-----------------------------------------------------
python_bindings : $(PYTHONBINDOBJ)
	@echo "Linking shared python library..."
	$(GXX) -shared -fPIC $(PYTHONBINDOBJ) -L$(PYTHONLIB) -L$(LIBDIR) -lFireDeamon $(PYTHONLDFLAGS) -o $(PYTHONDIR)/_cpp.so

.PHONY : python_swig
python_swig : python_cat
	$(SWIG)  -classic -c++ -I$(INCDIR) -I$(SWIGINC) -python -o $(PYTHONDIR)/FireDeamon_wrap.cxx $(PYTHONBINDCATDEST)

.PHONY : python_pseudoswig
python_pseudoswig : python_cat
	cp $(PYTHONBINDCXXCP) $(PYTHONBINDCXX)
	cp $(PYTHONBINDPYCP) $(PYTHONBINDPY)

.PHONY : python_cat
python_cat :
	cat $(PYTHONCAT) > $(PYTHONCATDEST)
	cat $(PYTHONBINDCAT) > $(PYTHONBINDCATDEST)
#-----------------------------------------------------
#              BUILD RULES FOR C++ LIBRARY
#-----------------------------------------------------
.PHONY : mainlib 
mainlib : $(MAINOBJ)
	@echo "Linking shared library..."
	$(GXX) -shared -fPIC $(MAINOBJ) $(OTHERLIB) -L$(CGALLIB) $(LDFLAGS) -Wl,-z,defs -o $(LIBDIR)/libFireDeamon.so
	@echo "Linking static library..."
	$(AR) rcs $(LIBDIR)/libFireDeamon.a $(MAINOBJ)
#-----------------------------------------------------
#              BUILD RULES FOR DOCUMENTATION
#-----------------------------------------------------
.PHONY : doc
doc: python_cat
	@echo "Creating documentation using doxygen..."
	@test "$(MINIMUMDEPS)" = "false" && echo "WARNING: not creating documentation for maximum functionaliry version" >&2
	doxygen doxyfile
	@echo "Compiling pdf using LaTeX..."
	$(MAKE) -C documentation/latex
	@mv documentation/latex/refman.pdf documentation/Documentation.pdf
#-----------------------------------------------------
#              CLEANING RULES
#-----------------------------------------------------
.PHONY : clean
clean : clean_bindings clean_mainlib clean_doc

.PHONY : clean_bindings
clean_bindings : 
	rm -f  $(PYTHONDIR)/*.o $(PYTHONDIR)/cpp.py  $(PYTHONDIR)/*.cxx $(PYTHONDIR)/*.so $(PYTHONCATDEST) $(PYTHONBINDCATDEST)

.PHONY : clean_mainlib
clean_mainlib : 
	rm -f $(MAINOBJ)
	rm -f $(LIBDIR)/*.so $(LIBDIR)/*.a
	rm -f $(TESTDIR)/*.o
	rm -f $(TESTDIR)/test*.exe

.PHONY : clean_doc
clean_doc : 
	rm -rf documentation
#-----------------------------------------------------
#    BUILD TEST PROGRAMME AND RUN TEST PROGRAMMES
#-----------------------------------------------------
.PHONY : test
test: $(TEST)
.PHONY : mainlib_test
mainlib_test : mainlib_test1 mainlib_test2 mainlib_test3 mainlib_test4
.PHONY : mainlib_test1
mainlib_test1 : $(TEST1OBJ)
	$(GXX) $(TEST1OBJ) $(OTHERLIB) -L$(CGALLIB) $(LDFLAGS) -o $(TESTDIR)/test1.exe
	@printf "Executing first test file "
	@test/test1.exe && echo "success" || echo "failed"
.PHONY : mainlib_test2
mainlib_test2 : $(TEST2OBJ) mainlib
	$(GXX) $(TEST2OBJ) $(OTHERLIB) -L$(LIBDIR) -L$(CGALLIB) $(LDFLAGS) -lFireDeamon -o $(TESTDIR)/test2.exe
	@printf "Executing second test file. Have you done 'make install' first? "
	@test/test2.exe && echo "success" || echo "failed"
.PHONY : mainlib_test3
#mainlib_test3 : $(TEST3OBJ) mainlib
mainlib_test3 : $(TEST3OBJ)
	$(GXX) $(TEST3OBJ) $(OTHERLIB) -L$(CGALLIB) $(LDFLAGS) -o $(TESTDIR)/test3.exe
	@printf "Executing third test file "
	@test/test3.exe && echo "success" || echo "failed"
.PHONY : mainlib_test4
mainlib_test4 : $(TEST4OBJ) mainlib
	$(GXX) $(TEST4OBJ) $(OTHERLIB) -L$(LIBDIR) -L$(CGALLIB) $(LDFLAGS) -lFireDeamon -o $(TESTDIR)/test4.exe
	@printf "Executing fourth test file. Have you done 'make install' first? "
	@test/test4.exe && echo "success" || echo "failed"
.PHONY : python_test
python_test : bindings
	@printf "Executing Python test file. Have you done 'make install' first?"
	@bash $(TESTDIR)/test.py.sh $(PYTHONDIR) $(PYTHON) && echo "success" || echo "failed"
#-----------------------------------------------------
#              GENERIC BUILD RULES
#-----------------------------------------------------
$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@printf "Compiling %-25s > %-25s\n" $< $@
	mkdir -p $(dir $@)
	$(GXX)  -std=c++0x -pedantic -Wall -Wextra -c -fPIC $(CPPFLAGS) $(OTHERINC) -I$(CGALINC) -I$(BOOSTINC) -I$(INCDIR) $< -o $@

# Make compilation rules for cxx files of bindings 
$(PYTHONDIR)/%.o : $(PYTHONDIR)/%.cxx
	@printf "Compiling %-25s > %-25s\n" $< $@
	mkdir -p $(dir $@)
	$(GXX) -std=c++0x -pedantic -Wall -Wextra -c -fPIC $(CPPFLAGS) $(OTHERINC) -I$(PYTHONINC) -I$(INCDIR) $< -o $@

# Make compilation rules for cpp files for test executables
$(TESTDIR)/%.o : $(TESTDIR)/%.cpp
	@printf "Compiling %-25s > %-25s\n" $< $@
	mkdir -p $(dir $@)
	$(GXX) -std=c++0x -pedantic -Wall -Wextra -c -fPIC $(CPPFLAGS) $(OTHERINC) -I$(CGALINC) -I$(BOOSTINC) -I$(INCDIR) $< -o $@
#-----------------------------------------------------
#              HEADER INSTALLATION RULES
#-----------------------------------------------------
$(PREFIX)/include/FireDeamon/%.h : $(INCDIR)/FireDeamon/%.h
	mkdir -p $(dir $@)
	cp $< $@
