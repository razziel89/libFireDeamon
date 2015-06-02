# makefile for libFireDeamon

#inport build parameters
include make.vars

#-----------------------------------------------------
#         SOURCE FOR MAIN C++ LIBRARY
#-----------------------------------------------------
MAINSRC  := $(SRCDIR)/skin_surface_deamon.cpp
MAINOBJ  := $(OBJDIR)/skin_surface_deamon.o
#-----------------------------------------------------
#     DEFAULT VARIABLES WITHOUT BINDINGS
#-----------------------------------------------------
SWIGRUN   :=
BINDRUN   :=
INSTALL   := main_install
UNINSTALL := main_uninstall
#-----------------------------------------------------
#    VARIABLE DEFINITIONS IF PYTHON IS REQUESTED
#-----------------------------------------------------
ifeq ($(USEPYTHON), true) 
PYTHONBINDSRC := $(PYTHONDIR)/FireDeamon_wrap.cxx
PYTHONBINDOBJ := $(PYTHONDIR)/FireDeamon_wrap.o
SWIGRUN   += python_swig
BINDRUN   += python_bindings
INSTALL   += python_install
UNINSTALL += python_uninstall
PYTHON_DEST_DIR := $(PREFIX)/lib/python$(shell $(PYTHON) -c 'import sys; print sys.version[:3]')/site-packages
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
main_install :
	mkdir -p $(PREFIX)/lib
	@mkdir -p $(PREFIX)/lib
	cp $(LIBDIR)/libFireDeamon.so $(PREFIX)/lib/libFireDeamon.so
	@cp $(LIBDIR)/libFireDeamon.so $(PREFIX)/lib/libFireDeamon.so
	cp $(LIBDIR)/libFireDeamon.a $(PREFIX)/lib/libFireDeamon.a
	@cp $(LIBDIR)/libFireDeamon.a $(PREFIX)/lib/libFireDeamon.a

.PHONY : main_uninstall 
main_uninstall : clean_mainlib
	rm -f $(PREFIX)/lib/libFireDeamon.so $(PREFIX)/lib/libFireDeamon.a
	@rm -f $(PREFIX)/lib/libFireDeamon.so $(PREFIX)/lib/libFireDeamon.a


.PHONY : python_install 
python_install :
	@echo "Will install in $(PYTHON_DEST_DIR)"
	mkdir -p $(PYTHON_DEST_DIR)
	@mkdir -p $(PYTHON_DEST_DIR)
	cp $(PYTHONDIR)/FireDeamon.py $(PYTHON_DEST_DIR)/FireDeamon.py
	@cp $(PYTHONDIR)/FireDeamon.py $(PYTHON_DEST_DIR)/FireDeamon.py
	cp $(PYTHONDIR)/_FireDeamon.so $(PYTHON_DEST_DIR)/_FireDeamon.so
	@cp $(PYTHONDIR)/_FireDeamon.so $(PYTHON_DEST_DIR)/_FireDeamon.so

.PHONY : python_uninstall 
python_uninstall :
	rm -rf $(PYTHON_DEST_DIR)/FireDeamon.py $(PYTHON_DEST_DIR)/_FireDeamon.so
	@rm -rf $(PYTHON_DEST_DIR)/FireDeamon.py $(PYTHON_DEST_DIR)/_FireDeamon.so
#-----------------------------------------------------
#              BUILD RULES FOR LANGUAGE BINDINGS
#-----------------------------------------------------
.PHONY : bindings
bindings : $(MAINOBJ) $(SWIGRUN)
	@if [[ "$(BINDRUN)" ]]; then $(MAKE) $(BINDRUN); fi
#-----------------------------------------------------
#              BUILD RULES FOR PYTHON MODULE
#-----------------------------------------------------
python_bindings : $(PYTHONBINDOBJ)
	@echo "Linking shared python library..."
	$(GXX) -shared -fPIC -I$(CGALINC) -I$(PYTHONINC) -I$(INCDIR) -L$(PYTHONLIB) -L$(CGALLIB) \
		$(LDFLAGS) $(PYTHONLDFLAGS) -o $(PYTHONDIR)/_FireDeamon.so $(MAINOBJ) $(PYTHONBINDOBJ)
	@$(GXX) -shared -fPIC -I$(CGALINC) -I$(PYTHONINC) -I$(INCDIR) -L$(PYTHONLIB) -L$(CGALLIB) \
		$(LDFLAGS) $(PYTHONLDFLAGS) -o $(PYTHONDIR)/_FireDeamon.so $(MAINOBJ) $(PYTHONBINDOBJ)

.PHONY : python_swig
python_swig :
	$(SWIG)  -classic -c++ -I$(INCDIR) -I$(SWIGINC) -python -o $(PYTHONDIR)/FireDeamon_wrap.cxx $(PYTHONDIR)/FireDeamon.i
	@$(SWIG) -classic -c++ -I$(INCDIR) -I$(SWIGINC) -python -o $(PYTHONDIR)/FireDeamon_wrap.cxx $(PYTHONDIR)/FireDeamon.i
#-----------------------------------------------------
#              BUILD RULES FOR C++ LIBRARY
#-----------------------------------------------------
.PHONY : mainlib 
mainlib : $(MAINOBJ)
	@echo "Linking shared library..."
	$(GXX) -shared -fPIC -I$(CGALINC) -I$(INCDIR) -L$(CGALLIB) $(LDFLAGS) -o $(LIBDIR)/libFireDeamon.so $(MAINOBJ)
	@$(GXX) -shared -fPIC -I$(CGALINC) -I$(INCDIR) -L$(CGALLIB) $(LDFLAGS) -o $(LIBDIR)/libFireDeamon.so $(MAINOBJ)
	@echo "Linking static library..."
	$(AR) rcs $(LIBDIR)/libFireDeamon.a $(MAINOBJ)
	@$(AR) rcs $(LIBDIR)/libFireDeamon.a $(MAINOBJ)
#-----------------------------------------------------
#              CLEANING RULES
#-----------------------------------------------------
.PHONY : clean
clean : clean_bindings clean_mainlib

.PHONY : clean_bindings
clean_bindings : 
	rm -f  $(PYTHONDIR)/*.o $(PYTHONDIR)/*.d.tmp $(PYTHONDIR)/*.py  $(PYTHONDIR)/*.cxx $(PYTHONDIR)/*.so
	@rm -f $(PYTHONDIR)/*.o $(PYTHONDIR)/*.d.tmp $(PYTHONDIR)/*.py  $(PYTHONDIR)/*.cxx $(PYTHONDIR)/*.so

.PHONY : clean_mainlib
clean_mainlib : 
	rm -f  $(OBJDIR)/*.o $(OBJDIR)/*.d.tmp
	@rm -f $(OBJDIR)/*.o $(OBJDIR)/*.d.tmp
	rm -f  $(LIBDIR)/*.so $(LIBDIR)/*.a
	@rm -f $(LIBDIR)/*.so $(LIBDIR)/*.a
#-----------------------------------------------------
#              GENERIC BUILD RULES
#-----------------------------------------------------
# Make compilation rules for cpp files of main library
$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@printf "Compiling %-25s > %-25s\n" $< $@
	@mkdir -p $(dir $@)
	$(GXX)	-c -fPIC -I$(CGALINC) -I$(INCDIR) $< -o $@
	@$(GXX) -c -fPIC -I$(CGALINC) -I$(INCDIR) $< -o $@

# Make compilation rules for cxx files of bindings 
$(PYTHONDIR)/%.o : $(PYTHONDIR)/%.cxx
	@printf "Compiling %-25s > %-25s\n" $< $@
	@mkdir -p $(dir $@)
	$(GXX)	-c -fPIC -I$(PYTHONINC) -I$(INCDIR) $< -o $@
	@$(GXX) -c -fPIC -I$(PYTHONINC) -I$(INCDIR) $< -o $@
