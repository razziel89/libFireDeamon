# makefile for libFireDeamon

default: clean
	$(MAKE) cppdoc pydoc

# ===== C++ documetnation =====
.PHONY : cppdoc
cppdoc:
	@which doxygen
	@echo "Creating documentation using doxygen..."
	doxygen Doxyfile

# ===== distribution rules =====

.PHONY : dist
dist: clean_dist
	@echo "Creating source distribution"
	@python setup.py sdist

.PHONY : publish
publish: dist
	@echo "Publishing to PyPI"
	@python -m twine upload dist/*

# ===== clean rules =====

.PHONY : clean
clean : clean_doc clean_dist

.PHONY : clean_doc
clean_doc : 
	rm -rf docs

.PHONY : clean_dist
clean_dist: 
	rm -rf build dist

# ===== Python documetnation =====
# Sphinx documentation
# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = FireDeamon
SOURCEDIR     = .
BUILDDIR      = docs/tmp_python

spinx-help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
pydoc: Makefile
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	mv -t docs/ $(BUILDDIR)/html/*
	rm -r $(BUILDDIR)
