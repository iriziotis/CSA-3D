# Minimal make file for Sphinx documentation

.PHONY: help clean html clean-html

SPHINXBUILD   ?= sphinx-build
BUILDDIR      ?= _build
SOURCEDIR     = .

help:
	@echo "Sphinx documentation targets:"
	@echo "  make html       - build HTML documentation"
	@echo "  make clean      - remove build directory"
	@echo "  make clean-html - clean HTML and rebuild"

clean:
	rm -rf $(BUILDDIR)

html:
	$(SPHINXBUILD) -b html $(SOURCEDIR) $(BUILDDIR)/html
	@echo "HTML documentation built in $(BUILDDIR)/html/"

clean-html: clean html
	@echo "Clean rebuild complete"
