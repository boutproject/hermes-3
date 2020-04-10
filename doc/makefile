# Some distros may provide sphinx-build-3 and sphinx-build-2 instead
# of sphinx-build. An automatic check would be helpful, but at least
# with this, it should be easy to change:
sphinx-build?=$(shell which sphinx-build-3 >/dev/null 2>&1 && echo sphinx-build-3 || echo sphinx-build)

all: sphinx
manual: all
# set some shorter names
pdf: sphinx-pdf
html: sphinx-html
man: sphinx-man
sphinx: sphinx-pdf

sphinx-pdf: doxygen
	$(sphinx-build) -b latex sphinx/ build/
	cd build && latexmk -pdf hermes3
	test -e hermes3.pdf || ln -s build/hermes3.pdf .
	@echo "Documentation is available in `pwd`/hermes3.pdf"

sphinx-html: doxygen
	$(sphinx-build) -b html sphinx/ html/
	@echo "Documentation available in file://`pwd`/html/index.html"

sphinx-man: doxygen
	$(sphinx-build) -b man sphinx/ man/
	@echo "Documentation available in `pwd`/man/hermes3.1"

# Run doxygen, ignore if it fails (leading '-')
doxygen:
	-cd doxygen && doxygen Doxyfile

# Run breathe-apidoc, ignore if it fails (leading '-')
breathe-autogen: doxygen
	-breathe-apidoc -f -o sphinx/_breathe_autogen doxygen/hermes3/xml

clean:
	@echo "Cleaning up..."
	@$(RM) -f *.pdf
	@$(RM) -rf html/
	@$(RM) -rf build/
	@$(RM) -rf man/
	@$(RM) -rf sphinx/_breathe_autogen/
