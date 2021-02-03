
DOCS = hermes-2-manual.pdf

.PHONY: all

all: $(DOCS)

%.pdf: %.tex 
	pdflatex $< 
	pdflatex $<

