IT = ata_squint_memo
TEXARGS = -interaction nonstopmode -halt-on-error -file-line-error
TEXDEPS = ATA.bib

default: $(IT).pdf

%.pdf: %.tex $(TEXDEPS)
	pdflatex $(TEXARGS) $(basename $<)
	bibtex $(basename $<)
	pdflatex $(TEXARGS) $(basename $<)
	pdflatex $(TEXARGS) $(basename $<)

clean:
	-rm -f $(IT).aux $(IT).bbl $(IT).blg $(IT).log $(IT).out $(IT).pdf $(IT).toc
