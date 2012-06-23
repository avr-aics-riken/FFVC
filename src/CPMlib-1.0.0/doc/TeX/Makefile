PDFS=manual.pdf 


.PHONY: all clean rename
.SUFFIXES: .pdf .dvi .tex 

all: $(PDFS) 

.tex.dvi:
	platex $*; \
	platex $*; \
	latex_count=5 ; \
	while egrep -s 'Rerun (LaTeX|to get cross-references right)' $*.log && [ $$latex_count -gt 0 ] ;\
	    do \
	    echo "Rerunning latex...." ;\
	    platex -kanji=utf8 $* ;\
	    latex_count=`expr $$latex_count - 1` ;\
	    done

.dvi.pdf:
	dvipdfmx $*
	cp $(PDFS) ..



clean:
	-rm -f *.pdf *.dvi *.log *.aux *.toc *.lof *.out section/*.aux
relclean:
	-rm  -f *.dvi *.log *.aux *.toc *.lot *.out section/*.aux

program_description_figures=directory_parse_mode.eps \
 parse_dependence_value.eps\
 leaf_parse_mode.eps program_structure.eps \
 parse_condition.eps vector_parse_mode.eps \
 parse_dependence.eps vector_termination_parse_mode.eps

manual.pdf:manual.dvi
manual.dvi:manual.tex
