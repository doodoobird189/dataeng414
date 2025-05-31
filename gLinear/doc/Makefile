
all:	coding_convention.ps programmer.ps



clean:
	for file in *~ *.aux core *.dvi *.ps; do \
		if [ -f $${file} ]; then \
			rm $${file}; \
		fi; \
	done;
	
	
%.dvi:	%.tex
	latex $< %% latex $<

%.ps:	%.dvi
	dvips -t A4 -d 600 -o $@ $<
