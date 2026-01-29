all: whale_measurements.out whale_shape.out
%.out: %.R
	R --no-save < $< > $@
clean:
	rm *~ *.out
install:
	cp whale_measurements.rda ../data
	
