ANALYSIS := QCD_INSTANTON_KKS
PLUGIN := Rivet$(ANALYSIS).so
SRC := $(ANALYSIS).cc

.PHONY: all clean

all: $(PLUGIN)

$(PLUGIN): $(SRC) $(ANALYSIS).info $(ANALYSIS).plot
	rivet-build $(PLUGIN) $(SRC)

clean:
	rm -f $(PLUGIN)
