TARGET=efficiency

FEDRALIBS := -lEIO -lEdb -lEbase -lEdr -lScan -lAlignment -lEmath -lEphys -lvt -lDataConversion

$(TARGET): $(TARGET).cpp
	g++ $(TARGET).cpp -w `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $(TARGET)

clean:
	$(RM) $(TARGET)
