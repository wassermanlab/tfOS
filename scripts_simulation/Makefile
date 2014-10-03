CC 		= g++
OBJS 		= sim_initiate.o ProteinBindingSimulator.o ProteinInteractions.o SequenceAccessibility.o BindingProtein.o
EXECUTABLE	= sim_initiate
CFLAGS		= -W -g

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@

clean:
	rm -f $(OBJS) $(EXECUTABLE)

.SUFFIXES: .o .cpp

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@
