
EXEC   = hydrostatic_column

CC       =  gcc
OPTIMIZE =  -O2

OPTIONS =  $(OPTIMIZE) $(OPT)


OBJS   = main.o hydrostatic_column.o block_decomposition.o

INCL   = hydrostatic_column.h block_decomposition.h

CFLAGS = $(OPTIONS)

LIBS   =  -lm

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC)
