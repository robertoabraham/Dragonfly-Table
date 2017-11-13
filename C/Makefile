CC = clang
CFLAGS = -std=c99 -g

MANAGER?=homebrew
ifeq ($(MANAGER),homebrew)
        LFLAGS = -L /usr/local/lib -lcfitsio -lgsl -lgslcblas 
        INCDIR = /usr/local/include
else
        LFLAGS = -L /opt/local/lib -lcfitsio -lgsl -lgslcblas 
        INCDIR = /opt/local/include
endif

DEPS = 
OBJ = 
PROGRAMS = tread tfitdist tfitpoly tfitsurf tablist tlowess

%.o: %.c $(DEPS)
	$(CC) -c $(CFLAGS) -I${INCDIR} -o $@ $< 

all: tread tfitdist tfitpoly tlowess

tread: tread.c table.o
	$(CC) -o $@ $^ -I$(INCDIR) $(FFLAGS) ${LFLAGS}

tfitdist: tfitdist.c table.o expfit.o gaussfit.o
	$(CC) -o $@ $^ -I$(INCDIR) $(FFLAGS) ${LFLAGS}

tfitpoly: tfitpoly.c table.o 
	$(CC) -o $@ $^ -I$(INCDIR) $(FFLAGS) ${LFLAGS}

tlowess: tlowess.c table.o 
	$(CC) -o $@ $^ -I$(INCDIR) $(FFLAGS) ${LFLAGS}

tfitsurf: tfitsurf.c table.o 
	$(CC) -o $@ $^ -I$(INCDIR) $(FFLAGS) ${LFLAGS}

tablist: tablist.c  
	$(CC) -o $@ $^ -I$(INCDIR) $(FFLAGS) ${LFLAGS}

clean:
	rm -f *.o

install: $(PROGRAMS)
	mv $(PROGRAMS) /usr/local/bin
