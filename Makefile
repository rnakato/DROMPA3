.SUFFIXES: .o .c

CC = gcc
TARGET = parse2wig drompa_peakcall drompa_draw

OBJS_UTIL = alloc.o readfile.o filehandle.o seq.o stringp.o argument.o outputfile.o compress.o 
HEADS_UTIL = alloc.h readfile.h filehandle.h seq.h stringp.h argument.h outputfile.h compress.h macro.h common.h memtest.h 

OBJS_PW = pw_main.o pw_init.o pw_param_new.o pw_readmapfile.o pw_filtering.o pw_makefile.o pw_complexity.o pw_gc.o readfasta.o pw_estimate.o pw_ccp.o $(OBJS_UTIL)
HEADS_PW = pw_init.h pw_filtering.h pw_param_new.h pw_readmapfile.h pw_gv.h pw_makefile.h pw_complexity.h pw_gc.h readfasta.h pw_estimate.h pw_ccp.h $(HEADS_UTIL)

OBJS_DR = drompa_param_new.o drompa_readfile.o dp_call.o gsl_func.o drompa_usage.o $(OBJS_UTIL)
HEADS_DR = drompa_gv.h drompa_param_new.h drompa_readfile.h dp_call.h gsl_func.h drompa_usage.h $(HEADS_UTIL)
OBJS_DP = dp_main.o dp_init.o $(OBJS_DR)
HEADS_DP = dp_init.h $(HEADS_DR)
OBJS_DD = dd_main.o dd_init.o dd_readannotation.o dd_stroke.o dd_profile.o dd_heatmap.o dd_otherfunc.o $(OBJS_DR)
HEADS_DD = dd_init.h dd_readannotation.h dd_gv.h dd_stroke.h dd_profile.h dd_heatmap.h dd_otherfunc.h $(HEADS_DR)

CFLAGS += -Wall -W -O3 $(FLAG)
LIBS += -lm -lz -lgsl -lgslcblas

ifdef READSV
DEBUG=1
CFLAGS += -DREADSV
endif
ifdef DEBUG
CFLAGS += -g -DDEBUG
endif
ifdef MEMTEST
CFLAGS += -DMEMTEST
LIBS += -L./lib -ldbg_malloc
endif
ifdef CLOCK
CFLAGS += -DCLOCK
endif
ifdef SHOW_WIGFILE
CFLAGS += -DSHOW_WIGFILE
endif

all: $(TARGET)

echo:
	@echo "CFLAGS = $(CFLAGS)"
	@echo "LIBS = $(LIBS)"

parse2wig: $(OBJS_PW)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
drompa_peakcall: $(OBJS_DP)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
drompa_draw: $(OBJS_DD)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) `pkg-config gtk+-2.0 --libs` 

dd_heatmap.o: dd_heatmap.c
	$(CC) -c $< $(CFLAGS) `pkg-config gtk+-2.0 --cflags`
dd_stroke.o: dd_stroke.c
	$(CC) -c $< $(CFLAGS) `pkg-config gtk+-2.0 --cflags`
.c.o:
	$(CC) -c $< $(CFLAGS) 

clean:
	rm -f $(OBJS_PW) $(OBJS_DP) $(OBJS_DD) $(TARGET)
	rm *~

$(OBJS_PW): $(HEADS_PW) Makefile
$(OBJS_DP): $(HEADS_DP) Makefile
$(OBJS_DD): $(HEADS_DD) Makefile
