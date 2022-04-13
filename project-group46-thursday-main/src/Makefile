# this makefile is intended for g++ on Linux

CC = g++
CFLAGS = -c -Wall -O3
# CFLAGS = -c
LDFLAGS = -pthread
INCLUDE = -I../include/
SOURCES = project.cpp iofunc.cpp filter.cpp fourier.cpp genfunc.cpp logfunc.cpp \
audiofunc.cpp mode.cpp mono_block.cpp RF_front_end.cpp stereo_Processing.cpp \
RDS.cpp RDSutility.cpp 
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = project

all: $(EXECUTABLE) clean

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDE) $< -o $@

clean:
	-rm $(OBJECTS)
