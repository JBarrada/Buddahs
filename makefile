CFLAGS=
LIBS=-Iinclude

OBJFILES = main.o BMP.o

all: clean main.exe

%.o: %.cpp
	g++ -c $(CFLAGS) $(LIBS) $< -o $@

main.exe: $(OBJFILES)
	g++ $(OBJFILES) $(LIBS) -o $@

clean:
	del *.o *.exe