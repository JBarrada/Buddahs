CFLAGS=-Ofast
LIBS=-Iinclude -pthread

OBJFILES = main.o BMP.o

all: clean main.exe

%.o: %.cpp
	C:\Program Files\mingw-w64\x86_64-7.1.0-posix-seh-rt_v5-rev1\mingw64\bin\g++ -c $(CFLAGS) $(LIBS) $< -o $@

main.exe: $(OBJFILES)
	C:\Program Files\mingw-w64\x86_64-7.1.0-posix-seh-rt_v5-rev1\mingw64\bin\g++ $(OBJFILES) $(LIBS) -o $@

clean:
	del *.o *.exe