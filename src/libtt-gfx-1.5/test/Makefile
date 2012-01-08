CXX = g++
CFLAGS = -g -O2 -I../include

all: test1.exe test2.exe

test1.exe: app.cxx test1.cxx ../src/libtt-gfx.a
	$(CXX) $(CFLAGS) -o $@ $^

test2.exe: app.cxx test2.cxx ../src/libtt-gfx.a
	$(CXX) $(CFLAGS) -o $@ $^

clean:
	rm -f *.exe

