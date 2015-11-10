test: test.cpp
	g++ test.cpp -o test -O2 -larmadillo

all: fbeuler_jb fbeuler_gs fbtrap_jb

fbeuler_jb: fbeuler_jb.cpp
	g++ fbeuler_jb.cpp -o fbeuler_jb -O2 -larmadillo

fbeuler_gs: fbeuler_gs.cpp
	g++ fbeuler_gs.cpp -o fbeuler_gs -O2 -larmadillo

fbtrap_jb: fbtrap_jb.cpp
	g++ fbtrap_jb.cpp -o fbtrap_jb -O2 -larmadillo

clean:
	rm test fbeuler_jb
