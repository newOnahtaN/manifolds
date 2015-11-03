test: test.cpp
	g++ test.cpp -o test -O2 -larmadillo

fbeuler_jb: fbeuler_jb.cpp
	g++ fbeuler_jb.cpp -o fbeuler_jb -O2 -larmadillo

clean:
	rm test fbeuler_jb
