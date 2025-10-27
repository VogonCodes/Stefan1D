all: plot.png

main: main.cpp
	g++ -o main main.cpp

data.dat anal.dat ym.dat &: main
	./main

plot.png stefan.png ym.png: data.dat anal.dat ym.dat plotting.gpi
	gnuplot plotting.gpi

clean:
	rm -f main data.dat anal.dat ym.dat plot.png stefan.png ym.png
