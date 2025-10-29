all: plot.png

main: main.cpp
	g++ -o main main.cpp

numSol.dat analSol.dat &: main
	./main

plot.png: numSol.dat analSol.dat plotting.gpi
	gnuplot plotting.gpi

clean:
	main numSol.dat analSol.dat plot.png
