all: plot.png

main: main.cpp
	g++ -o main main.cpp

numSol.dat analSol.dat ym.dat ymAnal.dat &: main
	for f in tmp/*; do rm "$$f"; echo "rm $$f"; done
	./main

plot.png: numSol.dat analSol.dat plotting.gpi
	gnuplot plotting.gpi

clean:
	main numSol.dat analSol.dat ym.dat ymAnal.dat plot.png
