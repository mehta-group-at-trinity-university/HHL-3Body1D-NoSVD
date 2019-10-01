
coup_plot.x: coup_plot.o
	ifort coup_plot.o -o coup_plot.x

coup_plot.o: coup_plot.f
	ifort -c coup_plot.f	
