/*
 * Conway's Game of Life
 *
 * A. Mucherino
 *
 * PPAR, TP5
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
//#include <mpi.h>

#define WORLD_SIZE 32

static int N = WORLD_SIZE;
static int M = WORLD_SIZE/(sizeof(unsigned int)*4);
static int itMax = 20;
typedef struct coord
{
	int k;
	int l;
} coord;

// allocation only
unsigned int* allocate()
{
	return (unsigned int*)calloc(N*M,sizeof(unsigned int));
}

// conversion cell location : 2d --> 1d
// (row by row)
coord code(int x, int y, int dx, int dy)
{
	int i = (x + dx)%N;
	int j = (y + dy)%N;
	if (i < 0)  i = N + i;
	if (j < 0)  j = N + j;
	coord c = {i*M + (j >> (sizeof(unsigned int))), j%(sizeof(unsigned int)*4)};
	return c;
}

// writing into a cell location 
void write_cell(int x,int y,unsigned int value,unsigned int *world)
{
	coord c;
	c = code(x,y,0,0);
	world[c.k] &= ~(0b11 << (c.l*2));
	world[c.k] |= (value << (c.l*2));
}

// random generation
unsigned int* initialize_random()
{
	int x,y;
	unsigned int cell;
	unsigned int *world;

	world = allocate();
	for (x = 0; x < N; x++)
	{
		for (y = 0; y < N; y++)
		{
			if (rand()%5 != 0)
				cell = 0;
			else if (rand()%2 == 0)
				cell = 1;
			else
				cell = 2;
			write_cell(x,y,cell,world);
		}
	}
	return world;
}

// dummy generation
unsigned int* initialize_dummy()
{
	int x,y;
	unsigned int *world;

	world = allocate();
	for (x = 0; x < N; x++)
		for (y = 0; y < N; y++)
			write_cell(x,y,x%3,world);
	return world;
}

// "glider" generation
unsigned int* initialize_glider()
{
	int x,y,mx,my;
	unsigned int *world;

	world = allocate();
	for (x = 0; x < N; x++)
		for (y = 0; y < N; y++)
			write_cell(x,y,0,world);

	mx = N/2 - 1; my = N/2 - 1;
	x = mx; y = my + 1;  write_cell(x,y,1,world);
	x = mx + 1; y = my + 2;  write_cell(x,y,1,world);
	x = mx + 2; y = my; write_cell(x,y,1,world);
				y = my + 1; write_cell(x,y,1,world);
				y = my + 2; write_cell(x,y,1,world);

	return world;
}

// "small exploder" generation
unsigned int* initialize_small_exploder()
{
	int x,y,mx,my;
	unsigned int *world;

	world = allocate();
	for (x = 0; x < N; x++)
		for (y = 0; y < N; y++)
			write_cell(x,y,0,world);

	mx = N/2 - 2; my = N/2 - 2;
	x = mx; y = my + 1; write_cell(x,y,2,world);
	x = mx + 1; y = my; write_cell(x,y,2,world);
				y = my + 1; write_cell(x,y,2,world);
				y = my + 2; write_cell(x,y,2,world);
	x = mx + 2; y = my; write_cell(x,y,2,world);
				y = my + 2; write_cell(x,y,2,world);
	x = mx + 3; y = my + 1; write_cell(x,y,2,world);

	return world;
}


// reading a cell
int read_cell(int x,int y,int dx,int dy,unsigned int *world)
{
	coord c = code(x,y,dx,dy);
	return ((world[c.k] >> (c.l*2)) & 0b11);
}

// updating counters
void update(int x,int y,int dx,int dy,unsigned int *world,int *nn,int *n1,int *n2)
{
	unsigned int cell = read_cell(x,y,dx,dy,world);
	if (cell != 0)
	{
		(*nn)++;
		if (cell == 1)
			(*n1)++;
		else
			(*n2)++;
	}
}

// looking around the cell
void neighbors(int x,int y,unsigned int *world,int *nn,int *n1,int *n2)
{
	int dx,dy;

	(*nn) = 0; (*n1) = 0; (*n2) = 0;

	// same line
	dx = -1; dy = 0; update(x,y,dx,dy,world,nn,n1,n2);
	dx = +1; dy = 0; update(x,y,dx,dy,world,nn,n1,n2);

	// one line down
	dx = -1; dy = +1; update(x,y,dx,dy,world,nn,n1,n2);
	dx =  0; dy = +1; update(x,y,dx,dy,world,nn,n1,n2);
	dx = +1; dy = +1; update(x,y,dx,dy,world,nn,n1,n2);	

	// one line up
	dx = -1; dy = -1; update(x,y,dx,dy,world,nn,n1,n2);
	dx =  0; dy = -1; update(x,y,dx,dy,world,nn,n1,n2);
	dx = +1; dy = -1; update(x,y,dx,dy,world,nn,n1,n2);
}

// computing a new generation
short newgeneration(unsigned int *world1,unsigned int *world2,int xstart,int xend)
{
	int x,y;
	int nn,n1,n2;
	unsigned int cell;
	short change = 0;

	// cleaning destination world
	for (x = 0; x < N; x++)
		for (y = 0; y < N; y++)
			write_cell(x,y,0,world2);

	// generating the new world
	for (x = xstart; x < xend; x++)
	{
		for (y = 0; y < N; y++)
		{
			cell = read_cell(x, y, 0, 0, world1);
			neighbors(x, y, world1, &nn, &n1, &n2);
			
			// Live cell
			if(cell != 0)
			{
				// Right number of neighbors, the cell lives on
				if((nn > 1) && (nn < 4))
					write_cell(x, y, cell, world2);
				// Else the cell dies, by over- or under-population (the cell is already set to 0)
				else
					change++;
			}
			// Reproduction, a new cell is born
			else if(nn == 3)
			{
				// Takes on the dominant genus
				if(n1 >= 2)
					write_cell(x, y, 1, world2);
				else
					write_cell(x, y, 2, world2);
				change++;
			}
			// In every other cases the cell remains empty
		}
	}
	return change;
}

// cleaning the screen
void cls()
{
	int i;
	for (i = 0; i < 10; i++)
		fprintf(stdout,"\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
}

// diplaying the world
void print(unsigned int *world)
{
	int i, j, n;
	cls();
	for (i = 0; i < N; i++)
		fprintf(stdout,"-");
	fprintf(stdout,"\n");
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			n = read_cell(i, j, 0, 0, world);
			if (n == 0)
				fprintf(stdout," ");
			if (n == 1)
				fprintf(stdout,"o");
			if (n == 2)
				fprintf(stdout,"x");
		}
		fprintf(stdout,"\n");
	}
	fprintf(stdout,"\n");

	for (i = 0; i < N; i++)
		fprintf(stdout,"-");
	fprintf(stdout,"\n");
	sleep(1);
}

// main
int main(int argc,char *argv[])
{
	int it,change;
	unsigned int *world1,*world2;
	unsigned int *worldaux;

	// getting started  
	//world1 = initialize_dummy();
	//world1 = initialize_random();
	//world1 = initialize_glider();
	world1 = initialize_small_exploder();
	world2 = allocate();
	print(world1);

	it = 0;  change = 1; 
	while (change && it < itMax)
	{
		change = newgeneration(world1,world2,0,N);
		worldaux = world1; world1 = world2; world2 = worldaux;
		print(world1);
		it++;
	}

	// ending
	free(world1); free(world2);
	return 0;
}

