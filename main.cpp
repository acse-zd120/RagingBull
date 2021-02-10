#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <omp.h>

using namespace std;

//Note that this is a serial implementation with a periodic grid
vector<vector<bool>> grid, new_grid;
int imax, jmax;
int max_steps = 100;

int num_neighbours(int ii, int jj)
{
	int ix, jx;
	int cnt = 0;
	for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++)
			if (i != 0 || j != 0)
			{
				ix = (i + ii + imax) % imax;
				jx = (j + jj + jmax) % jmax;
				if (grid[ix][jx]) cnt++;
			}
	return cnt;
}

//Overload the function in order to calculate the neighbours of divided parts
int num_neighbours(int ii, int jj, int xmax, int ymax, vector<vector<bool>> mat)
{
	int ix, jx;
	int cnt = 0;
	for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++)
			if (i != 0 || j != 0)
			{
				ix = (i + ii + xmax) % xmax;
				jx = (j + jj + ymax) % ymax;
				if (mat[ix][jx]) cnt++;
			}
	return cnt;
}

void grid_to_file(int it)
{
	stringstream fname;
	fstream f1;
	fname << "output" << "_" << it << ".dat";
	f1.open(fname.str().c_str(), ios_base::out);
	for (int i = 0; i < imax; i++)
	{
		for (int j = 0; j < jmax; j++)
			f1 << grid[i][j] << "\t";
		f1 << endl;
	}
	f1.close();
}

void do_iteration(void)
{  

#pragma omp parallel for num_threads(2)  
	for (int i = 0; i < imax; i++)
		for (int j = 0; j < jmax; j++)
		{
			new_grid[i][j] = grid[i][j];
			int num_n = num_neighbours(i, j);
			if (grid[i][j])
			{
				if (num_n != 2 && num_n != 3)
					new_grid[i][j] = false;
			}
			else if (num_n == 3) new_grid[i][j] = true;
		}
	grid.swap(new_grid);
}

//Overload the iteration fuction in order to iterate the divided parts
void do_iteration(int xmax, int ymax, vector<vector<bool>>& mat)
{
	//initialize a new matrix storing the new grid
	vector<vector<bool>> new_grid1;
	new_grid1.resize(xmax, vector<bool>(ymax));
	
#pragma omp parallel for num_threads(4)  
	for (int i = 0; i < xmax; i++)
		for (int j = 0; j < ymax; j++)
		{
			new_grid1[i][j] = mat[i][j];
			int num_n = num_neighbours(i, j, xmax, ymax, mat);
			if (mat[i][j])
			{
				if (num_n != 2 && num_n != 3)
					new_grid1[i][j] = false;
			}
			else if (num_n == 3) new_grid1[i][j] = true;
		}
	mat.swap(new_grid1);
}

void do_split_method(int n)
{
	vector<vector<vector<bool>>> mat(n);

	//first set the size of two-dimension matrix
	for (int i = 0; i < n; i++)
	{
		mat[i].resize(imax);
	}

	//first set the size of three-dimension matrix
	//or the memory will go wrong 
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < imax; j++)
		{
			mat[i][j].resize(jmax);
		}
	}

	//divide the grid into n parts
	for (int p = 0; p < n; p++)
	{
		if (p == 0)
		{
			for (int i = 0; i < imax; i++)
			{
				mat[p][i].push_back(grid[i][jmax-2]);
				mat[p][i].push_back(grid[i][jmax-1]);

				for (int j = 0; j < jmax / n + 2; j++)
				{
					mat[p][i].push_back(grid[i][j]);
					cout << mat[p][i][j];
				}
			}
		}
		else if (p == n - 1)
		{
			cout << endl;
			for (int i = 0; i < imax; i++)
			{
				for (int j = jmax * p / n - 2; j < jmax; j++)
				{
					mat[p][i].push_back(grid[i][j]);
					cout << mat[p][i][j];
				}
				mat[p][i].push_back(grid[i][0]);
				mat[p][i].push_back(grid[i][1]);
			}
		}
		else
		{
			cout << endl;
			for (int i = 0; i < imax; i++)
				for (int j = jmax * p / n - 2; j < jmax * (p + 1) / n + 2; j++)
				{
					mat[p][i].push_back(grid[i][j]);
					cout << mat[p][i][j];
				}
		}
	}

	//using fuction to iterate all parts 
	for (int iter=0; iter<2;iter++)
#pragma omp parallel for num_threads(4) 
		for (int p = 0; p < n; p++)
			{
				do_iteration(imax, jmax / n + 4, mat[p]);
			}

	//merge all parts together
	for (int i = 0; i < imax; i++)
		for (int j = 0; j < jmax; j++)
		{
			int p = j / (jmax / n) - 1;
			grid[i][j] = mat[p][i][j+2];
		}
}

int main(int argc, char *argv[])
{
	srand(time(NULL));
	imax = 100;
	jmax = 100;
	grid.resize(imax, vector<bool>(jmax));
	new_grid.resize(imax, vector<bool>(jmax));

	//set an initial random collection of points - You could set an initial pattern
	for (int i = 0; i < imax; i++)
		for (int j = 0; j < jmax; j++)
			grid[i][j] = (rand() % 2);

	//clock start
	float startTime = omp_get_wtime();

	//using dividing method to decrease computational cost 
	for (int n = 0; n < max_steps/4; n++)
	{
		cout << "it: " << n << endl;
		do_split_method(4);
		grid_to_file(n);
	}	

	//method without split up
	//for (int n = 0; n < max_steps; n++)
	//{
	//	cout << "it: " << n << endl;
	//	do_iteration();
	//	grid_to_file(n);
	//}	

	//clock end, print out runnning time
	float endTime = omp_get_wtime();
	printf("running time =: %f\n", endTime - startTime);
	return 0;
}
