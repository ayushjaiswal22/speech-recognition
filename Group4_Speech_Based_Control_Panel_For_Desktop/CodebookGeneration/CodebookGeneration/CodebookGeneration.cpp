// CodebookGeneration.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include<iostream>
#include<fstream>
#include<cmath>

#define CBS 32
using namespace std;

double bin[CBS][20000][12];
double data[30000][12];      //array that will be having data of the universe
double TokhureDistance(double a[], double b[])
{
	double weight[12] = { 1, 3, 7, 13, 19, 22, 25, 33, 42, 50, 56, 61 };   //Tokhure's Weight
	double sum = 0;
	//Calculation of Tokhure's Distance
	for (int i = 0; i < 12; ++i)
	{
		sum += weight[i] * (a[i] - b[i])*(a[i] - b[i]);
	}
	return sum;
}

void centroid(double a[], double b[][12], int size)
{
	//a is the pointer to the array in which centroid vector of dimension 12 is to be stored
	//b[][] is the array of rows equal to size(vectors of dimension 12) whose centroid is to be calculated


	//Initialization of the array where data is to be stored

	for (int i = 0; i < 12; ++i)
		a[i] = 0;

	//Calculations of the sum of the columns

	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < 12; ++j)
		{
			a[j] += b[i][j];
		}
	}

	//Calculation of the average (means the centroid)
	for (int i = 0; i < 12; ++i)
		a[i] = a[i] / size;
}


void copy(double a[], double b[])
{
	//copying vector b[] to vector a[]
	for (int i = 0; i < 12; ++i)
		a[i] = b[i];
}

double distortion(double a[],double b[][12], int size)
{
	//calculate distortion using Tokhure's distance
	double sum = 0;
	for (int i = 0; i < size; ++i)
	{
		sum += TokhureDistance(a, b[i]);
	}
	return sum;
}

void Kmeans(double cb[][12],int cb_size,double data[][12],int datasize)
{
	//This function uses k means algo to find best set of centroid for data[][] w.r.t. cb[] vectors
	//change
	   //bin can store upto 1000 vectors of dimension 12 in each bin (bins atmost 8)
	int bin_size[CBS];  //contains size of the  bins



	//var cb_size = vectors in code book presently
	//var datasize =  total vectors to be classifed

	double old_distortion = 9999999999;   //Initial arbitrary value
	double new_distortion;
	double difference=99999;    //contains the value of difference of old and new distortion

	while (difference>0.00000001)
	{
		difference = 0;
		//Initialization of the bin_size
		for (int i = 0; i < CBS; ++i)
			bin_size[i] = 0;
		new_distortion = 0;

		double min_distance;    //minimum distance from the vectors of the codebook
		int index_bin;               //index of the bin where vector is to be placed

		//classification into the bins
		for (int i = 0; i < datasize; ++i)
		{
			min_distance = 999999999;  //Initialization of min_distance
			index_bin = 0;   //Initializaton of the index_bin

			for (int j = 0; j < cb_size; ++j)
			{
				double distance = TokhureDistance(cb[j], data[i]);
				if (distance < min_distance)
				{
					min_distance = distance;
					index_bin = j;
				}
			}

			//copy data to bin
			copy(bin[index_bin][bin_size[index_bin]], data[i]);   //copying
			bin_size[index_bin] += 1;    //size increased

		}

		//Calculation of the new centroid
		for (int i = 0; i < cb_size; ++i)
		{
			centroid(cb[i], bin[i], bin_size[i]);
		}
		
		//calculation of the distortion
		for (int i = 0; i < cb_size; ++i)
		{
			new_distortion += distortion(cb[i], bin[i],bin_size[i]);
		}
		cout << "\nBins and the elements in those bins are:" << endl;
		for (int i = 0; i < cb_size; ++i)
		{
			cout << "bin" << i + 1 << ": " << bin_size[i] <<" elements"<< endl;
		}
		new_distortion /= cb_size;
		cout.precision(40);
		cout << "Distortion :" << new_distortion << endl;
		difference = old_distortion - new_distortion;
		old_distortion = new_distortion;
		cout <<"Difference(Old Distortion  - New distortion) :"<< difference << endl;
	}

}


int _tmain(int argc, _TCHAR* argv[])
{
	
	int size_data = 0;    //contains the size of the array data[]
	ifstream in;
	in.open("..\\HMM_universe.txt", ios::in);

	//Saving data into the array
	while (!in.eof())
	{
		for (int i = 0; i < 12; ++i)
		{
			in >> data[size_data][i];
		}
		size_data++;
	}
	in.close();

	cout << "Total vectors in the universe are:" << size_data << endl;
	cout << endl;

	double CodeBook[CBS][12];   //Our Required Codebook                
	double temp_cb[CBS][12];   //In temp_cb splitting will be there (temporarily) at last data will be copied back to CodeBook[][]
	int codebook_size = 0;  // Keeps track of the current size of the code book
	int tempcb_size = 0;           //Size of the temp_cb[] array


	int choice;
	cout << "Enter the choice number:\n1. Make code Book by K means\n2. Make code book by LBG" << endl;
	cin >> choice;

	switch (choice)
	{
	case 1:
	{
		//Initialization of the codebook for k means with a vector of each vowel utterence randomly
		for (int i = 0; i < CBS; ++i)
		{
			if (i < 5)     //five from different vowels
			{
				for (int j = 0; j < 12; ++j)
				{
					CodeBook[i][j] = data[50 * i][j];
				}
			}
			else
			{
				for (int j = 0; j < 12; ++j)    //sarting element with respect to index
				{
					CodeBook[i][j] = data[i][j];
				}
			}

		}
		//Applyting K means
		Kmeans(CodeBook, CBS, data, size_data);

		//Printing the code book made by K means
		ofstream out;
		out.open("..\\codebook_kmeans.txt", ios::out);
		cout << "Codebook is as follows:" << endl;
		for (int i = 0; i < CBS; ++i)
		{
			cout << "Row:" << i + 1 << endl;
			for (int j = 0; j < 12; ++j)
			{
				cout << CodeBook[i][j] << endl;
				out << CodeBook[i][j] << "\t";
			}
			cout << endl;
			out << endl;
		}
		out.close();

		break;
	}
	case 2:
	{
		//Calculation of centroid and iniialization of the codebook
		centroid(CodeBook[0], data, size_data);
		codebook_size++;
		cout << "Initial codebook vector is as follows:" << endl;
		for (int i = 0; i < 12; ++i)
			cout << CodeBook[0][i] << endl;
		cout << endl;


		//Now sum array contains the centroid of the from where we will start LBG

		//LBG  Algorithm 
		while (codebook_size < CBS)
		{
			tempcb_size = 0;   //initialization of tempcb_size

			//Splitting the code book with splitting parameter 0.03
			for (int i = 0; i < codebook_size; ++i)
			{
				for (int j = 0; j < 12; ++j)
				{
					temp_cb[tempcb_size][j] = CodeBook[i][j] + 0.03;
					temp_cb[tempcb_size + 1][j] = CodeBook[i][j] - 0.03;
				}
				tempcb_size += 2;
			}
			cout << "Temp code book size after splitting:" << tempcb_size << endl;

			//Applying K means algorithm
			Kmeans(temp_cb, tempcb_size, data, size_data);


			//copying temp_cb to actual codebook
			for (int i = 0; i < tempcb_size; ++i)
			{
				copy(CodeBook[i], temp_cb[i]);
			}
			//updating of original codebook size
			codebook_size = tempcb_size;

		}


		//Printing of the codebook and sending it to the file
		ofstream out;
		out.open("..\\codebook_lbg.txt", ios::out);
		cout << "Codebook is as follows:" << endl;
		for (int i = 0; i < CBS; ++i)
		{
			cout << "Row:" << i + 1 << endl;
			for (int j = 0; j < 12; ++j)
			{
				cout << CodeBook[i][j] << endl;
				out << CodeBook[i][j] << "\t";
			}
			cout << endl;
			out << endl;
		}
		out.close();
		break;
	}
	default:
		cout << "Please enter appropriate choice" << endl;
		break;
	}


	return 0;
}

