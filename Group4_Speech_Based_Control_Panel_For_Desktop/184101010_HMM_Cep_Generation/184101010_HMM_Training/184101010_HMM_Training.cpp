// 184101010_HMM_Training.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#define MAX_SAMPLES 20000

int total_num_samples = 0;
long double arr[MAX_SAMPLES]= {0};
long double sample;
double long sum = 0;
long double mean = 0;
long double maximum = 0;
int index = 0;
int index2 = 0;
long double Ham_Win[320];
long double R[13];
long double work_arr[320] = {0};
long double test_data[5][12] = {0};
long double Tokhura_Weights[12] = {1,3,7,13,19,22,25,33,42,50,56,61};
long double lp1,lp2;
int T = 0;
bool flag = 0;

string names[7] = {"all","browser","close","home","music","settings","vlc"};

void Compute_DC_Shift()
{
	//Computation of DC_Shift
   ifstream infile1;
   ifstream infile2;
   string location;
   location = "WordData\\";
   string temp;
	if(lp1 == 0)
		temp = "all";
	else if(lp1==1)
		temp = "browser";
	else if(lp1 ==2)
		temp = "close";
	else if(lp1 == 3)
		temp = "home";
	else if (lp1 == 4)
		temp = "music";
	else if (lp1 == 5)
		temp = "settings";
	else if (lp1 == 6)
		temp = "vlc";
   cout<<temp<<endl;
   location = location + temp + "\\" + temp + "_" +to_string(lp2) + ".txt";
   infile1.open("silence.txt");
   cout<<"after silence"<<endl;
   infile2.open(location.c_str());
   ofstream outfile1;
   outfile1.open("DC_Shifted.txt");

   cout << "Reading from the file" << endl;

   int a=0;

   while (!infile1.eof())							
   {
			infile1 >> sample;
			sum = sum + sample;
			total_num_samples++;
    }
    mean = (sum/total_num_samples) ;														//Calculating the mean
    cout << "DC Shift value is: " << mean << endl;

	int flag = 0;
	string line;
	total_num_samples = 0;
   	/*for(int i =0 ; i <5 ; i++)
	{
		infile2.ignore(256,'\n');
	}*/
	while (!infile2.eof())
    {
	
        infile2 >> sample;
		total_num_samples++;
        sample = sample - mean;									//Subtracting the DC_Shift value and storing it in seperate file
        outfile1 << sample << endl;
        if(flag == 0)
        {
            flag =1;
            maximum = sample;
			index = total_num_samples;
        }
        if(abs(sample) > maximum)									//Calculating the maximum amplitude for Normalization
        {
            maximum = abs(sample);
			index = total_num_samples;
        }

    }

	cout << "Total number of samples: " << total_num_samples << endl;
	T = (((total_num_samples-320)/80)+1);
	cout<< "T: "<<T<<endl;
    cout << "maximum: " << maximum << endl;
	cout << "Index: " << index << endl; 

   // close the opened file.
   infile1.close();
   infile2.close();
   outfile1.close();
}

void Normalization()
{
	ifstream infile3;
   ofstream outfile2;
   infile3.open("DC_Shifted.txt");
   outfile2.open("Normalised.txt");
   while(!infile3.eof())
   {
        infile3 >> sample ;
        sample = 5000/maximum*sample;										//Normalization Process
        outfile2 << sample << endl;
	}
	cout << "Normalized" << endl;

	infile3.close();
	outfile2.close();

	ifstream infile4;
	infile4.open("Normalised.txt");
	int i = 0;

	while(!infile4.eof())
	{
		infile4 >> sample;
		arr[i++] = sample;
	}
	
	infile4.close();
}

void Compute_Hamming_Window()
{	// Applying Hamming Window
	int i = 0;

	for(i=0;i<320;i++)
	{
		Ham_Win[i] = (0.54 - 0.46*cos((2*3.14*i)/319));
	}


	// End of Hamming Window Function

}

void Calculate_Ci_Prime(long double Ci[])
{
	//Calculating Ci Prime Values

	//Procedure for Liftering
	string location,location2;
	int i =0,j=0,P=12;
	ofstream outfile7;
	ofstream outfile8;
	location = "Cepstral_coefficients_file\\digit";
	location = location + to_string(lp1) +"\\" + "184101010_" + to_string(lp1) + "_" + to_string(lp2) + "cep.txt";
	outfile7.open(location.c_str(),std::ios::app);
	location2 = "Cepstral_coefficients_file\\HMM_universe.txt";
	outfile8.open(location2.c_str(),std::ios::app);
	if(flag == 0)
	{
		flag = 1;
		outfile7 << T << endl;
	}
	long double raised_window[12];
	long double ci_primes[13] = {0};

	for(i=1;i<P+1;i++)
	raised_window[i-1] = P/2 * sin((3.14*i)/P) + 1;

	for(j=1;j<P+1;j++)
	{
		ci_primes[j] = raised_window[j-1] * Ci[j];
	}

	for(j=1;j<P+1;j++)
	{
		test_data[index2][j-1] = ci_primes[j];
		outfile7 << ci_primes[j] <<endl;
		outfile8 << ci_primes[j] <<endl;
	}
	index2++;

	outfile7.close();
	outfile8.close();

	//End of Ci Prime Calculation
}

void Calculate_Ci(long double ai[])
{
	//Start of Calculation of Ci Values

	long double C[13] = {0};
	int i =0, j=0;
	int P =12;
	ofstream outfile6;
	outfile6.open("Ci.txt");
	C[0] = logl(R[0]);					//R[0] is the Energy Value.

	for(i=1;i<=P;i++)
	{
		sum = 0;
		for(j=1;j<=(i-1);j++)
		{ 
			sum = sum + (j*C[j]*ai[i-j]/i);			//Computing the Ci Values
		}
		C[i] = ai[i] + sum;
	}

	for(i=1;i<=P;i++)
		outfile6 << C[i] << endl;

	outfile6.close();

	//End of calculation if Ci Values
	Calculate_Ci_Prime(C);
}

void Calculate_Ai(long double R[])
{
		//Start of Ai Calculation
	int i = 0,j=0;
	int no_of_Frames = 5;
	int P = 12;						//No of coefficients
	int loop =0;
	long double E[13] = {0},k[13] = {0};
	long double alpha[13][13] = {0};
	long double alpha_Final[13] = {0};

	ofstream outfile5;
	outfile5.open("Ai.txt");
	
		E[0] = R[0];
		for(i=1;i<=P;i++)
		{
			sum = 0;
			for(j=1;j<=i-1;j++)
			{

				sum = sum + (alpha[i-1][j]*R[i-j]);
			}
			k[i] = (R[i] - sum)/E[i-1];
				
			alpha[i][i] = k[i];

			for(j=1;j<=i-1;j++)
			{
				alpha[i][j] = alpha[i-1][j] - k[i]*alpha[i-1][i-j];
			}

			E[i] = (1-(k[i]*k[i]))*E[i-1];
		}

		for(j=1;j<=P;j++)
		{
			outfile5 << alpha[P][j] << " ";
									//Printing the Ai Values
			alpha_Final[j] = alpha[P][j];
		}

		outfile5 << endl; 
	outfile5.close();
	
	//End of Ai calculation function

	Calculate_Ci(alpha_Final);
}

void Calculate_Ri(long double work_arr[])
{
		//Now calculating the frame wise Ri's

	//Initializing and declaring the variables being used for the computation

	int i = 0,j=0,k=0;
	int offset = 0;			//For Frame differentiation
	ofstream outfile4;
	outfile4.open("Ri.txt",std::ios::app);
	
	sum = 0;

	for(i = 0;i < 13;i++)
	{
		k=0;j=i;
		while(j < 320)
		{
			sum = sum + (work_arr[k]*work_arr[j]);
			k++;
			j++;
		}
		R[i] = sum;
		outfile4 << R[i] << endl;
		sum = 0;

	}
	outfile4.close();
	//End of Ri Calculation
	Calculate_Ai(R);

}

void Hamming(long double work_arr[])
{
    for(int i=0;i<320;i++)
    {
        work_arr[i]=Ham_Win[i]*work_arr[i];
    }
}

void initialize_global_variables()
{
	int i=0,j=0;
		for(i=0;i<320;i++)
		{
			work_arr[i] = 0;
			Ham_Win[i] = 0;
		}
	total_num_samples = 0;
	
	for(i=0;i<MAX_SAMPLES;i++)
	arr[MAX_SAMPLES]= 0;
	
	flag = 0;
	sample=0;
	sum = 0;
	mean = 0;
	maximum = 0;
	index = 0;
	index2 = 0;

	for(i=0;i<13;i++)
	{
		R[i] = 0;
	}

	for(i=0;i<5;i++)
	{
		for(j=0;j<12;j++)
		{
			test_data[i][j] = 0;
		}
	}
}

int main () 
{
	int initial = 0;
	int i=0,k=0,index=0;
	int loop1,loop2;
	for(loop1=0;loop1<=6;loop1++)
	{
		lp1=loop1;
		for(loop2=1;loop2<=20;loop2++)
		{
			cout<<"here,,,,"<<loop1<<" "<<loop2<<endl;
			initialize_global_variables();
			initial = 0;
			i=0;
			k=0;
			index=0;
			lp2=loop2;
			Compute_DC_Shift();
			Normalization();
			Compute_Hamming_Window();
			for(k=0;k<T;k++)
			{
			    for(i=0;i<320;i++)
				{
					work_arr[i]=arr[initial+i];
				}
				Hamming(work_arr);
				Calculate_Ri(work_arr);
				initial+=80;
			}
			T=0;

		}
	}
}


