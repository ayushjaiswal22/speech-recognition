// 184101010_Testing_HMM_Model.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#define T_MAX 300
#define MAX_SAMPLES 40000
#define M 32	//Number of observational symbols
#define N 5		//Number of unknown states

long double Prob[10] = {0};
int ind=0;
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
long double Tokhura_distance[32] = {0};
int Observation_Sequence[T_MAX] = {0};
long double codebook[32][12] = {0};
long double work_arr2[T_MAX][12] = {0};
long double alpha[T_MAX][N+1] = {0};		//forward variable
long double A[N+1][N+1] =	{0};	//Probabilities of transitions
long double B[N+1][M+1] =	{0};	//Probabilities of obsevational symbols
long double pi[N+1] =	{0};		//Probabbilties of states being the initial state

int T = 0;
bool flag = 0;
bool flag1 = 0;
void initialize_work_arr()
{
	int i=0,j=0;
	for(i=0;i<T_MAX;i++)
	{
		for(j=0;j<12;j++)
		{
			work_arr2[i][j] = 0;
		}
	}
}

void fileread()
{
	int i=0,j=0;
	ifstream infile;
	string location;
	location = "Cepstral_coefficients_file\\";
	location = location + "184101010_cep.txt";
	infile.open(location.c_str());
	infile >> T;
	while(!infile.eof())
	{
		for(i=0;i<T;i++)
		{
			for(j=0;j<12;j++)
			{
				infile >> work_arr2[i][j];
				//cout << work_arr2[i][j] << " ";
			}
			//cout<<endl;
		}
	}

	infile.close();
}

int  min_Tokhura()
{
	int i = 0, index = 0;
	long double min = 999999;
	for(i =0 ; i < 32 ; i++)
	{
		if(Tokhura_distance[i] < min)
		{
			index = i;
			min = Tokhura_distance[i];
		}
	}
	return index+1;

}


void calculate_tokura_distance()
{
	ofstream outfile;
	string location;
	location = "Cepstral_coefficients_file\\";
	location = location + "184101010_Obs_seq.txt";
	outfile.open(location.c_str());
	int i=0,j=0,k=0;
	long double sum =0;
	outfile << T << endl;
	for(k=0;k<T;k++)
	{
		for(i=0;i<32;i++)
		{
			sum =0;
			for(j=0;j<12;j++)
			{
				sum = sum + (Tokhura_Weights[j]*((codebook[i][j]-work_arr2[k][j])*(codebook[i][j]-work_arr2[k][j])));
			}
			Tokhura_distance[i] = sum;
		}
		Observation_Sequence[k] = min_Tokhura();
		outfile << Observation_Sequence[k] << endl;
	}
	outfile.close();
}


void Compute_DC_Shift(int f)
{
	//Computation of DC_Shift
   ifstream infile1;
   ifstream infile2;
   string location;
   if(f==1)
   {
	   cout<<"Recording Silence!!!"<<endl<<endl;
	   system("Recording_Module.exe 5 Silence.wav Silence.txt");
   }
   cout<<"Recording Signal"<<endl<<endl;
   system("Recording_Module.exe 4 Sound.wav 184101010_digit_full.txt");
   ifstream infile3;
   ofstream outfile;
   outfile.open("184101010_digit.txt");
   infile3.open("184101010_digit_full.txt");

   for(int i=0;i<20000;i++)
	   infile3 >> sample;

   for(int i=0;i<30000;i++)
   {
	   infile3 >> sample;
	   outfile << sample << endl;
   }

   location = "184101010_digit.txt";
   infile1.open("Silence.txt");
   infile2.open(location.c_str());

   ofstream outfile1;
   outfile1.open("DC_Shifted.txt");

  // cout << "Reading from the file" << endl;
   while (!infile1.eof())							
    {
        infile1 >> sample;
        sum = sum + sample;
        total_num_samples++;
    }
    mean = (sum/total_num_samples) ;														//Calculating the mean
  //  cout << "DC Shift value is: " << mean << endl;

	int flag = 0;

	total_num_samples = 0;
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

	//cout << "Total number of samples: " << total_num_samples << endl;
	//T = ((total_num_samples-320)/80)+1;
	//cout<< "T: "<<T<<endl;
    //cout << "maximum: " << maximum << endl;
	//cout << "Index: " << index << endl; 

   // close the opened file.

   infile1.close();
   infile2.close();
   outfile1.close();
}

void Normalization()
{
	long double d[40000] = {0};
	long double data[40000] = {0};
	int data_size=0;
	int ds = 0;
	ifstream infile3;
   ofstream outfile2;
   infile3.open("DC_Shifted.txt");
   outfile2.open("Normalised.txt");
   while(!infile3.eof())
   {
        infile3 >> sample ;
		ds++;
        sample = 5000/maximum*sample;										//Normalization Process
        outfile2 << sample << endl;
	}
	//cout << "Normalized" << endl;

	infile3.close();
	outfile2.close();

	ofstream oufile2;
	
	infile3.open("Normalised.txt");
	int l=0;
	while(!infile3.eof())
	{
		infile3 >> d[l];
		l++;
	}
	ds = l;
	//cout<<"Done storing the data...."<<endl;
	long double frame_energy[400];
	int frame_zcr[400];

	int start_point = 0;
	int end_point = 0;
	ofstream out;
	ofstream out1;
	out1.open("energy.txt", ios::out);
	out.open("zcr.txt", ios::out);

	int frm_counter = 0;     //stores total number of frames

	for (int i = 0; i < ds; i = i + 80)
	{
		//0 is for negative and 1 is for positive
		int sign1=1;   //+ve by default
		int sign2 = 1;   //+ve by default
		//initializeing sign1
		if (d[i] < 0)
			sign1 = 0;

		int zcr = 0;
		double energy=0;

		for (int j = i; j < 320 + i && j<ds; ++j)
		{
			//code for zcr
			if (d[j] < 0)
				sign2 = 0;    //-ve
			else if (d[j]>0)
				sign2 = 1;    //+ve
			if (sign2 != sign1)
				zcr++;
			//sign updation
			sign1 = sign2;
			energy += (d[j] * d[j]);

		}

		//storing frmae wise analysis to array
		frame_energy[frm_counter] = energy;
		frame_zcr[frm_counter] = zcr;
		frm_counter++;
		out << zcr << endl;
		out1 << energy << endl;
	}

	out.close();
	out1.close();
		
	//finding starting point

	for (int i = 0; i < frm_counter; ++i)
	{
		int flag1 = 0;
		int flag2 = 0;
		if (frame_energy[i] > 9.0e+7 /* || frame_zcr[i] > 50 */)
		{
			for (int j = 0; j < 8; ++j)
			{
				if (frame_energy[i + j] <= 9.0e+7 )
				{
					flag1 = 1;
					break;
				}
			}

			/*for (int j = 0; j < 8; ++j)
			{
				if (frame_zcr[i + j] <= 50)
				{
					flag2 = 1;
					break;
				}
			}*/
			
			if (flag1 == 0 /* || flag2==0*/)
			{
				start_point = i * 80;
				break;
			}
		}
	}


	//ending point
	for (int i = frm_counter-1; i >= 0; --i)
	{
		int flag1 = 0;
		int flag2 = 0;
		if (frame_energy[i] > 9.0e+7/* || frame_zcr[i]>50*/)
		{
			for (int j = 0; j > (-8); --j)
			{
				if (frame_energy[i + j] <= 9.0e+7)
				{
					flag1 = 1;
					break;
				}
			}
			/*	for (int j = 0; j > (-8); --j)
			{
				if (frame_zcr[i + j] <= 50)
				{
					flag2 = 1;
					break;
				}
			}*/
				if (flag1 == 0 /*|| flag2 == 0*/)
			{
				end_point = i * 80;
				break;
			}
		}
	}
//	cout << "Starting point :" << start_point << endl;
	//cout << "End point :" << end_point << endl;
		
	for (int i = start_point; i<=end_point; ++i)
	{
		data[data_size] = d[i];
		data_size++;
	}
		//applying the normalization
	
	out.open("trimmed.txt", ios::out);
	for (int j = 0; j < data_size; ++j)
	{
		out << data[j] << endl;
	}
	out.close();
	ifstream infile4;
	infile4.open("trimmed.txt");
	int i = 0;

	while(!infile4.eof())
	{
		infile4 >> sample;
		arr[i++] = sample;
	}
	
	infile4.close();
	T = ((i-320)/80)+1;
//	cout<< "T: "<<T<<endl;
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
	location = "Cepstral_coefficients_file\\";
	location = location + "184101010_cep.txt";
	outfile7.open(location.c_str(),ios::app);
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
	}
	index2++;

	outfile7.close();

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
	/*int i=0,j=0;
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
	}*/
	int i=0,j=0;
	for(i=0;i<9;i++)
		Prob[i] = 0;
	ind = 0;
	total_num_samples = 0;
	for(i=0;i<MAX_SAMPLES;i++)
		arr[i] = 0;
	sample=0;
	sum = 0;
	mean = 0;
	maximum = 0;
	index = 0;
	index2 = 0;
	for(i=0;i<320;i++)
	{
		work_arr[i] = 0;
		Ham_Win[i] = 0;
	}
	for(i=0;i<13;i++)
		R[i] = 0;

	for(i=0;i<5;i++)
		for(j=0;j<12;j++)
			test_data[i][j] = 0;
	for(i=0;i<32;i++)
		Tokhura_distance[i] = 0;
	for(i=0;i<T_MAX;i++)
		Observation_Sequence[i] = 0;
	for(i=0;i<32;i++)
		for(j=0;j<12;j++)
			codebook[i][j] = 0;
	for(i=0;i<T_MAX;i++)
		for(j=0;j<12;j++)
			work_arr2[i][j] = 0;
	for(i=0;i<T_MAX;i++)
		for(j=0;j<N+1;j++)
			alpha[i][j] = 0;		//forward variable

	for(i=0;i<N+1;i++)
		for(j=0;j<N+1;j++)
			A[i][j] = 0;	//Probabilities of transitions
	for(i=0;i<N+1;i++)
		for(j=0;j<M+1;j++)
			B[i][j] = 0;	//Probabilities of obsevational symbols
	for(i=0;i<N+1;i++)
	pi[i] = 0;		//Probabbilties of states being the initial state

	T = 0;
	flag = 0;
	flag1=0;
}
void initialize_code_book_vector()
{
	ifstream infile;
	infile.open("codebook.txt");
	int i=0,j=0;

	for(i=0;i<32;i++)
	{
		for(j=0;j<12;j++)
		{
			infile >> codebook[i][j];
		}
	}
	infile.close();
}

void initialize_A_B_Pi_Matrix(long double i1)
{
	ifstream infile1;
	ifstream infile2;
	ifstream infile3;
	string location1,location2,location3;
	location1 = "Average_Model\\184101010_";
	location2 = "Average_Model\\184101010_";
	location3 = "Average_Model\\184101010_PI.txt";
	
	location1 = location1 + to_string(i1) + "A.txt";
	location2 = location2 + to_string(i1) + "B.txt";

	int i=0,j=0;//variables used for loop
	
	//Initializing the A matrix
	infile1.open(location1);
	for(i=0;i<=N;i++)
	{
		for(j=0;j<=N;j++)
		{
			infile1>>A[i][j];
			//cout<<A[i][j]<<"  ";
		}
		//cout<<endl;
	}
	//Initializing the B matrix
	infile2.open(location2);

	for(i=0;i<=N;i++)
	{
		for(j=0;j<=M;j++)
		{
			infile2>>B[i][j];
			//cout<<B[i][j]<<"  ";
		}
		//cout<<endl;
	}
	//Initializing the PI matrix
	infile3.open(location3);

	for(i=0;i<=N;i++)
	{
		infile3>>pi[i];
	}

	infile1.close();
	infile2.close();
	infile3.close();
}

void forward_procedure()
{
	int i=0,j=0,t=0;		//varaibles used for loop

	//Initialization Step
	for(i=1;i<=N;i++)
	{
		alpha[1][i] = pi[i]*B[i][Observation_Sequence[1]];
	}

	//Induction Step
	long double sum = 0;		//variable to store local SUM
	for(t=1;t<T;t++)
	{
		for(j=1;j<=N;j++)
		{
			sum=0;
		
			for(i=1;i<=N;i++)
			{
				sum+=alpha[t][i]*A[i][j];
				//cout<<sum<<" ";
			}
			//cout<<endl;
			alpha[t+1][j] = sum*B[j][Observation_Sequence[t+1]];
		}
	}
	sum=0;
	//Termination Step
	for(i=1;i<=N;i++)
	{
		sum += alpha[T][i];
	}

	/*for(i=1;i<=T;i++)
	{
		for(j=1;j<=N;j++)
		{
			cout<<alpha[i][j]<<"\t";
		}
		cout<<endl;
	}*/
	
//	cout<<"the probability is "<<sum<<endl;
	Prob[ind++] = sum;

}
void read_observation_sequence()
{
	int i=0;
	string location;
	location = "Cepstral_coefficients_file\\";
	location = location + "184101010_Obs_seq.txt";
	ifstream infile;
	infile.open(location);
	infile >> T;
	for(i=1;i<=T;i++)
	{
		infile >> Observation_Sequence[i];
	}
	infile.close();
}
int Test (int f) 
{
		string names[7] = {"All","Browser","Close","Home","Music","Settings","Vlc"};
		int status = remove("Cepstral_coefficients_file\\184101010_cep.txt");
		//cout<<"file removal status: "<<status<<endl;
		int initial = 0;
		int i=0,k=0,index=0;
		initialize_global_variables();
		initial = 0;
		i=0;
		k=0;
		index=0;
		Compute_DC_Shift(f);
		Normalization();
		if(T > 400)
		{
			cout<< "Cannot analyze the data..." <<endl;
			return 0;
		}
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
		i=0;
		initialize_code_book_vector();
		T=0;
		initialize_work_arr();
		fileread();
		calculate_tokura_distance();
		read_observation_sequence();

		for(i=0;i<=6;i++)
		{
			initialize_A_B_Pi_Matrix(i);
			forward_procedure();
		}
		long double max = 0;
		int max_ind = 0;
		for(i=0;i<=6;i++)
		{
			if(max < Prob[i])
			{
				max=Prob[i];
				max_ind = i;
			}
		}
		cout<< "WORD Recognized is: "<<names[max_ind]<<endl;
		return max_ind;
}
void main(){
		int ind1,ind2,flg=1;int input=0;
		while(!input){
			
		ind1 = Test(flg);
		flg=0;

		switch (ind1)
		{
			case 0:
				cout<<"Not in the dictionary"<<endl;
				break;
			case 1:
				cout<<"Opening Browser"<<endl;
				system("start Chrome.exe");
				break; 
			case 2:
				cout<<"Close What"<<endl;
				flag1=1;
				ind2 = Test(flg);
				switch (ind2){
				case 0:
					cout<<"Closing All Apps"<<endl;
					break;
				case 1:
					cout<<"Closing Browser"<<endl;
					system("taskkill /im chrome.exe");
					break;
				case 3:
					cout<<"Closing Home"<<endl;
					break;
				case 4:
					cout<<"Closing Music"<<endl;
					system("taskkill /im wmplayer.exe");
					break;
				case 5:
					cout<<"Cannot Close Settings"<<endl;
					
					break;
				case 6:
					cout<<"Closing VLC"<<endl;
					system("taskkill /im vlc.exe");
					break;
				}
				break;
			case 3:
				cout<<"Opening Home"<<endl;
				system("explorer =");
					
				break;
			case 4:
				cout<<"Opening Music"<<endl;
				system("start wmplayer.exe");
				break;
			case 5:
				cout<<"Opening Settings"<<endl;
				system("start control panel");
				break;
			case 6:
				cout<<"Opening VLC"<<endl;
				system("start vlc.exe");
				break;
		}
		
		cout<<"Press 0 to continue...."<<endl;
		cin>>input;
	}
}