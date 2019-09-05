// 184101010_training_HMM_model.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
int T;
#define M 32	//Number of observational symbols
#define N 5		//Number of unknown states
#define T_MAX 300	//Max time of observations 

int Observation_Sequence[T_MAX] = {0};

long double new_pstar = 0.00005;
long double old_pstar = 1.0e-200;

long double work_arr[T_MAX][12] = {0};
long double codebook[32][12] = {0};
long double Tokhura_Weights[12] = {1,3,7,13,19,22,25,33,42,50,56,61};
long double Tokhura_distance[32] = {0};
long double A[N+1][N+1] =	{0};	//Probabilities of transitions
long double B[N+1][M+1] =	{0};	//Probabilities of obsevational symbols
long double pi[N+1] =	{0};		//Probabbilties of states being the initial state
int obser[T_MAX] =	{0};			//Observation sequence
long double alpha[T_MAX][N+1] = {0};		//forward variable
long double beta[T_MAX][N+1] = {0};		//backward variable
long double gamma[T_MAX][N+1];			//gamma variable
long double xi[T_MAX][N+1][N+1];			//xi variable

//for calculating the average model
long double A_Avg[N+1][N+1] =	{0};	//Probabilities of transitions
long double B_Avg[N+1][M+1] =	{0};	//Probabilities of obsevational symbols

void fileread(long double i1,long double i2)
{
	int i=0,j=0;
	ifstream infile;
	string location;
	location = "..\\HMM_data\\Cepstral_coefficients_file\\digit";
	location = location + to_string(i1) + "\\" + "184101010_" + to_string(i1) + "_" + to_string(i2) + "cep.txt";
	infile.open(location.c_str());
	infile >> T;
	while(!infile.eof())
	{
		for(i=0;i<T;i++)
		{
			for(j=0;j<12;j++)
			{
				infile >> work_arr[i][j];
				//cout << work_arr[i][j] << " ";
			}
			//cout<<endl;
		}
	}

	infile.close();
}

void initialize_work_arr()
{
	int i=0,j=0;
	for(i=0;i<T_MAX;i++)
	{
		for(j=0;j<12;j++)
		{
			work_arr[i][j] = 0;
		}
	}
}

void initialize_code_book_vector()
{
	ifstream infile;
	infile.open("..\\HMM_data\\codebook.txt");
	int i=0,j=0;

	for(i=0;i<32;i++)
	{
		for(j=0;j<12;j++)
		{
			infile >> codebook[i][j];
			//cout << codebook[i][j] << " ";
		}
		//cout<<endl;
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

void calculate_tokura_distance(long double i1,long double i2)
{
	ofstream outfile;
	string location;
	location = "..\\HMM_data\\Observation_Sequence_file\\digit";
	location = location + to_string(i1) + "\\" + "184101010_" + to_string(i1) + "_" + to_string(i2) + "_Obs_seq.txt";
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
				sum = sum + (Tokhura_Weights[j]*((codebook[i][j]-work_arr[k][j])*(codebook[i][j]-work_arr[k][j])));
			}
			Tokhura_distance[i] = sum;
		}
		Observation_Sequence[k] = min_Tokhura();
		outfile << Observation_Sequence[k] << endl;
	}
	outfile.close();
}

void storing_A_B_Pi_Matrix(long double i1,long double i2)
{

	string A_location,B_location,Pi_location;			//To store the locations of the respective files
	ofstream outfile1,outfile2,outfile3;				//File pointers to all the matrix files
	A_location = "..\\HMM_data\\A_Matrix_file\\digit";
	B_location = "..\\HMM_data\\B_Matrix_file\\digit";
	Pi_location = "..\\HMM_data\\Pi_Matrix_file\\digit";
	A_location = A_location + to_string(i1) + "//" + "184101010_" + to_string(i1) + "_" + to_string(i2) + "_A.txt";
	B_location = B_location + to_string(i1) + "//" + "184101010_" + to_string(i1) + "_" + to_string(i2) + "_B.txt";
	Pi_location = Pi_location + to_string(i1) + "//" + "184101010_" + to_string(i1) + "_" + to_string(i2) + "_Pi.txt";
	outfile1.open(A_location);
	outfile2.open(B_location);
	outfile3.open(Pi_location);

	int i =0,j=0;				//looping variables
	
	for(i=0;i<=N;i++)
	{
		for(j=0;j<=N;j++)
		{
			outfile1<<A[i][j]<<"\t";
		}
		outfile1<<endl;
	}

	for(i=0;i<=N;i++)
	{
		for(j=0;j<=M;j++)
		{
			outfile2<<B[i][j] << "\t";
		}
		outfile2<<endl;
	}

	for(i=0;i<=N;i++)
	{
		outfile3<<pi[i] << "\t";
	}

	outfile1.close();
	outfile2.close();
	outfile3.close();
}



void backward_procedure()
{
	int i=0,j=0,t=0;		//varaibles used for loop

	//Initialization Step
	for(i=1;i<=N;i++)
		beta[T][i] = 1;


	//Induction Step
	long double sum = 0;		//variable to store local SUM
	for(t=T-1;t>=1;t--)
	{
		for(i=1;i<=N;i++)
		{
			sum=0;

			for(j=1;j<=N;j++)
			{
				sum+=A[i][j]*B[j][obser[t+1]]*beta[t+1][j];
			}
			beta[t][i] = sum;
		}
	}
	//cout<<"BETA values/n/n"<<endl;
	/*for(t=1;t<=T;t++)
	{
		for(j=1;j<=N;j++)
		{
			cout<<beta[t][j]<<"\t";
		}
	cout<<endl<<endl;
	}*/
}


void initialize_A_B_Pi_Matrix()
{
	ifstream infile1;
	ifstream infile2;
	ifstream infile3;
	ifstream infile4;

	
	int i=0,j=0;//variables used for loop
	
	//Initializing the A matrix
	infile1.open("..\\HMM_data\\Initial_Model\\A.txt");
	for(i=0;i<=N;i++)
		for(j=0;j<=N;j++)
			infile1>>A[i][j];

	//Initializing the B matrix
	infile2.open("..\\HMM_data\\Initial_Model\\B.txt");

	for(i=0;i<=N;i++)
		for(j=0;j<=M;j++)
			infile2>>B[i][j];
	
	//Initializing the PI matrix
	infile3.open("..\\HMM_data\\Initial_Model\\PI.txt");

	for(i=0;i<=N;i++)
		infile3>>pi[i];

	infile1.close();
	infile2.close();
	infile3.close();
}

void forward_procedure()
{
	int i=0,j=0,t=0;		//varaibles used for loop

	//Initialization Step
	for(i=1;i<=N;i++)
		alpha[1][i] = pi[i]*B[i][obser[1]];


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
			}

			alpha[t+1][j] = sum*B[j][obser[t+1]];
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
	
	//cout<<"the probability is "<<sum<<endl;

}

void EM_procedure()
{
	int i=0,j=0,t=0;						//looping variables
	long double sum=0;						//local sum variable


	//cout<<"Printing Xiiii"<<endl;
	for(i=1;i<=N;i++)

	{
		for(j=1;j<=N;j++)
		{
			for(t=1;t<=T-1;t++)
			{
				sum=0;
				int ii=0,jj=0;
				for(ii=1;ii<=N;ii++)
				{
					for(jj=1;jj<=N;jj++)
					{
						sum += alpha[t][ii]*A[ii][jj]*B[jj][obser[t+1]]*beta[t+1][jj];
					}
				}

				xi[t][i][j] = (alpha[t][i]*A[i][j]*B[j][obser[t+1]]*beta[t+1][j])/sum;
				//cout<<xi[j][i][j]<<"\t";
			}
			//cout<<endl<<endl;
		}
		//cout<<endl;
	}

/*	for(t=1;t<T;t++)
	{
		for(i=1;i<=N;i++)
		{
			for(j=1;j<=N;j++)
			{
				cout<<xi[t][i][j]<<"\t";
			}
			cout<<endl<<endl;
		}
		cout<<endl;
	}
	*/
	for(t=1;t<=T;t++)
	{
		for(i=1;i<=N;i++)
		{
			sum =0;
			for(j=1;j<=N;j++)
			{
				sum += xi[t][i][j];
			}
			gamma[t][i] = sum; 
		}
	}
	//cout<<"Gamma Values"<<endl;
	for(t=1;t<=T;t++)
	{
		for(i=1;i<=N;i++)
		{
			//cout<<gamma[t][i]<<"\t";
		}
		//cout<<endl;
	}
	//cout<<endl<<endl;
	//cout<<"NEW PI\n"<<endl;
	for(i=1;i<=N;i++)
	{
		pi[i] = gamma[1][i];
		//cout<<pi[i]<<" ";
	}
	//cout<<endl<<endl;

	//cout<< "NEW A:"<<endl;
	long double temp =0;
	for(i=1;i<=N;i++)
	{
		for(j=1;j<=N;j++)
		{
			sum =0;
			for(t=1;t<=T-1;t++)
			{
				sum += xi[t][i][j];
			}
			temp =0;
			for(t=1;t<=T-1;t++)
			{
				temp += gamma[t][i];
			}

			A[i][j] = sum/temp;
			//cout<<A[i][j]<<"\t";
		}
		//cout<<endl;
	}
	int k=0;
//	cout<<"New B:"<<endl;
	long double threshold = pow((double)10,(double)-25);
	long double sum2 = 0,max = 0.0, stochastic_check=0;
	int max_index= 0;
	for(j=1;j<=N;j++)
	{
		sum2 = 0;
		for(k=1;k<=M;k++)
		{
			sum =0;
			for(t=1;t<=T;t++)
			{
				if(obser[t] == k)
				{
					sum += gamma[t][j];
				}
			}
			temp =0;
			for(t=1;t<=T;t++)
			{
				temp += gamma[t][j];	
			}
			B[j][k] = sum/temp;

			if(B[j][k] < threshold)
			{
				sum2 += B[j][k];
				B[j][k] = threshold;
			}
			else
			{
				if(B[j][k] > max)
				{
					max = B[j][k];
					max_index = k;
				}
			}
			//cout<<B[j][k]<<"\t";
		}
		//cout<<endl;
		B[j][max_index] -= sum2;
	}

	//cout<<"Stochastic Check in B Matrix"<<endl;
	for(j=1;j<=N;j++)
	{
		stochastic_check = 0;
		for(k=1;k<=M;k++)
		{
			stochastic_check += B[j][k];
		}
		//cout<<"Stochastc_Row"<<j<<": "<<stochastic_check<<endl;
	}
	//cout<<endl;
}
void viterbi() 
{
	long double delta[T_MAX][N+1] = {0};
	int q[T_MAX] = {0};
	long double si[T_MAX][N+1] = {0};

	int i=0,j=0,t=0;//variables used for looping

	//Initialization Step
	for(i=1;i<=N;i++)
	{
		delta[1][i] = pi[i]*B[i][obser[1]];
		si[1][i] = 0;
	}

	

//	for(i=1;i<=N;i++)
//		cout<<si[1][i]<<"\t";
//	cout<<endl<<endl;

	//Recursion Step
	long double max = 0;//variable used to determine the MAX variable
	int index = 0;		//max index
	for(t=2;t<=T;t++)
	{
		for(j=1;j<=N;j++)
		{
			max=0;
			for(i=1;i<=N;i++)
			{
				if(max<delta[t-1][i]*A[i][j])
				{
					max = delta[t-1][i]*A[i][j];
					index = i;
				}
			}
			si[t][j] = index;
			delta[t][j] = max*B[j][obser[t]];
		}
	}

/*	cout<<"delta"<<endl;
	for(j=1;j<=T;j++)
	{
		for(i=1;i<=N;i++)
			cout<<delta[j][i]<<"\t";
		cout<<endl;
	}

	cout<<"si"<<endl;
	for(j=1;j<=T;j++)
	{
		for(i=1;i<=N;i++)
			cout<<si[j][i]<<"\t";
		cout<<endl;
	}

	cout<<endl<<endl;*/
	//Termination Step
	max=0;
	long double P;
	index=0;
	for(i=1;i<=N;i++)
	{
		if(max<delta[T][i])
		{
			max=delta[T][i];
			index = i;
		}
	}
	P = max;
	new_pstar = P;
	q[T] = index;

	//State Sequence -> path backtracking

	for(t=T-1;t>=1;t--)
	{
		q[t] = si[t+1][q[t+1]];
	}
	
	//cout<<endl<<endl<<"State Sequence is:"<<endl;
	for(t=1;t<=T;t++)
	{
		//cout<<q[t]<<" ";
	}
	//cout<<endl;

	//cout<<"P star probability is: "<<P<<endl;

}

void read_observation_sequence(long double i1, long double i2)
{
	int i=0;
	string location;
	location = "..\\HMM_data\\Observation_Sequence_file\\digit";
	location = location + to_string(i1) + "\\184101010_" + to_string(i1) + "_" + to_string(i2) + "_Obs_seq.txt";
	ifstream infile;
	infile.open(location);
	infile >> T;
	for(i=1;i<=T;i++)
	{
		infile >> obser[i];
	}
	infile.close();
}

void average_model(long double i1)
{
	string location1,location2;
	location1 = "..\\HMM_data\\Average_Model\\";
	location2 = location1;
	location1 = location1 + "184101010_" + to_string(i1) + "A.txt";
	location2 = location2 + "184101010_" + to_string(i1) + "B.txt";

	ofstream outfile1,outfile2;
	outfile1.open(location1);
	outfile2.open(location2);
	int i=0,j=0;
	
	for(i=0;i<=N;i++)
	{
		for(j=0;j<=N;j++)
		{
			A_Avg[i][j] /= 20;
			outfile1 << A_Avg[i][j] <<"\t";
		}
		outfile1<<endl;
	}
	for(i=0;i<=N;i++)
	{
		for(j=0;j<=M;j++)
		{
			B_Avg[i][j] /= 20;
			outfile2 << B_Avg[i][j]<<"\t";
		}
		outfile2<<endl;
	}

	outfile1.close();
	outfile2.close();
}

void initial_avg_model(long double i1)
{
	string location1,location2;

	location1 = "..\\HMM_data\\Average_Model\\";
	location2 = location1;
	location1 = location1 + "184101010_" + to_string(i1) + "A.txt";
	location2 = location2 + "184101010_" + to_string(i1) + "B.txt";

	ifstream infile1,infile2;
	infile1.open(location1);
	infile2.open(location2);

	int i=0,j=0;

	for(i=0;i<=N;i++)
	{
		for(j=0;j<=N;j++)
		{
			infile1 >> A[i][j];
		}
	}

	for(i=0;i<=N;i++)
	{
		for(j=0;j<=M;j++)
		{
			infile2 >> B[i][j];
		}
	}

	infile1.close();
	infile2.close();
}

void add_A_B_Matrix()
{
	int i=0,j=0;
	
	for(i=0;i<=N;i++)
	{
		for(j=0;j<=N;j++)
		{
			A_Avg[i][j] += A[i][j];
			//cout<<A_Avg[i][j]<<" ";
		} 
		//cout<<endl;
	}
	for(i=0;i<=N;i++)
	{
		for(j=0;j<=M;j++)
		{
			B_Avg[i][j] += B[i][j];
			//cout<<B_Avg[i][j]<<" ";
		} 
		//cout<<endl;
	}
}

void initialize_A_B_Matrix()
{
	int i=0,j=0;
	
	for(i=0;i<=N;i++)
	{
		for(j=0;j<=N;j++)
		{
			A_Avg[i][j] = 0;
		}
	}
	for(i=0;i<=N;i++)
	{
		for(j=0;j<=M;j++)
		{
			B_Avg[i][j] = 0;
		}
	}
}

int main()
{
	int i=0,j=0;
	initialize_code_book_vector();
	initialize_A_B_Pi_Matrix();
	for(i=0;i<=6;i++)
	{
		for(j=1;j<=20;j++)
		{
			T=0;
			initialize_work_arr();
			fileread(i,j);
			calculate_tokura_distance(i,j);
			storing_A_B_Pi_Matrix(i,j);
		}
	}
	initialize_A_B_Pi_Matrix();
	for(i=0;i<=6;i++)
	{
		for(j=1;j<=20;j++)
		{
			old_pstar = 1.0e-250;
			initialize_A_B_Pi_Matrix();
			read_observation_sequence(i,j);
			int k=0;
			for(k=0;k<40;k++)
			{
				forward_procedure();
				viterbi();
				backward_procedure();
				EM_procedure();
				if((new_pstar/old_pstar) < 1.00005)
					break;
				old_pstar = new_pstar;
			}
			//cout<< "Converging @ iteration: "<<k<<endl;
			//forward_procedure();
		//	viterbi();
			storing_A_B_Pi_Matrix(i,j);
			add_A_B_Matrix();

		}
		T=0;
		average_model(i);
		initialize_A_B_Matrix();
	}
	for(i=0;i<=6;i++)
	{
		initial_avg_model(i);
		for(j=1;j<=20;j++)
		{
			old_pstar = 1.0e-250;
			read_observation_sequence(i,j);
			int k=0;
			for(k=0;k<40;k++)
			{
				forward_procedure();
				viterbi();
				backward_procedure();
				EM_procedure();
				if((new_pstar/old_pstar) < 1.00005)
					break;
				old_pstar = new_pstar;
			}
			//cout<< "Converging @ iteration: "<<k<<endl;
			//forward_procedure();
			//viterbi();
			storing_A_B_Pi_Matrix(i,j);
			add_A_B_Matrix();
			initialize_A_B_Pi_Matrix();
		}
		T=0;
		average_model(i);
		initialize_A_B_Matrix();
	}
}