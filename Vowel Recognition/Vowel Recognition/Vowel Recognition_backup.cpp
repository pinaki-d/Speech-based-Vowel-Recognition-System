// Vowel Recognition.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdio.h"
#include "fstream"
#include "iostream"
#include "conio.h"
#include "string"
#include "math.h"
#include "vector"
#include "sstream"
#include "stdlib.h"
#include "iomanip"
using namespace std;

#define totalBlocks 5
#define p 12
#define Q 12
#define N 320
#define PI 3.14

long double inputarr[90000],x,DCShift=0, maxAmp;
long double steadyInput[1605], STE[50000];//steadyInput array contains amp value of 5 steady frames for each input signal text file
string line,silenceAmp;
ifstream silence,input;
ofstream output;
int m,i,j=0,index=0,n;
long double w[N] = {0}, sin_w[Q+1] = {0} ;//hamming window and raised sine window
long double tokhura_w[12]={1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};//tokhura weights provided in the question
long double avgCepstral[totalBlocks][p+1] = {0};//storing average of ci's for each input file oncthe basis of frame wise
long double referenceCi[totalBlocks][p+1] = {0};//fetching the avg ci's from reference text file and storing in the array
long double c[totalBlocks][p+1]={0};
long double tokhuraDistance[5] = {0};
int predictedVowelIndex;

/*Preprocessing the input signal*/
/*Opens theinput file and stores the amp values into an array, corrects the DC shift, calculates the normlisation factor, performs normalization,
finds out the 5 steady frames and stores those 5 steady frams in an array for further processing*/
void preProcess(string filename)
{
	//cout<<endl<<filename<<endl;
	input.open(filename);
	/*Opening the input file and fetching the amp values & storing the values into array*/
	if(!input)
	{
		cout<<"Could not open the file\n";
		input.clear();
	}

	index=0;
	while(getline(input,line))
	{
		inputarr[index] = 0;//initializing the inputarr to zero
		size_t sz;
		x = stod(line, &sz);//string to double conversion
		inputarr[index++] = x;//storing the amplitudes from text file to array
		
	}
	int totalSamples = index;
	//cout<<"Total samples = "<<totalSamples<<endl;
	input.close();

	/*Correcting DC Shift of whole input amplitude values*/
	if(DCShift != 0.0)
	for(i=0;i<totalSamples;i++)
	{
		inputarr[i] -= DCShift;
	}

	/*Calculating the maximum aplitude value of the speech signal*/
	maxAmp=0;
	for(i=0;i<totalSamples;i++)
	{
		if(maxAmp<abs(inputarr[i]))//storing the maximum amplitude that will be used later for normaliation
		maxAmp = abs(inputarr[i]);
	}
	//cout<<"Before normalization, maxm amp = "<<maxAmp<<endl;

	/*Normalization of whole input speech signal*/
	long double normFactor = 10000/maxAmp;
	//cout<<"Normalization factor = "<<normFactor<<endl;
	for(i=0;i<totalSamples;i++)
	{
		inputarr[i] *= normFactor;
	}

	/*Applying hamming weights to the samples*/
	for(i=0;i<totalSamples;i++)
	{
		inputarr[i] *= w[i%N];
	}
	/*Finding STE value of each frame of the input signal*/
	long double energy=0;
	int steIndex = 0;
	for(i=0;i<totalSamples;i++)
	{
		energy = energy +(inputarr[i]*inputarr[i]);
		if((i+1)%320==0)
		{
			STE[steIndex] = 0;
			energy = energy/N;
			STE[steIndex++] = energy;
			energy = 0;
		}
	}

	/*Finding max STE value for determining 5 steady frames*/
	int maxSTEIndex = 0;
	long double maxSTEValue = 0;
	for(i=0;i<steIndex;i++)
	{
		if(STE[i]>maxSTEValue)
		{
			maxSTEValue = STE[i];
			maxSTEIndex = i;
		}
	}

	/*Consider the interval: 2 frames before and 2 frames after the frame having maximum STE value for calculation of coefficients*/
	int start = (maxSTEIndex-2)*N;//starting sample number of steady frame
	int end = (maxSTEIndex+2)*N + (N-1);//ending sample number of steady  frame
	//cout<<"Start sample no. of steady input : "<<start<<endl;
	//cout<<"End sample no. of steady input : "<<end<<endl;
	for(i=start,j=0;i<=end;i++,j++)//copying all the samples of interval under consideration to steadyInput array for further processing
	{
		steadyInput[j] = 0;
		steadyInput[j] = inputarr[i];
		//steadyInput[j] = steadyInput[j]*w[i%N];//multiplying hamming weight with each amplitude of steady input
	}
}

void computationOfCoefficients()
{
	/*-------------------------------------------Computation of Ri's, ai's and ci's---------------------------------------------------------*/
	j=0;
	long double r[totalBlocks][p+1]={0};
	for(i=0; i<totalBlocks; i++)
	{
		int q =0;
		while(q<=p)
		{
			for(m=0; m <= N-1-q ; m++)
			{
				r[i][q] += steadyInput[j+m]*steadyInput[j+m+q];
			}
			q++;
		}
		//cout<<"j value for frame "<<i<<" is : "<<j<<endl;
		j+=N;
	}
	//cout<<"\nValues of Ri's are: \n";
	for(i=0;i<totalBlocks;i++)
	{
		//cout<<"Frame : "<<i<<endl;
		for(j=0;j<=p;j++)
		{
			//printf("R%d : %6lf\n",j, r[i][j]);
		}
	}
	 
	long double E[p+1] = {0};
	long double alpha[p+1][p+1]={0};
	long double a[totalBlocks][p+1]={0};
	long double k[p+1] = {0};
	long double threshold = 0;
	 
	for(n=0;n<totalBlocks; n++)
	{
		if(r[n][0] > threshold)
	 	{
	 		E[0] = r[n][0];
	 		for (i=1; i<=p; i++)
	 		{
	 			long double summationTerm = 0;
	 			for(j=1;j<=i-1;j++)
	 			{
	 				summationTerm += alpha[i-1][j] * r[n][i-j];
				}
				k[i] = (r[n][i] - summationTerm) / E[i-1];
				alpha[i][i] = k[i];
				
				if(i>1)
				{
					for(j=1;j<=i-1;j++)
					{
						alpha[i][j] = alpha[i-1][j] - (k[i]*alpha[i-1][i-j]);
					}
				}
				E[i] = (1-(k[i]*k[i]))*E[i-1];
			}
			
			for(j=1;j<=p;j++)
			{
				a[n][j] = alpha[p][j];
			}
		}
		
	}

	//cout<<"The AUTOCORELATION COEFFICIENTS are: \n";
	n=0;
	while(n<totalBlocks)
	{
		//cout<<"For frame "<<n<< " : "<<endl;
		for(j=1;j<=p;j++)
		{
			//printf("a%d : %6lf\n",j,a[n][j]);
		}
		n++;
	}

	/*Inverting the sign of ai's before calculating ci's*/
	n=0;
	while(n<totalBlocks)
	{
		for(j=1;j<=p;j++)
		{
			a[n][j] = a[n][j] * (-1);
		}
		n++;
	}




	/*=================================Calculating Cepstral Coefficients for all the frame of steady input file=======================================*/
	for(int n=0;n<totalBlocks; n++)
	{
		c[n][0] = log(r[n][0]*r[n][0]);
		for(m=1;m<=p;m++)
		{
			c[n][m] = 0;
			int k;
			long double summationTerm = 0;
			for(k=1;k<=m-1;k++)
			{
				summationTerm += ((float)k/m)*c[n][k]*a[n][m-k];
			}	
			c[n][m] = a[n][m] + summationTerm;
			c[n][m] = c[n][m] * (sin_w[m]);
		}
	}

	//cout<<"The CEPSTRAL COEFFICIENTS are: \n";
	n=0;
	while(n<totalBlocks)
	{
		//cout<<"For frame "<<n<< " : "<<endl;
		for(j=1;j<=p;j++)
		{
			//printf("c%d before sine wave : %6lf\n",j,c[n][j]);
			//c[n][j] = c[n][j] * (sin_w[j]);
			//printf("c%d : %6lf\n",j,c[n][j]);
			//cout<<"a"<<j<<"	"<<a[n][j]<<endl;
			//printf("c%d after sine wave : %6lf\n",j,c[n][j]);
		}
		n++;
	}

	/*This will add all ci's framewise of each input signal*/
	n=0;
	while(n<totalBlocks)
	{
		for(j=1;j<=Q;j++)
		{
			avgCepstral[n][j] += c[n][j];
			//printf("%6lf\t",avgCepstral[n][j]);
		}
		//printf("\n");
		n++;
	}
}

/*==============================Training the data set========================================*/
void training(string filename)
{
	preProcess(filename);
	computationOfCoefficients();
}

/*======================================Storing the reference cepstral coefficients into text file=====================================*/
void saveCiToText(string outputFileName)
{
	output.open(outputFileName);
	for(i=0;i<totalBlocks;i++)
	{
		for(j=1;j<=Q;j++)
		{
			avgCepstral[i][j] = avgCepstral[i][j]/10;//taking the average of ci's after adding framewise from all input signal  of a particular vowel
			
			output<<std::fixed << std::setprecision(8)<<avgCepstral[i][j]<<"\t\t\t";
			avgCepstral[i][j] = 0;
		}
		
		output<<endl;
	}
	output.close();
}

void readFromReferenceFile(string refFilePath)
{
	input.open(refFilePath);
	if(!input)
	{
		cout<<"Could not open the file\n";
		input.clear();
	}
	int row = 0, column;
	vector<string> cepstral;
	string ci;
	while(getline(input,line))
	{
		cepstral.push_back(line);
		for(size_t i = 0; i < cepstral.size(); i++)
		{
			column = 1;
			stringstream ss(cepstral[i]);
			while(ss >> ci)
			{
				size_t sz;
				x = stod(ci, &sz);//string to double conversion
				referenceCi[row][column++] = x;
			}
			row++;
		}
		cepstral.pop_back();
	}
	input.close();
	/*printf("Reference vowel cepstral matrix:\n");
	for(int i1=0;i1<5;i1++)
	{
		for(int j1=1;j1<=12;j1++)
		{
			printf("%6lf\t", referenceCi[i1][j1]);
			//cout<<referenceCi[i][j]<<std::setprecision(8);
		}
		cout<<"\n";
	}*/
}
void testing(string testfile)
{
	preProcess(testfile);
	computationOfCoefficients();
	char ref_vowel[] = {'a','e','i','o','u'};
	string refFile = "cepstral_";
	for(i=0;i<5;i++)
	{
		tokhuraDistance[i] = 0;
		string refFilePath = refFile + ref_vowel[i] + "_ref.txt";
		readFromReferenceFile(refFilePath);
		//cout<<"In Testing for vowel: "<<ref_vowel[i]<<endl;
		for(j=0;j<totalBlocks;j++)
		{
			for(int k=1; k<=p;k++)
			{
				tokhuraDistance[i] = (long double)tokhuraDistance[i] + ((float)tokhura_w[k-1])*((long double)c[j][k]-(long double)referenceCi[j][k])*((long double)c[j][k]-(long double)referenceCi[j][k]);
				//cout<<"Ci : "<<c[j][k]<<endl;
				//cout<<"reference Ci : "<<referenceCi[j][k]<<endl;
			}
		}
	}
	//debugging
	/*for(i=0;i<5;i++)
	{
		printf("Distance from %c is %6lf\n",ref_vowel[i], tokhuraDistance[i]);
	}*/
	double long minTokhuraDistance = tokhuraDistance[0];
	predictedVowelIndex = 0;
	for(i=0;i<5;i++)
	{
		if(minTokhuraDistance>tokhuraDistance[i])
		{
			minTokhuraDistance = tokhuraDistance[i];
			predictedVowelIndex = i;
		}
	}
}
int _tmain(int argc, _TCHAR* argv[])
{
	silence.open("silence.txt");//opening the silence file to calculate DC Shift value
	int silenceSampleCount = 0;//to store the number of samples in silence file
	if(!silence)
	{
		cout<<"Could not open the silence text file\n";
		input.clear();
	}

	/*Calculating DC Offset value by reading a silence file*/
	while(getline(silence,silenceAmp))
	{
		size_t sz;
		x = stod(silenceAmp, &sz);
		DCShift += x;
		silenceSampleCount++;
	}
	silence.close();
	//cout<<silenceSampleCount<<endl;

	DCShift = DCShift/silenceSampleCount; //will add this DCShift value to all samples of input speech
	//printf("DC Shift value = %6lf\n",DCShift);

	/*Hamming window calculation*/
	for(i=0;i<N;i++)
	{
		w[i] = 0.54 - (0.46*cos((2*PI*i)/(N-1)));
	}

	/*Raised sine window calculation*/
	for(i=1;i<=Q;i++)
	{
		sin_w[i] = 1 + (((float)Q/2)*(sin((PI*i)/Q)));
		//printf("sin_w at %d : %6lf\n",i,sin_w[i]);
	}

	string filename="",outputfile="", in_path = "204101038_", out_path = "cepstral_";
	char vowels[] = {'a','e','i','o','u'};
	char filenumber[] = {'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'}; 
	int filecount,vowelcount;

	/*==================================Populating the training data set===========================*/
	/*for(vowelcount=0;vowelcount<5;vowelcount++)
	{
		for(filecount=1;filecount<=10;filecount++)
		{
			if(filecount == 10)
			filename = in_path + vowels[vowelcount] + "_" + filenumber[filecount/10] + filenumber[filecount%10] + ".txt";
			else filename = in_path + vowels[vowelcount] + "_" + filenumber[filecount] + ".txt";
			training(filename);
		}
		outputfile = out_path + vowels[vowelcount] + "_ref.txt";

		saveCiToText(outputfile);
	}*/
	/*==================================Testing the data set===========================*/
	for(vowelcount=0;vowelcount<5;vowelcount++)
	{
		for(filecount=11;filecount<=20;filecount++)
		{
			filename = in_path + vowels[vowelcount] + "_" + filenumber[filecount/10] + filenumber[filecount%10] + ".txt";
			//cout<<"Testing for "<<vowels[vowelcount]<<" using file  "<<filename<<endl;
			testing(filename);
			cout<<"Actual   = "<<vowels[vowelcount]<<endl;
			cout<<"Predited = "<<vowels[predictedVowelIndex]<<endl;
		}
		cout<<endl<<endl;
	}
	/*for(filecount=0;filecount<20;filecount++)
	{
		int randomVowelIndex = rand() % 5;
		int randomFileCount = rand() % 11 + 10;
		filename = in_path + vowels[randomVowelIndex] + "_" + filenumber[randomFileCount/10] + filenumber[randomFileCount%10] + ".txt";
		testing(filename);
		cout<<"File tested = "<<filename<<endl;
		cout<<"Uttered   = "<<vowels[randomVowelIndex]<<endl;
		cout<<"Predicted = "<<vowels[predictedVowelIndex]<<endl;
		cout<<"\n\n";
	}*/

	return 0;
}

