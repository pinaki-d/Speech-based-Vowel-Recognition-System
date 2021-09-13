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
long double c[totalBlocks][p+1]={0};//to store cepstral coefficients of 5 steady frames after computation of ci
long double tokhuraDistance[5] = {0};//to store tokhura distance
int predictedVowelIndex;

/*Preprocessing the input signal*/
/*Opens the input file and stores the amp values into an array, corrects the DC shift, calculates the normlization factor, performs normalization,
finds out the 5 steady frames and stores those 5 steady frams in an array for further processing*/
void preProcess(string filename)
{
	input.open(filename);
	/*Opening the input file and fetching the amp values & storing the values into array*/
	if(!input)
	{
		cout<<"Could not open the file "<<filename<<"\n";
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
	input.close();

	/*Correcting DC Shift of whole input amplitude values if DCShift value is not equal to 0*/
	if(DCShift != 0.0)
	for(i=0;i<totalSamples;i++)
	{
		inputarr[i] -= DCShift;
	}

	/*Calculating the maximum aplitude value of the speech signal*/
	maxAmp=0;
	for(i=0;i<totalSamples;i++)
	{
		if(maxAmp<abs(inputarr[i]))//finding the maxm amplitude and storing the maximum amplitude that will be used later for normalization
		maxAmp = abs(inputarr[i]);
	}

	/*Normalization of whole input speech signal*/
	long double normFactor = 10000/maxAmp;
	for(i=0;i<totalSamples;i++)
	{
		inputarr[i] *= normFactor;
	}

	/*Applying hamming weights to the samples*/
	for(i=0;i<totalSamples;i++)
	{
		inputarr[i] *= w[i%N];//w[] array contains hamming weights which are being multiplied with the input signal amplitude
	}

	/*Finding STE value of each frame of the input signal*/
	long double energy=0;
	int steIndex = 0;
	for(i=0;i<totalSamples;i++)
	{
		energy = energy +(inputarr[i]*inputarr[i]);
		if((i+1)%320==0)//if no. of samples becomes 320, it means we have traversed one frame and now we should store the STE value of that frame into STE[ ] array
		{
			STE[steIndex] = 0;//initialising STE[ ] array
			energy = energy/N;//dividing the energy value by no. of samples, i.e., 320
			STE[steIndex++] = energy;
			energy = 0;
		}
	}

	/*Finding max STE value for determining 5 steady frames and storing the max STE value , frame number which has the maximum STE in variables maxSTEValue and maxSTEIndex*/ 
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

	/*Considering the interval: 2 frames before and 2 frames after the frame having maximum STE value for calculation of coefficients*/
	int start = (maxSTEIndex-2)*N;//starting sample number of steady frame
	int end = (maxSTEIndex+2)*N + (N-1);//ending sample number of steady  frame
	for(i=start,j=0;i<=end;i++,j++)//copying all the samples of interval under consideration to steadyInput array for further processing
	{
		steadyInput[j] = 0;
		steadyInput[j] = inputarr[i];
	}
}

/*This function omputes the Ri's, ai's and cepstral cofficient ci's for each of the 5 steady frames of an input speech signal*/
void computationOfCoefficients()
{
	/*==== Computation of Ri's ======*/
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
		j+=N;
	}

/*===========Computation of ai's using Durbin's algorithm=============*/	 
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

	/*Inverting the sign of ai's before calculating ci's*/
	/*As discussed with sir in class, I am getting better results without inverting ai's before calculating ci's. So I have commented out this part which inverts the sign
	of ai's*/
	/*n=0;
	while(n<totalBlocks)
	{
		for(j=1;j<=p;j++)
		{
			a[n][j] = a[n][j] * (-1);
		}
		n++;
	}*/

	/*===================================Calculating Cepstral Coefficients for all the frame of steady input file=======================================*/
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
		}
	}

	/*Applying raised sine window weights to calculated cepstral coefficients*/
	for(int n=0;n<totalBlocks; n++)
	{
		for(m=1;m<=p;m++)
		{
			c[n][m] = c[n][m] * (sin_w[m]);
		}
	}

	/*This will add all ci's framewise of each input signal(steady frames) and store in a 2D array named avgCepstral[ ][ ] which will later be dumped into a text file*/
	n=0;
	while(n<totalBlocks)
	{
		for(j=1;j<=Q;j++)
		{
			avgCepstral[n][j] += c[n][j];
		}
		n++;
	}
}

/*This fucntion is used to train the system with the file whose file name is passed as argument to this function*/
void training(string filename)
{
	preProcess(filename);//Preprocessing the input files
	computationOfCoefficients();//computing Ri's, ai's and ci's
}

/*======================================Storing the reference cepstral coefficients into text file=====================================*/
/*After the computation of average of ci's which are stored in a 2D array named avgCepstral, this fucntion is called which dumps the values into a text file.
For every vowel, one text file is created and the name of the file is passed on to this function as a parameter*/
void saveCiToText(string outputFileName)
{
	output.open(outputFileName);
	for(i=0;i<totalBlocks;i++)
	{
		for(j=1;j<=Q;j++)
		{
			avgCepstral[i][j] = avgCepstral[i][j]/10;//taking the average of ci's after adding framewise from all input signal  of a particular vowel
			
			output<<std::fixed << std::setprecision(8)<<avgCepstral[i][j]<<"\t\t\t";
			avgCepstral[i][j] = 0;//resetting it to zero so that it can store reault for the next file input signal correctly
		}
		
		output<<endl;
	}
	output.close();
}

/*=====================Reading from reference file and storing into array for testing the test file from reference cepstral coefficients of vowels===================*/
/*After the system is trained, when we move on to perform testing of program on some input speech signals, we should read the cepstral coefficients of reference vowels 
from the files which were created earlier. This fucntion does that, it reads the reference vowel file and then stores the values back in a 2D array and now this 2D array is 
used to compare the ci's of test vowel and reference vowel*/
void readFromReferenceFile(string refFilePath)
{
	input.open(refFilePath);
	if(!input)
	{
		cout<<"Could not open the file "<<refFilePath<<"\n";
		input.clear();
	}
	int row = 0, column;
	vector<string> cepstral;
	string ci;
	while(getline(input,line))//reading one row from text file containing the ci's of reference vowel
	{
		cepstral.push_back(line);
		for(size_t i = 0; i < cepstral.size(); i++)
		{
			column = 1;
			stringstream ss(cepstral[i]);
			while(ss >> ci)//reading one ci's after other in a particular row
			{
				size_t sz;
				x = stod(ci, &sz);//string to double conversion
				referenceCi[row][column++] = x;//storing the ci into 2D array referenceCi
			}
			row++;
		}
		cepstral.pop_back();
	}
	input.close();
}

/*This function does testing of a vowel. It first calls preProcess() to perform preproceesing of the test file and then it calls computationOfCoefficients() to compute
ci's of the test file. After that it calls readFromReferenceFile() for each reference vowel and compute tokhura distance from that reference vowel. The distance computed
is stored in an array named tokhuraDistance[] like: distance of test vowel from 'a' is stored in tokhuraDistance[0], from 'e' is stored in tokhuraDistance[1], from 'i' is
stored in tokhuraDistance[2] and so on.*/
void testing(string testfile)
{
	preProcess(testfile);//preprocessing of test file
	computationOfCoefficients();//computations of coefficients of test vowel
	char ref_vowel[] = {'a','e','i','o','u'};
	string refFile = "cepstral_";//filename of reference file starts with this prefix
	for(i=0;i<5;i++)
	{
		tokhuraDistance[i] = 0;
		string refFilePath = refFile + ref_vowel[i] + "_ref.txt";//forming the file name of reference file so that it can be passed to function readFromReferenceFile( )
		readFromReferenceFile(refFilePath);
		//cout<<"In Testing for vowel: "<<ref_vowel[i]<<endl;
		for(j=0;j<totalBlocks;j++)
		{
			for(int k=1; k<=p;k++)
			{
				tokhuraDistance[i] = (long double)tokhuraDistance[i] + ((float)tokhura_w[k-1])*((long double)c[j][k]-(long double)referenceCi[j][k])*((long double)c[j][k]-(long double)referenceCi[j][k]);
			}
		}
	}
	//printing the tokhura's distance of test vowel from each of reference vowel
	for(i=0;i<5;i++)
	{
		printf("Distance from %c is %6lf\n",ref_vowel[i], tokhuraDistance[i]);
	}
	/*Finding the minimum value of tokhura distance*/
	double long minTokhuraDistance = tokhuraDistance[0];
	predictedVowelIndex = 0;
	for(i=0;i<5;i++)
	{
		if(minTokhuraDistance>tokhuraDistance[i])
		{
			minTokhuraDistance = tokhuraDistance[i];
			predictedVowelIndex = i;//the vowel from which tokhura distance is minimum is stored using predictedVowelIndex. If it is least distant from a, predictedVowelIndex=0
			//, if least distant from e then predictedVowelIndex=1, if least distant from i then predictedVowelIndex=2, and so on.
		}
	}
}
/*The main() function. Program execution starts from here*/
int _tmain(int argc, _TCHAR* argv[])
{
	silence.open("silence.txt");//opening the silence file to calculate DC Shift value
	int silenceSampleCount = 0;//to store the number of samples in silence file
	if(!silence)
	{
		cout<<"Could not open the silence text file\n";
		silence.clear();
	}

	/*Calculating DC Offset value by reading a silence file*/
	while(getline(silence,silenceAmp))
	{
		size_t sz;
		x = stod(silenceAmp, &sz);
		DCShift += x;
		silenceSampleCount++;//counting number of samples in silence file
	}
	silence.close();

	DCShift = DCShift/silenceSampleCount; //will add this DCShift value to all samples of input speech to correct DC shift of input speech signal

	/*Hamming window calculation*/
	for(i=0;i<N;i++)
	{
		w[i] = 0.54 - (0.46*cos((2*PI*i)/(N-1)));
	}

	/*Raised sine window calculation*/
	for(i=1;i<=Q;i++)
	{
		sin_w[i] = 1 + (((float)Q/2)*(sin((PI*i)/Q)));
	}
	/*storing the common prefix of filenames of input and reference file in variables in_path and out_path*/
	string filename="",outputfile="", in_path = "204101042_", out_path = "cepstral_";
	char vowels[] = {'a','e','i','o','u'};//will append this to prefix of file name 
	char filenumber[] = {'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'};//will append this to the prefix of file name
	int filecount,vowelcount;

	/*==================================Populating the training data set===========================*/
	/*this part opens training file one by one and calls training() for each training file to train the system and then calls saveCiToText() for each vowel to dump the
	reference ci's values into a text file*/
	cout<<"Please wait while the program is being trained..."<<endl;
	for(vowelcount=0;vowelcount<5;vowelcount++)
	{
		for(filecount=1;filecount<=10;filecount++)
		{
			if(filecount == 10)
			filename = in_path + vowels[vowelcount] + "_" + filenumber[filecount/10] + filenumber[filecount%10] + ".txt";
			else filename = in_path + vowels[vowelcount] + "_" + filenumber[filecount] + ".txt";
			training(filename);
		}
		outputfile = out_path + vowels[vowelcount] + "_ref.txt";

		saveCiToText(outputfile);//dumps the ci's into a text file whose filename is saved in outputfile variable
	}
	/*==================================Testing the data set===========================*/
	while(true)
	{
		int accuracy = 0;
		int option;
		cout<<"Enter any one option: "<<endl;
		cout<<"1. Testing some existing file(s) choosen randomly from a set of 50 test files (10 files for each vowel)."<<endl;
		cout<<"2. Testing all existing file."<<endl;
		cout<<"3. Real time testing."<<endl;
		cout<<"4. Exit."<<endl;
		cout<<endl;
		cout<<"Your choice : ";
		cin>>option;
		switch(option)
		{
			case 1: //Testing some files by seleting randomly
				
				cout<<"\nEnter no. of files to be chosen randomly for testing : ";
				int f_count;
				cin>>f_count;
				
				for(filecount=0;filecount<f_count;filecount++)
				{
					int randomVowelIndex = rand() % 5;//generating a index for vowel selection randomly
					int randomFileNumber = rand() % 11 + 10;//generating a index number for  randomly between 11 and 20 
					filename = in_path + vowels[randomVowelIndex] + "_" + filenumber[randomFileNumber/10] + filenumber[randomFileNumber%10] + ".txt";//forming the test file name
					cout<<"File tested = "<<filename<<endl;
					testing(filename);
					cout<<"Actual    = "<<vowels[randomVowelIndex]<<endl;//printing the vowel selected for testing
					cout<<"Predicted = "<<vowels[predictedVowelIndex]<<endl;//printing the vowel based on minimum tokhura distance obtained from all reference vowels
					cout<<"\n\n";
					if(randomVowelIndex == predictedVowelIndex)
						accuracy++;//if actual==predicted, increment the accuracy
				}
				cout<<"Accuracy = "<<accuracy<<"/"<<f_count<<endl<<endl;//printing the accuracy stats
				break;
			case 2:
				for(vowelcount=0;vowelcount<5;vowelcount++)
				{
					for(filecount=11;filecount<=20;filecount++)
					{
						filename = in_path + vowels[vowelcount] + "_" + filenumber[filecount/10] + filenumber[filecount%10] + ".txt";
						cout<<"Testing the file "<<filename<<endl;
						testing(filename);
						cout<<"Actual   = "<<vowels[vowelcount]<<endl;//printing the vowel selected for testing
						cout<<"Predicted = "<<vowels[predictedVowelIndex]<<endl;//printing the vowel based on minimum tokhura distance obtained from all reference vowels
						cout<<"\n\n";
						if(vowelcount == predictedVowelIndex)
						accuracy++;//if actual==predicted, increment the accuracy
					}
					cout<<endl<<endl;
				}
				cout<<"Accuracy = "<<accuracy<<"/ 50"<<endl<<endl;//printing the accuracy stats
				break;
			case 3:
				system("Recording_Module.exe 3 input_file.wav input_file.txt");//calling the recording module to receive input speech
				filename = "input_file.txt";
				testing(filename);//calling the testing() and passing the input speech file
				cout<<"\nPredicted = "<<vowels[predictedVowelIndex]<<endl<<endl;//displaying the predicted vowel on the basis of minimum tokhura distance
				break;
			case 4:
				exit(1);

			default:
				cout<<"Invalid option entered"<<endl;
				break;
		}
	}

	return 0;
}

