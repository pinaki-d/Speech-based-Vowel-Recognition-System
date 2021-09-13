========================================================================
    CONSOLE APPLICATION : Vowel Recognition Project Overview
========================================================================

The input files (used for training) are provided in this project. The utterance 1-10 for each vowel are being used to train the program. 
So total 50 files are used for training.

The files used for testing are the ones with utterance numbered 11-20 for each vowel.

The program uses some functions whose funtionalities are decribed briefly in the following:
1. preProcess()-It opens the input speech file, corrects DC shift, performs normalization, finds out the steady frames 
   (5 frames in total for each input file).
2. computationOfCoefficients() - computes the ri's, ai's and ci's. It also applies raised sine window to the calculated ci's.
3. training() - trains the system by calling above 2 functions in order.
4. saveCiToText() -dumps the average ci values of a vowel into text file.
5. readFromReferenceFile() - reads ci values from the above text file and forms a 2D array containing the reference cepstral coefficients
6. testing() - calls the first 2 functions for preprocessing and computation of coefficients of test vowel. Then for each vowel, it calls 
   5th function to form the reference ci's matrix, finds tokhura distance from that reference vowel. Finally it finds the minimum tokhura 
   distance among all 5 reference vowels.


The user is given 3 choices to select from:
1. Input no of files to be tested and test that many files selected randomly by using rand() function to select any vowel and any 
   utterance number between 11-20. Finally it shows the accuracy value.
2. to test all the existing test files (total 50 files; 10 files for each vowel).
3. real time testing where speech is recorded by calling the recording module from within the program.

After training the program, 5 text files are generated namely cepstral_a_ref, cepstral_e_ref, cepstral_i_ref, cepstral_o_ref, cepstral_u_ref
These files contains the average cepstral coeffiients of for each steady frame of each vowel.
These files are used as reference files while testing any vowel.