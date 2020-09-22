#include <math.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <time.h>
#define _USE_MATH_DEFINES
using namespace std;

//all variables
float mean = 0;
float sum = 0;
float numLines = 0;
float varNum = 0;
string line = "";
float variance = 0;
float stDev = 0;
bool analyzing = true;
char nuc1;
char nuc2;
int AA = 0;
int AC = 0;
int AT = 0;
int AG = 0;
int CC = 0;
int CA = 0;
int CT = 0;
int CG = 0;
int TT = 0;
int TA = 0;
int TC = 0;
int TG = 0;
int GG = 0;
int GA = 0;
int GC = 0;
int GT = 0;
int numA = 0;
int numC = 0;
int numT = 0;
int numG = 0;
float bigrams = 0;
float a;
float b;
float C;
float D ;
int reqA;
int reqC;
int reqT;
int reqG;
int selection;
ofstream outputStream;
ifstream inputStream;

//This method will calculate the mean and then will output the result into a file
int dnaMean(){
  mean = sum/numLines;
  //Output results to a file
  outputStream << "Mean: " << mean << endl;
  return mean;
}


//Creating a method that calculates the variance and standard deviation
int dnaVar(string fileName){
  inputStream.open(fileName);
  ifstream inputFile;
  inputFile.open(fileName);
  //We are iterating through the file and also calculating the numerator of the variance formula
  while(getline(inputFile, line)){
    varNum += pow((int)line.size() - mean, 2);
  }
  //now we calculate the full variance and standard deviation and output the results
  variance = varNum / (numLines - 1);
  outputStream << "Variance: " << variance << endl;
  return variance;
  inputFile.close();
}
//this is a method for calculating the standard deviation
int dnaStdev(){
  stDev = sqrt(variance);
  outputStream << "Standard Deviation: " << stDev << endl;
  return stDev;
}


//this is a method for calculating the guassion value
float dnaGaussian(){
  //generating the random numbers
  float a = ((float) rand())/(RAND_MAX);
  float b = ((float) rand())/(RAND_MAX);
  //using  Box-Muller transform to compute a random variable
  float C = (sqrt(-2 * log(a))) * (cos(2* M_PI * b));
  //converting the standard Gaussian to a normal random variable
  float D = sqrt(variance) * C + mean;
  return D;
}


//here is our main method
int main(int argc, char ** argv){
  //file name being accepted as a command line argument
  string fileName = argv[1];
  //here we are opening the file in the inputStream
  inputStream.open(fileName);
  //making the output of the desciptive statistics mattmead.out
  outputStream.open("mattmead.out");

  //we will now check to make sure the inputStream file was done correctly
  if(!inputStream){
    cerr << "Could not open file: " << fileName << endl;
    //exiting the file due to the failure of a correct file
    exit(1);
  }
  //This will constantly change the seed of the random numbers while keeping it randomized between [0,1)
  srand((unsigned)time(NULL));
  //here i am making it so that the user can continue analyzing files until they choose to stop
  while(analyzing){
    //going through each line of the inputed file
    while(getline(inputStream, line)){
      //incrementing numLines in order to use for future calculations
      numLines++;
      //we will now iterate through each letter in a line
      int i = 0;
      while(i < line.size()){
        nuc1 = toupper(line[i]);
        nuc2 = toupper(line[i+1]);
        //checking all the possible bigrams that start with nucleotide A
        if(nuc1 == 'A'){
          numA++;
          //switch statement to see what bigram is being formed
          switch(nuc2){
            case 'A':
              AA++;
              bigrams++;
              break;

            case 'C':
              AC++;
              bigrams++;
              break;

            case 'T':
              AT++;
              bigrams++;
              break;

            case 'G':
              AG++;
              bigrams++;
              break;

            default:
              break;
          }
        }
        //checking all the possible bigrams that start with nucleotide C
        if(nuc1 == 'C'){
          numC++;
          //switch statement to see what bigram is being formed
          switch(nuc2){
            case 'C':
              CC++;
              bigrams++;
              break;

            case 'A':
              CA++;
              bigrams++;
              break;

            case 'T':
              CT++;
              bigrams++;
              break;

            case 'G':
              CG++;
              bigrams++;
              break;

            default:
              break;
          }
        }
        //checking all the possible bigrams that start with nucleotide T
        if(nuc1 == 'T'){
          numT++;
          //switch statement to see what bigram is being formed
          switch(nuc2){
            case 'T':
              TT++;
              bigrams++;
              break;

            case 'A':
              TA++;
              bigrams++;
              break;

            case 'C':
              TC++;
              bigrams++;
              break;

            case 'G':
              TG++;
              bigrams++;
              break;

            default:
              break;
          }
        }
        //checking all the possible bigrams that start with nucleotide G
        if(nuc1 == 'G'){
          numG++;
          //switch statement to see what bigram is being formed
          switch(nuc2){
            case 'G':
              GG++;
              bigrams++;
              break;

            case 'A':
              GA++;
              bigrams++;
              break;

            case 'C':
              GC++;
              bigrams++;
              break;

            case 'T':
              GT++;
              bigrams++;
              break;

            default:
              break;
          }
        }
        i++;
      }
      //we are tracking how many nucleotides are in each DNA string
      sum += line.size();
    }

    //now that we have gone through every line in the input file, we can close it up
    inputStream.close();

    //We will start to write out the results to the outputStream file now that we have all of the information need
    outputStream << "Matthew Pascual-Mead" << endl;
    outputStream << "2300074" << endl;
    outputStream << "CPSC-350-02" << endl;
    outputStream << "Prof. Rene German" << endl;
    outputStream << "Assignment 1" << endl;
    outputStream << endl;

    //Lets start outputting some of the summary statistics
    dnaMean();
    dnaVar(fileName);
    dnaStdev();
    outputStream << "Gaussian Value: " << dnaGaussian() << endl;
    outputStream << endl;

    //nucleotide probability calculations
    float probA = (numA/sum);
    float probC = (numC/sum);
    float probT = (numT/sum);
    float probG = (numG/sum);

    //now lets do nucleotide probabilities
    outputStream << "Probability  of nucleotide A: " << probA * 100 << "%" << endl;
    outputStream << "Probability  of nucleotide C: " << probC * 100 << "%" << endl;
    outputStream << "Probability  of nucleotide T: " << probT * 100 << "%" << endl;
    outputStream << "Probability  of nucleotide G: " << probG * 100 << "%" << endl;
    outputStream << endl;

    //probabilities for the bigrams starting with A
    outputStream << "Probability  of bigram AA: " << AA / bigrams * 100 << "%" << endl;
    outputStream << "Probability  of bigram AC: " << AC / bigrams * 100 << "%" << endl;
    outputStream << "Probability  of bigram AT: " << AT / bigrams * 100 << "%" << endl;
    outputStream << "Probability  of bigram AG: " << AG / bigrams * 100 << "%" << endl;
    outputStream << endl;

    //probabilities for the bigrams starting with C
    outputStream << "Probability  of bigram CC: " << CC / bigrams * 100 << "%" << endl;
    outputStream << "Probability  of bigram CA: " << CA / bigrams * 100 << "%" << endl;
    outputStream << "Probability  of bigram CT: " << CT / bigrams * 100 << "%" << endl;
    outputStream << "Probability  of bigram CG: " << CG / bigrams * 100 << "%" << endl;
    outputStream << endl;

    //probabilities for the bigrams starting with T
    outputStream << "Probability  of bigram TT: " << TT / bigrams * 100 << "%" << endl;
    outputStream << "Probability  of bigram TA: " << TA / bigrams * 100 << "%" << endl;
    outputStream << "Probability  of bigram TC: " << TC / bigrams * 100 << "%" << endl;
    outputStream << "Probability  of bigram TG: " << TG / bigrams * 100 << "%" << endl;
    outputStream << endl;

    //probabilities for the bigrams starting with G
    outputStream << "Probability  of bigram GG: " << GG / bigrams * 100 << "%" << endl;
    outputStream << "Probability  of bigram GA: " << GA / bigrams * 100 << "%" << endl;
    outputStream << "Probability  of bigram GC: " << GC / bigrams * 100 << "%" << endl;
    outputStream << "Probability  of bigram GT: " << GT / bigrams * 100 << "%" << endl;
    outputStream << endl;

    /*this is appending the 1000 DNA strings that follow the Guassian distribution
    while also following the same relative frequency of nucleotides
    **/
    int l = 0;
    while(l < 1000){
      D = dnaGaussian();

      //the string we will append nucleotides to
      string appendNuc = "";

      //A  nucleotides required
      reqA = D * probA;
      //Appending the correct amount of A nucleotides
      for(int q = 0; q < reqA; ++q){
        appendNuc += "A";
      }

      //C  nucleotides required
      reqC = D * probC;
      //Appending the correct amount of C nucleotides
      for(int w = 0; w < reqC; ++w){
        appendNuc += "C";
      }

      //T  nucleotides required
      reqT = D * probT;
      //Appending the correct amount of T nucleotides
      for(int e = 0; e < reqT; ++e){
        appendNuc += "T";
      }

      //G  nucleotides required
      reqG = D * probG;
      //Appending the correct amount of G nucleotides
      for(int r = 0; r < reqG; ++r){
        appendNuc += "G";
      }
      outputStream << appendNuc << endl;
      l++;
    }

    cout << "DNA was analyzed! Would you like to analyze another file? Yes(1) or No(2): " << endl;
    //asking for the users input
    cin >> selection;
    //if they press 1 they will be prompted to enter another file
    if(selection == 1){
      cout << "Enter the name of the file you would like to analyze next. " << endl;
      cin >> fileName;
      ifstream inputStream;
      inputStream.open(fileName);
      //if the file does not exist it will once again exit the program
      if(!inputStream){
        cerr << "Could not open file: " << fileName << endl;
        exit(1);
      }
    }
    //if the user enter a number other than one it will exit the loop and the program will end
    else{
      //will set the while loop to false stopping the analyzing of anymore DNA
      analyzing = false;
    }
    outputStream << endl;
    outputStream << endl;
    //resetting all of the variables back to the original values so that that we can analyze multiple files
    line = "";
    sum = 0;
    numLines = 0;
    numA = 0;
    numC = 0;
    numT = 0;
    numG = 0;
    bigrams = 0;
    AA = 0;
    AC = 0;
    AT = 0;
    AG = 0;
    CC = 0;
    CA = 0;
    CT = 0;
    CG = 0;
    TT = 0;
    TA = 0;
    TC = 0;
    TG = 0;
    GG = 0;
    GA = 0;
    GC = 0;
    GT = 0;
    varNum = 0;
  }
  //closing the output stream
  outputStream.close();
}
