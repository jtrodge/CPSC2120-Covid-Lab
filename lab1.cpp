/*
Jason Rodgers
CPSC 2121
jtrodge@clemson.edu
Brian Dean
*/

#include <iostream>
#include <string>
#include <fstream>
//Runtime: O( ); N : # of sequnces in the input; M : average lenght of a sequence; K : lenghth of probe sequence

int main() {
  int N = 0; //Num of lines
  std::string label; //labels in txt file
  std::string dnaseq; //dna sequence in txt file
  int delta = 0; //Num of delta sequnces in our sequence
  int not_delta = 0; //Num of none delta sequences in our sequence

    //Create/Open File
    std::ifstream inputfile("covid-smaller.txt"); // head -n 20 covid.txt > covid-smaller.txt

    //Read in total number of lines in txt file
    while (inputfile >> label >> dnaseq) {
      N++;
    }
    bool * is_delta = new bool[N];
    std::string * sequence = new std::string[N];

    //Close file and reopen
    inputfile.clear();
    inputfile.seekg(0, std::ios_base::beg);
    //Loop to find all delta variants
    for (int i = 0; i < N; i++) {
      inputfile >> label >> sequence[i];
      is_delta[i] = label == "delta_variant";
      if (is_delta[i]) {
        delta++;
      }
      else {
        not_delta++;
      }
    }

    //Close file and reopen
    inputfile.clear();
    inputfile.seekg(0, std::ios_base::beg);

    int K = 100; //lenghth of probe sequence - k
    std::string result, tempResult, tempString;
    bool match = false; //true if it is not a match, false if it is a match
    int tempNum = 0; //Keeps track if the dna seq and probe is not a match
    int falseP = 0; //Num of false positives
    int tempFP = 0;
    int falseN = 0; //Num of false
    int tempFN = 0;
    double FNR = 0; //False Negative rate
    double FPR = 0; //False Positive rate
    double errorR = 0; //Error Rate
    double tempER = 10000; //Temp Error Rate to find best K-probe

    //Compare Probe to the best K-Probe Sequence
    for (int i = 0; i < N; i++) {
      if (!is_delta[i]) {
        continue;
      }

      //Use substr to find the find all of the best probe sequences
      for (unsigned int j = 0; j < sequence[i].length() - K; j++) {
        result = sequence[i].substr(j, K);

        //Reinitialization during loops
        falseN = 0;
        falseP = 0;
        errorR = 0;

        //Compare probe chars to find best delta sequence
        for (int a = 0; a < N; a++) {

          //Use substring to get the best 100-probe delta sequence
          for (unsigned int b = 0; b < sequence[a].length() - K; b++) {
            tempString = sequence[a].substr(b, K);

            //Compare strings character by characterx
            tempNum = 0;
            match = false;
            for (int c = 0; c < K; c++) {
              if(tempString[c] != result[c]) {
                tempNum++;
              }
              if (tempNum > 1) {
                break;
              }
            }
              if (tempNum <= 1) {
                match = true;
                break;
              }
            }
            if (is_delta[a] && match == false) {
              falseN++;
            }
            else if (!is_delta[a] && match == true) {
              falseP++;
            }
          }
          FPR = static_cast<double>(falseP) / not_delta;
          FNR = static_cast<double>(falseN) / delta;
          errorR = 2 * FPR + FNR;

          //Compute Best Probe Sequence
          if (tempER > errorR) {
            tempFN = falseN;
            tempFP = falseP;
            tempER = errorR;
            tempResult = result;
          }
        }
      }
      FPR = static_cast<double>(tempFP) / not_delta;
      FNR = static_cast<double>(tempFN) / delta;
      std::cout << "Best Probe Sequence: " << tempResult << '\n';
      printf("False Positives: %d False positive rate: %lf\n", tempFP, FPR);
      printf("False Negatives: %d False negative rate: %lf\n", tempFN, FNR);
      printf("The error rate: %lf\n", tempER);
      std::cout << "The worst-case running time of my approach: O(N^8 M^5 K)\n"; //Case in which no break statements would be added
      std::cout << "The anticipated running time of my approach: O(N^2 M^2 K)/X (X = Num of break statements)\n";
      //Deallocate array memory
      delete [] sequence;
      delete [] is_delta;

  return 0;
}
