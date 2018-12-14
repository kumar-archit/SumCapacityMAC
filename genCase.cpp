/* A program to find the capacity of MAC channel Multi user case (each alphabet is binary) by implementing Blahut Arimoto Algorithm
 *
 *  Created on: 10-Oct-2018
 *     
 */

/*
Input[Optional]: 
1. input probability distribution of input alphabet. 
2. transition matrix for the channel
[3]. Number of iterations for Expectation-Maximisation to run.
[4]. Accepted error in Capacity.
Output:
Capacity of the Channel
*/

#include <bits/stdc++.h>
using namespace std;
//Initialization of inputAlphabetsDistribution vector and transmissionMatrix 2d vector by taking input from user/ file
void init(vector<vector<double> > &inputAlphabetsDistribution, vector<vector<double> > &transmissionMatrix){
  int siz=inputAlphabetsDistribution.size(), i,j;
  for(i=0; i<siz; i++){
    for(j=0; j<2; j++)
      cin>>inputAlphabetsDistribution[i][j];
  }
  siz=transmissionMatrix.size();
  for(i=0; i<siz; i++){
    for(j=0; j<transmissionMatrix[0].size(); j++)
      cin>>transmissionMatrix[i][j];
  }
}

/*
 * a function to convert a decimal number to binary number, the size of binary number is equal to inputAlphabetsDisSize.
 * this function is required to find the row corresponding to the input vector in the transmission matrix. Here decToBin  
 * is useful only because the input alphabets are binary. In general case xM method is used. See also: binToDec
 */
void decToBin(int x, int inputAlphabetsDisSize, vector<int>& rep) {
  for(int i = 0; i < inputAlphabetsDisSize; ++i) rep.push_back(0);
  int i = 0;
  while(x > 0) {
    rep[i] = x % 2;
    ++i;
    x = x/2;
  }
  reverse(rep.begin(), rep.end());
}

/*
 * a function to calculate capacity  by first calculating probability distribution of output alphabet
 * noOfPermutation	:	total number of possible input vectors 
 * probY		:	vector to store output alphabet probability distribution
 */
double capacity(vector<double> &outputDis, vector<vector<double> > &inputAlphabetsDistribution, vector<vector<double> > &transmissionMatrix){
  //calc P(Y)
  int noOfPermutation = transmissionMatrix.size();
  int ysize = transmissionMatrix[0].size();
  int inpAlpSize = inputAlphabetsDistribution.size();
  double cumPro, probY;
  for(int y = 0; y < ysize; ++y) {
		probY = 0;
    for(int i = 0; i < noOfPermutation; ++i) {
      vector<int> rep;
      decToBin(i, inpAlpSize, rep);  
      cumPro = 1;
      for(int j = 0; j < rep.size(); ++j) {
        cumPro *= inputAlphabetsDistribution[j][rep[j]];
      }
      probY += cumPro*transmissionMatrix[i][y];
   }
    outputDis[y] = probY;
  }
  double capacity=0;
  for(int y=0; y<ysize; y++){
    for(int i=0; i< noOfPermutation; ++i){
      vector<int> rep;
      decToBin(i, inpAlpSize, rep);
      cumPro = 1;
      for(int j = 0; j < rep.size(); ++j) {
        cumPro *= inputAlphabetsDistribution[j][rep[j]];
      }
      capacity += cumPro*transmissionMatrix[i][y]*log(transmissionMatrix[i][y]/outputDis[y]);
    }
  }

  return capacity;
}
/*
 * a function to change the probability distributions of input alphabets using 
 * mutual information. This updation of probability distribution happens in each step of EM method.
 * see also: prob_adjustment_updated
 */
void prob_adjustment(vector<vector<double> > &transitionMatrix, vector<vector<double> >& inputAlphabetsDistribution, vector<vector<double> > &mutualIm) {
  int M  = inputAlphabetsDistribution.size();
  for(int i=0; i<M; i++){
    for(int j=0; j<2; j++){
      double num; 
      double deno=0;
      num = mutualIm[i][j];
      for(int k=0; k < 2 ; k++) {
		  deno += (inputAlphabetsDistribution[i][k]*mutualIm[i][k]);
	  }
      if(deno == 0) cout << "Error: Denominator zero" << endl;
	  inputAlphabetsDistribution[i][j] *= num/deno;
    }
  }
}

/*
 * similar to prob_adjustment, implementation slightly different due to the fact of using a different function monotonically  
 * increasing function than exponential.
 */
void prob_adjustment_updated(vector<vector<double> > &transitionMatrix, vector<vector<double> >& inputAlphabetsDistribution, vector<vector<double> > &mutualIm) {
  int M  = inputAlphabetsDistribution.size();
  for(int i=0; i<M; i++){
    for(int j=0; j<2; j++){
      double num; 
      double deno=0;
	  // cout << "vm: " << verti_marginalVals[i][j] << endl;
      num = exp(mutualIm[i][j]) - 1;
      for(int k=0; k < 2 ; k++) {
		  // cout << "inp Dis: " << inputAlphabetsDistribution[i][k] << endl;
		  deno += (inputAlphabetsDistribution[i][k])*(exp(mutualIm[i][k]) - 1);
	  }
      if(deno == 0) cout << "Chutiya: " << endl;
	  // cout << "Num : " << num << " " << "Deno : " << deno << endl;
	  inputAlphabetsDistribution[i][j] *= num/deno;
    }
  }
}

/*
 * converts a binary number to decimal number
 */
int binToDec(vector<int> rep){
  int temp = rep.size()-1;
  int siz = rep.size();
  int ans=0;
  for(; temp>=0;temp--){
    ans += rep[temp]*1<<(siz-temp-1);
  }
  return ans;
}


/*
 * a function to calculate mutual information of each alphabets' every symbol in a 2d vector mutualIm
 * tempDistribution is a 2d vector used to store the probability distribution of all alphabets other than the 
 * currently focussed one in Marginalisation algorithm. This matrix is of size (noOfInputAlphabets-1)*2.
 *
 * a temporary vector xm_dash (one less dimension) is generated which can be used to in-turn generate many more  
 * (equal to size of omitted alphabet).
 *
 */
void marginalisation( vector<vector<double> > &mutualIm, vector<double> &outputDis, vector<vector<double> > &inputAlphabetsDistribution, vector<vector<double> > &transmissionMatrix){
  int M = inputAlphabetsDistribution.size();
  int count;
  double product;
  vector<vector<double> > tempDistribution(M-1, vector<double> (2,0));
	int lim = pow(2, M-1);
	for(int i = 0; i < M; i++) {
      count=0;
      for(int j = 0; j < M; j++) if(j != i){
        for(int k=0; k<2; k++){
          tempDistribution[count][k] = inputAlphabetsDistribution[j][k];
        }
        count++;
      } 
      for(int k=0; k<2; k++){  
        product = 1;
        for(int j = 0; j < lim; ++j) {
          vector<int> xm_dash;
          decToBin(j, M - 1, xm_dash);
          vector<int> xm = xm_dash;
          xm.emplace(xm.begin() + i, k);
          int pos = binToDec(xm);
          double p_xm_dash=1 ;
          for(int p=0; p<M-1; p++){
            p_xm_dash *= tempDistribution[p][xm_dash[p]];
          }
          for(int y = 0; y < outputDis.size(); ++y) {
            double p_y_given_x_M = transmissionMatrix[pos][y];
            double p_y = outputDis[y];
            product *= pow((p_y_given_x_M)/(p_y), p_xm_dash*p_y_given_x_M);
          }
        }
        mutualIm[i][k] = product;
      }
    }
}

/*
 * marginalisation_updated is a similar function to marginalisation with slightly different implementation due to
 * use of a different monotonically increasing function than exponential function.
 *
 */
void marginalisation_updated( vector<vector<double> > &mutualIm, vector<double> &outputDis, vector<vector<double> > &inputAlphabetsDistribution, vector<vector<double> > &transmissionMatrix){
  int M = inputAlphabetsDistribution.size();
  int count;
  double product;
  vector<vector<double> > tempDistribution(M-1, vector<double> (2,0));
	int lim = pow(2, M-1);
	for(int i = 0; i < M; i++) {
      count=0;
      for(int j = 0; j < M; j++) if(j != i){
        for(int k=0; k<2; k++){
          tempDistribution[count][k] = inputAlphabetsDistribution[j][k];
        }
          count++;
      } 
      for(int k=0; k<2; k++){  
	      double sum = 0;
      	 for(int j = 0; j < lim; ++j) {
           vector<int> xm_dash;
           decToBin(j, M - 1, xm_dash);
           vector<int> xm = xm_dash;
           xm.emplace(xm.begin() + i, k);
           int pos = binToDec(xm);
           double p_xm_dash=1 ;
           for(int p=0; p<M-1; p++){
             p_xm_dash *= tempDistribution[p][xm_dash[p]];
           }
           double p_xm = 1;
           for(int p=0; p<M; p++){
             p_xm *= inputAlphabetsDistribution[p][xm[p]];
           }
           double term = 0;
           //double contri = 1;
           for(int y = 0; y < outputDis.size(); ++y) {
             double p_y_given_x_M = transmissionMatrix[pos][y];
             double p_y = outputDis[y];
             term += p_y_given_x_M*log(p_y_given_x_M/p_y);

            //  product *= pow((p_y_given_x_M)/(p_y), p_xm_dash*p_y_given_x_M);
           }
		   term *= p_xm_dash;
		   sum += term;
         }
        mutualIm[i][k] = sum;
       }
  		
     }
}

/*
 *
 *
 */
void maximisation(vector<vector<double> > &transmissionMatrix, vector<vector<double> > &inputAlphabetsDistribution, vector<vector<double> > &mutualIm){
// 	prob_adjustment(transmissionMatrix, inputAlphabetsDistribution, mutualIm);	
	prob_adjustment_updated(transmissionMatrix, inputAlphabetsDistribution, mutualIm);	

}
/*
 * Implementation of Expectation-Maximisation paradigm
 * Expectation step is performed by marginalisation function
 * Maximisation step is performed by Maximisation function
 */
double em(vector<vector<double> > &inputAlphabetsDistribution, vector<vector<double> > &transmissionMatrix, double error = 0.01, int iterations = 1){
  double capa, old_c = 0;
  int M = inputAlphabetsDistribution.size();
  vector<double> outputDis(transmissionMatrix[0].size(),0);
  vector<vector <double> >mutualIm(M, vector<double> (2,0));
  int LIM = iterations; 
  while(iterations--){
    capa = capacity(outputDis, inputAlphabetsDistribution, transmissionMatrix);
    if(abs(capa - old_c) < error) return capa;
    old_c = capa;
	  cout << LIM - iterations << " " << capa << endl;
//     marginalisation(mutualIm, outputDis, inputAlphabetsDistribution, transmissionMatrix);
    marginalisation_updated(mutualIm, outputDis, inputAlphabetsDistribution, transmissionMatrix);
    maximisation(transmissionMatrix, inputAlphabetsDistribution, mutualIm);
  }
  return capa;
}

/* noOfInputAlphabets:	total input users
 * outputAlphabetSize:	no of symbols in output alphabet
 * inputAlphabetDistribution:	a 2d vector to store probability distributions of all the input alphabets
 * tranmissionMatrix: a matrix of size (product of size of all input alphabets) * size of output alphabet to store
 *
 */
int main(){
    int noOfInputAlphabets, outputAlphabetSize, siz;
    cin>>noOfInputAlphabets >> outputAlphabetSize;
    siz = 1 << noOfInputAlphabets;
    vector <vector<double> > inputAlphabetsDistribution(noOfInputAlphabets, vector<double> (2, 0));
    vector <vector<double> > transmissionMatrix(siz, vector<double> (outputAlphabetSize, 0));
    init(inputAlphabetsDistribution, transmissionMatrix);
	
    em(inputAlphabetsDistribution, transmissionMatrix, 0.000001, 10000);

	return 0;
}
