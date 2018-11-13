//binary alphabet multi-user
#include <bits/stdc++.h>
using namespace std;
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
vector<int> decToBin(int x, int inputAlphabetsDisSize) {
  vector<int> ans(inputAlphabetsDisSize, 0);
  int i = 0;
  while(x > 0) {
    ans[i] = x % 2;
    ++i;
    x = x/2;
  }
  reverse(ans.begin(), ans.end());
  return ans;
}

double capacity(vector<vector<double> > &inputAlphabetsDistribution, vector<vector<double> > &transmissionMatrix){
  //calc P(Y)
  int noOfPermutation = transmissionMatrix.size();
  int ysize = transmissionMatrix[0].size();
  double cumPro, probY;
  vector<double> outputDis(ysize, 0);
  for(int y = 0; y < ysize; ++y) {
		probY = 0;
    for(i = 0; i < noOfPermutation; ++i) {
      vector<int> rep = decToBin(i, noOfPermutation);
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
      vector<int> rep = decToBin(i, noOfPermutation);
      cumPro = 1;
      for(int j = 0; j < rep.size(); ++j) {
        cumPro *= inputAlphabetsDistribution[j][rep[j]];
      }
      capacity += cumPro*transmissionMatrix[i][y]*log10(transmissionMatrix[i][y]/outputDis[y])/log10(2);
    }
  }
  
  return capacity;
}
void expectation(vector<vector<double> > &mutualIm, vector<vector<double> > &inputAlphabetsDistribution, double capa){
  //mutualIm: no. of rows = no. of users
  //No. of cols = Size of input alphabet: Here 2
  double cumPro, cumSum;
  for(int m = 0; m < mutualIm.size(); m++){
    	cumSum = 0;
    	for(int i = 0; i < inputAlphabetsDistribution[0].size(); ++i) {
        cumSum += inputAlphabetsDistribution[m][i];
      }
      for(int i = 0; i < inputAlphabetsDistribution[0].size(); ++i) {
        
      }
  }
}
void maximisation(vector<vector<double> > &mutualIm, vector<vector<double> > &inputAlphabetsDistribution, double capa){
  double cumPro, cumSum = 0;
  //Calculate Normalization factor: cumsum
  for(int i = 0; i < mutualIm.size(); i++){
      cumPro = 1;
      vector<int> rep = decToBin(i, noOfPermutation);
      for(int j = 0; j < rep.size(); ++j) {
        cumPro *= inputAlphabetsDistribution[j][rep[j]];
      	}
      cumSum += cumPro*exp(mutualIm[i]);
  }
  
  for(int i=0; i<mutualIm.size(); i++){
    
  }
}
double em(vector<vector<double> > &inputAlphabetsDistribution, vector<vector<double> > &transmissionMatrix, double error = 0.01, int iterations = 1){
  double capa, old_c = 0;
	vector<vector<double> > mutualIm(inputAlphabetsDistribution.size(),vector<double> (2,0));
  while(iterations--){
	  capa = capacity(inputAlphabetsDistribution, transmissionMatrix);
  	if(abs(capa - old_c) < error) return capa; 
    old_c = capa;
    expectation(mutualIm, inputAlphabetsDistribution, capa);
    maximisation(mutualIm, inputAlphabetsDistribution, capa);
  }
  return capa;
}
int main(){
  int noOfInputAlphabets, outputAlphabetSize,	 siz;
  cin>>noOfInputAlphabets;
  siz = 1 << noOfInputAlphabets;
  vector <vector<double> > inputAlphabetsDistribution(noOfInputAlphabets, vector<double> (2, 0));
  vector <vector<double> > transmissionMatrix(siz, vector<double> (outputAlphabetSize, 0));
  init(inputAlphabetsDistribution, transmissionMatrix);
	em(inputAlphabetsDistribution, transmissionMatrix, 0.01, 1);
	return 0;
}