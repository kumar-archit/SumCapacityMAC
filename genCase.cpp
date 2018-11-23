#include <bits/stdc++.h>
using namespace std;
//Initialization
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
//Utility function
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
int binToDec(vector<int> rep){
  int temp = rep.size()-1;
  int siz = rep.size();
  int ans=0;
  for(; temp>=0;temp--){
    ans += rep[temp]*1<<(siz-temp-1);
  }
  return ans;
}
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

void maximisation(vector<vector<double> > &transmissionMatrix, vector<vector<double> > &inputAlphabetsDistribution, vector<vector<double> > &mutualIm){
	prob_adjustment(transmissionMatrix, inputAlphabetsDistribution, mutualIm);	
}

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
    marginalisation(mutualIm, outputDis, inputAlphabetsDistribution, transmissionMatrix);
    maximisation(transmissionMatrix, inputAlphabetsDistribution, mutualIm);
  }
  return capa;
}
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
