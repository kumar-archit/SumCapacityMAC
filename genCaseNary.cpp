//general case arbitrary alphabet size
//MAC
#include <bits/stdc++.h>
using namespace std;
void xM(int x, vector<vector<double> > &inputAlphabetsDistribution, vector<int>& inpVector) {
  int i, siz;
  for(i=0; i<inpVector.size(); i++)inpVector[i]=0;
	i=0;
  while(x > 0) {
  	siz = (inputAlphabetsDistribution[i].size());
    inpVector[i] = x % siz;
    ++i;
    x /= siz;
  }
  reverse(inpVector.begin(), inpVector.end());
}
int nPerm(vector<int> inpVector, vector<vector<double> > &inputAlphabetsDistribution){
  int siz = inpVector.size();
  int ans=0;
  for(int temp = siz - 1; temp>=0;temp--){
  	//here we mush denote ordinal values to symbols i.e. they must start from 0 and be sequential integers
    ans += inpVector[temp]*(inputAlphabetsDistribution[temp].size());
  }
  return ans;
}

double capacity(vector<double> &outputDis, vector<vector<double> > &inputAlphabetsDistribution, vector<vector<double> > &transmissionMatrix){
  //calc P(Y)
	int noOfPermutation = transmissionMatrix.size();
	int outputAlphabetSize = transmissionMatrix[0].size();
	int noOfInputAlphabets = inputAlphabetsDistribution.size();
	double cumPro, probY, capacity=0;
	vector<int> inpVector(noOfInputAlphabets,0);
	int i, j, y;
	for(y = 0; y < outputAlphabetSize; ++y) {
		probY = 0;
		for(i = 0; i < noOfPermutation; ++i) {
		  xM(i, inputAlphabetsDistribution, inpVector);  
		  cumPro = 1;
		  for(j = 0; j < inpVector.size(); ++j) {
		    cumPro *= inputAlphabetsDistribution[j][inpVector[j]];
		  }
		  probY += cumPro*transmissionMatrix[i][y];
		}
		outputDis[y] = probY;
	}

	for(y=0; y<outputAlphabetSize; y++){
		for(i=0; i< noOfPermutation; ++i){
		  xM(i, noOfInputAlphabets, inpVector);
		  cumPro = 1;
		  for(j = 0; j < inpVector.size(); ++j) {
		    cumPro *= inputAlphabetsDistribution[j][inpVector[j]];
		  }
		  capacity += cumPro*transmissionMatrix[i][y]*log(transmissionMatrix[i][y]/outputDis[y]);
		}
	}

	return capacity;
}
void marginalisation( vector<vector<double> > &mutualIm, vector<double> &outputDis, vector<vector<double> > &inputAlphabetsDistribution, vector<vector<double> > &transmissionMatrix){
	int noOfInputAlphabets = inputAlphabetsDistribution.size();
	int count, siz, pos;
	double product, p_xm_dash;
	vector<vector<double> > tempDistribution;
	tempDistribution.resize(noOfInputAlphabets-1);
	vector<int> xm_dash(noOfInputAlphabets-1);
	vector<int> xm;
	// int lim = pow(2, M-1);
	int noOfPermutation = 1, lim;
	for(int i=0; i<noOfInputAlphabets; i++) 
		noOfPermutation *= inputAlphabetsDistribution.size();

	for(int i = 0; i < noOfInputAlphabets; i++) {
		count=0;
		for(int j = 0; j < noOfInputAlphabets; j++) if(j != i){
			siz = inputAlphabetsDistribution[j].size();
			tempDistribution[count].resize(siz);
			for(int k=0; k<siz; k++){
				tempDistribution[count][k] = inputAlphabetsDistribution[j][k];
			}
			count++;
		}
		siz = inputAlphabetsDistribution[i].size();
		lim = noOfPermutation/siz; 
		for(int k=0; k<siz; k++){  
			product = 1;
			 for(int j = 0; j < lim; ++j) {
				xM(j, tempDistribution, xm_dash);
				xm = xm_dash;
				xm.emplace(xm.begin() + i, k);
				pos = nPerm(xm);
				p_xm_dash=1 ;
				for(int p=0; p<noOfInputAlphabets-1; p++){
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
void prob_adjustment(vector<vector<double> >& inputAlphabetsDistribution, vector<vector<double> > &mutualIm) {
	int noOfInputAlphabets  = inputAlphabetsDistribution.size();
	double num, deno;
	for(int i=0; i<noOfInputAlphabets; i++){
		for(int j=0; j<inputAlphabetsDistribution[i].size(); j++){
			deno=0;
			num = mutualIm[i][j];
			for(int k=0; k < inputAlphabetsDistribution[i].size() ; k++) {
				deno += (inputAlphabetsDistribution[i][k]*mutualIm[i][k]);
			}
			if(deno == 0) 
				cout << "Error: Denominator zero" << endl;
			inputAlphabetsDistribution[i][j] *= num/deno;
		}
	}
}

double em(vector<vector<double> > &inputAlphabetsDistribution, vector<vector<double> > &transmissionMatrix, double error = 0.01, int iterations = 1){
	double capa, old_c = 0;
	int noOfInputAlphabets = inputAlphabetsDistribution.size(), outputAlphabetSize=transmissionMatrix[0].size();
	vector<double> outputDis(outputAlphabetSize,0);
	vector<vector <double> >mutualIm;
	mutualIm.resize(noOfInputAlphabets);
	for(int i=0;i<noOfInputAlphabets; i++){
		mutualIm[i].resize(inputAlphabetsDistribution[i].size(),0);
	}
	int LIM = iterations; 
	while(iterations--){
		capa = capacity(outputDis, inputAlphabetsDistribution, transmissionMatrix);
		if(abs(capa - old_c) < error) return capa;
		old_c = capa;
		cout << LIM - iterations << " " << capa << endl;
		marginalisation(mutualIm, outputDis, inputAlphabetsDistribution, transmissionMatrix);
		prob_adjustment(inputAlphabetsDistribution, mutualIm) {

	}
}

void init(int noOfInputAlphabets, int outputAlphabetSize, vector<int> alphabetSize, vector<vector<double> > &inputAlphabetsDistribution, vector<vector<double> > &transmissionMatrix){
	int  i,j, siz=1, tSiz;
	alphabetSize.resize(noOfInputAlphabets);
	for(i=0; i<noOfInputAlphabets; i++){
		cin>>tSiz;
		alphabetSize[i]=tSiz;
		siz *= tSiz;
		inputAlphabetsDistribution[i].resize(tSiz);
		for(j=0; j<tSiz; j++)
			cin>>inputAlphabetsDistribution[i][j];
	}
	transmissionMatrix.resize(siz);
	for(i=0; i<siz; i++)
		transmissionMatrix[i].resize(outputAlphabetSize);

	for(i=0; i<siz; i++){
		for(j=0; j<outputAlphabetSize; j++)
		  cin>>transmissionMatrix[i][j];
	}
}

int main(){
    int noOfInputAlphabets, outputAlphabetSize;
    cin>>noOfInputAlphabets >> outputAlphabetSize;
    //the transmission matrix is still of size 2^(Sum n_i)*outputAlphabetSize
	vector <vector<double> > inputAlphabetsDistribution;
	vector<int> alphabetSize;
    vector <vector<double> > transmissionMatrix;
    init(noOfInputAlphabets, outputAlphabetSize, alphabetSize, inputAlphabetsDistribution, transmissionMatrix);
    em(inputAlphabetsDistribution, transmissionMatrix, 0.000001, 10000);
	return 0;
}
