/* A program to find the capacity of MAC channel Single user case by implementing Blahut Arimoto Algorithm
 *
 * implement.cpp
 *
 *  Created on: 28-Aug-2018
 *      Author: Rajanish Upadhyay
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

#include <vector>
#include <iostream>
#include <bits/stdc++.h>
using namespace std;

//Input r and p from user/ file
void init_user(vector<double>& r, vector<vector<double> >& p) {
	for(int i = 0; i < r.size(); ++i)
		cin >> r[i];
	for(int i = 0; i < p.size(); ++i)
		for(int j = 0; j < p[i].size(); ++j)
			cin >> p[i][j];
}

//calculate q using r and p
void q_(vector < double > &r,vector < vector < double > > &p,vector < vector < double > > &q ){
// 	 cout<<"Calculating q"<<endl;
	int i, j, k;
	double c=0;
	for(i=0; i<r.size(); i++){
		for(j=0; j< p[i].size(); j++){
			c=0;
			for(k=0; k< r.size(); k++){
				c+= r[k]*p[k][j];
			}
			q[j][i] = r[i]*p[i][j]/c;
		}
	}
// 	cout<<"Updated value of q matrix : "<<endl;
// 	for(i=0; i<r.size(); i++){
// 		for(j=0; j<p[i].size(); j++){
// 			cout << q[i][j] << " ";
// 		}
// 		cout << endl;
// 	}

// 	 cout<<"Done c q"<<endl;
}

//update r using p and q
void r_(vector <double> &r,vector<vector<double> > &p, vector<vector<double> > &q){
// 	 cout<<"Calculating r"<<endl;
	double sum =0, prod =1;
	int i, j;
	for(i=0; i<r.size(); i++){
		prod = 1;
		for(j=0; j< p[i].size(); j++){
			prod *= pow(q[j][i], p[i][j]);
		}
		r[i] = prod;
		sum += prod;
	}
	for(i=0; i<r.size(); i++){
		r[i] /= sum;
	}
// 	cout<<"Updated value of r : "<<endl;
// 	for(i=0; i< r.size(); i++){
// 		cout<<" "<<r[i];
// 	}
// 	cout<<endl;
}

// function to calculate capacity
double c(vector<double> &r,vector<vector<double> > &p,std::vector<std::vector<double> > &q ){
// 	 cout<<"Calculating c"<<endl;
	double capacity=0;
	int i,j;
	for(i=0; i< r.size(); i++){
		for(j=0; j< p[i].size(); j++){
			if(p[i][j] == 0)
				continue;
			capacity += r[i]*p[i][j]*log2(q[j][i]/r[i]);
		}
	}
// 	cout<<"Updated value of c ; "<<capacity<<endl;
	return capacity;
}

//Expectation-Maximisation
double em(std::vector<double> &r, std::vector<std::vector<double> > &p, std::vector< vector <double> > &q, int iters = 100, double error  = 0.01){
	//cout<<"Expectation Maximisation"<<endl;
	double oldc,capacity=0.0 , diff = 1;
	while(iters-- && diff >= error){
		cout<<"In iteration : "<<iters<<endl;
		q_(r,p,q);
		r_(r,p,q);
		oldc = capacity;
		cout << "oldc " << oldc << endl;
		capacity = c(r,p,q);
		diff = abs(oldc - capacity);
	}
	//cout<<"Expectation Maximisation complete"<<endl;
	return capacity;
}
/* I_ALP_SIZE	:	input alphabet size
 * O_ALP_SIZE	:	output alphabet size
 * r		:	a vector to store the input probability distribution 
 * p		:	transition matrix for the channel
 * q		:	reverse transition matrix
 * init_user()	:	a function which initialises r, p & q
 * em()		:	a function implementing EM method to calculate capacity 
 */
int main() {
	int I_ALP_SIZE, O_ALP_SIZE;
	cin >> I_ALP_SIZE >> O_ALP_SIZE;
	vector<double> r(I_ALP_SIZE);						//input distribution
	vector<vector<double> > p(I_ALP_SIZE, std::vector<double> (O_ALP_SIZE,0)); 			//transition matrix p(y/x)
	vector<vector<double> > q(O_ALP_SIZE, std::vector<double> (I_ALP_SIZE,0));
	init_user(r,p);			//reverse transition matrix q(x/y)
	double capacity;
	capacity  = em(r,p,q, 10);
	cout<<capacity<<endl;
}



