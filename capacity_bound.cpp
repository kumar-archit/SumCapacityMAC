#include <iostream>
#include <bits/stdc++.h>
using namespace std;

int binToDec(vector<int>& rep){
  int temp = rep.size()-1;
  int siz = rep.size();
  int ans=0;
  for(; temp>=0;temp--){
    ans += rep[temp]*1<<(siz-temp-1);
  }
  return ans;
}
double lim_r1_sum_r2(vector<vector<double> >& transition_matrix, vector<double>& p1, vector<double>& p2) {
    double rate = 0.0;
    vector<int> rep;
    for(int y = 0; y < 3; ++y) {
            for(int x1 = 0; x1 < 2; ++x1) {
                for(int x2 = 0; x2 < 2; ++x2) {
                    double term = 1.0;
                    rep.clear();
                    rep.push_back(x1);
                    rep.push_back(x2);
                    // double term = 1.0;
                    term *= p1[x1];
                    term *= p2[x2];
                    int idx = binToDec(rep);
                    term *= transition_matrix[idx][y];
                    double num = transition_matrix[idx][y];
                    double deno = 0.0;
                    for(int k1 = 0; k1 < 2; ++k1) {
                        double contri = 1.0;
                        contri *= p1[k1];
                        for(int k2 = 0; k2 < 2; ++k2) { 
                            contri *= p2[k2];
                            rep.clear();
                            rep.push_back(k1);
                            rep.push_back(k2);
                            int index = binToDec(rep);
                            contri *= transition_matrix[index][y];
                            deno += contri;
                        }
                    }
                    term *= log(num/deno)/log(2);
                    rate += term;
                }
            }
        }
    return rate;
}
double calc_lim(vector<vector<double> >& transition_matrix, vector<double>& p1, vector<double>& p2, int opt) {
    double rate = 0.0;
    vector<int> rep;
    if(opt == 0) {
        for(int y = 0; y < 3; ++y) {
            for(int x1 = 0; x1 < 2; ++x1) {
                for(int x2 = 0; x2 < 2; ++x2) {
                    double term = 1.0;
                    rep.clear();
                    rep.push_back(x1);
                    rep.push_back(x2);
                    // double term = 1.0;
                    term *= p1[x1];
                    term *= p2[x2];
                    int idx = binToDec(rep);
                    term *= transition_matrix[idx][y];
                    double num = transition_matrix[idx][y];
                    double deno = 0.0;
                    for(int k = 0; k < 2; ++k) {
                        double contri = 1.0;
                        contri *= p1[k];
                        rep.clear();
                        rep.push_back(k);
                        rep.push_back(x2);
                        int index = binToDec(rep);
                        contri *= transition_matrix[index][y];
                        deno += contri;
                    }
                    term *= log(num/deno)/log(2);
                    rate += term;
                }
            }
        }
    }
    else {
        for(int y = 0; y < 3; ++y) {
            for(int x1 = 0; x1 < 2; ++x1) {
                for(int x2 = 0; x2 < 2; ++x2) {
                    double term = 1.0;
                    rep.clear();
                    rep.push_back(x1);
                    rep.push_back(x2);
                    // double term = 1.0;
                    term *= p1[x1];
                    term *= p2[x2];
                    int idx = binToDec(rep);
                    term *= transition_matrix[idx][y];
                    double num = transition_matrix[idx][y];
                    double deno = 0.0;
                    for(int k = 0; k < 2; ++k) {
                        double contri = 1.0;
                        contri *= p2[k];
                        rep.clear();
                        rep.push_back(x1);
                        rep.push_back(k);
                        int index = binToDec(rep);
                        contri *= transition_matrix[index][y];
                        deno += contri;     
                    }
                    term *= log(num/deno)/log(2);
                    rate += term;
                }
            }
        }
    }
    return rate;
}
int main() {
//	 vector<vector<double> > transition_matrix = {
//       {0.2, 0.3, 0.5},
//       {0.7, 0.2, 0.1}, 
//       {0.5, 0.1, 0.4}, 
//       {0.3, 0.4, 0.3}
//     };
     vector<vector<double> > transition_matrix = {
       {0.4, 0.1, 0.5}, 
       {0.3, 0.2, 0.5},
       {0.5, 0.4, 0.1}, 
       {0.2, 0.799, 0.001}
     };
//    vector<vector<double> > transition_matrix = {
//      {0.1, 0.2, 0.7},
//      {0.3, 0.5, 0.2},
//      {0.3, 0.4, 0.3},
//      {0.8, 0.1, 0.1}
//    };
    
	
    // for(int i = 0; i < transition_matrix.size(); ++i) {
    //     for(int j = 0; j < transition_matrix[i].size(); ++j)
    //         cout << transition_matrix[i][j] << " ";
    //     cout << endl;
    // }
    vector<double> p1, p2;
    int ctr = 0;
    for(double alpha = 1e-2; alpha <= 1.0 - 1e-2; alpha += 1e-2) {
        for(double beta = 1.0 - 1e-2; beta >= 1e-2; beta -= 1e-2) {
            p1.clear();
            p2.clear();
            // ++ctr;
            p1.push_back(alpha); p1.push_back(1 - alpha);
            p2.push_back(beta); p2.push_back(1-beta);
            double r1 = calc_lim(transition_matrix, p1, p2, 0);
            double r2 = calc_lim(transition_matrix, p1, p2, 1);
            double r1_add_r2 = lim_r1_sum_r2(transition_matrix, p1, p2);
            // cout << r1 << " " << r2 << " " << r1_add_r2 << endl;
 //            if(r1 + r2 <= r1_add_r2) {
 //               if(r1<0.01 && r2 < 0.01) r1=r2 =0;
 //               else if(r1 < 0.05 && r1 < r2) r1=0;
 //               else if(r2< 0.05 && r2 < r1) r2=0;
            //     // cout << ctr << endl;
 //                cout << r1 << " " << r2 << endl;
 //           }
 //           else{
                if(r1<0.01 && r2 < 0.01) r1=r2 =0;
                else if(r1 < 0.05 && r1 < r2) r1=0;
                else if(r2< 0.05 && r2 < r1) r2=0;
                cout << r1 << " " << r2 << endl;
            
//            }
        }
    }
    // double x = 1.0 - 1e-3;
    // cout << x << endl;
    return 0;
}
