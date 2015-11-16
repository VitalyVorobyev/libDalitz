#ifndef DRAWBDPARAMS_H
#define DRAWBDPARAMS_H

#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TCanvas.h"

#include <fstream>
#include <vector>
#include <sstream>
#include <string>

using namespace std;

class DrawBDParams{
public:
  DrawBDParams();
  void DrawBinsmABmAC(const string& infile,const string& outname){
    DrawBDP(infile,outname,1);
  }
  void DrawBinsmABmBC(const string& infile,const string& outname){
    DrawBDP(infile,outname,2);
  }
  void DrawBinsmACmBC(const string& infile,const string& outname){
    DrawBDP(infile,outname,3);
  }
  void DrawBDP(const string& infile,const string& outname,const int type);
  void DrawCS(const vector<double>& C,const vector<double>& S,const string& fname);
  void DrawK(const vector<double>& K,const vector<double>& Kb,const string& fname);

  void SetCSRef(vector<double>& C,vector<double>& S) {Cref = &C; Sref = &S; m_csrf = true; return;}
  void SetKRef(vector<double>& K,vector<double>& Kb) {Kref = &K; Kbref = &Kb; m_krf = true; return;}
  void RemoveCSRef(void) {m_csrf = false; return;}
  void RemoveKRef(void) {m_krf = false; return;}

private:
//  int nbins;
  vector<double>* Kref;
  vector<double>* Kbref;
  vector<double>* Cref;
  vector<double>* Sref;
  bool m_krf;
  bool m_csrf;
};

#endif // DRAWBDPARAMS_H
