#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include <vector>
#include "TMath.h"

std::pair<double,double> findIntersectionXV(double x,double v);
std::pair<double,double> findIntersectionXU(double x,double u);
std::pair<double,double> findIntersectionUV(double u,double v);

void plotter(){

  TFile *f = new TFile("../tracks.root");
  TTree *t = (TTree*) f->Get("tracks");
  Int_t nEntries = t->GetEntries();

  TH2F *h_first_hits = new TH2F("h_first_hits","",332,0,332,332,0,332);
  TH2F *h_last_hits = new TH2F("h_last_hits","",332,0,332,332,0,332);
  TH1F *h_track_lenght = new TH1F("h_track_length","",100,0,100);
  TH1F *h_track_angle = new TH1F("h_track_angle","",100,-3.15,3.15);
 
  cout<<nEntries<<endl;
  
  std::vector<double> *Xhit = 0;
  std::vector<double> *Vhit = 0;
  std::vector<double> *Uhit = 0;

  std::vector<double> *Xhit2 = 0;
  std::vector<double> *Vhit2 = 0;
  std::vector<double> *Uhit2 = 0;

  std::vector<bool> *beam = 0;
  
  t->SetBranchAddress("first_hit_X", &Xhit);
  t->SetBranchAddress("first_hit_V", &Vhit);
  t->SetBranchAddress("first_hit_U", &Uhit);

  t->SetBranchAddress("last_hit_X", &Xhit2);
  t->SetBranchAddress("last_hit_V", &Vhit2);
  t->SetBranchAddress("last_hit_U", &Uhit2);

  t->SetBranchAddress("beam", &beam);
  
  for (int i=0; i< nEntries; i++) {
    t->GetEntry(i);
    cout<<"Event="<<i+1<<endl;
    for (int j = 0; j<Xhit->size(); j++) {

      if ((beam->at(j))) continue;

      //cout<<"XVU="<<Xhit->at(j)<<" "<<Vhit->at(j)<<" "<<Uhit->at(j)<<endl;
      
      std::pair<double,double> xv(findIntersectionXV(Xhit->at(j),Vhit->at(j)));
      std::pair<double,double> xu(findIntersectionXU(Xhit->at(j),Uhit->at(j)));
      std::pair<double,double> uv(findIntersectionUV(Uhit->at(j),Vhit->at(j)));
      
      std::pair<double,double> xv2(findIntersectionXV(Xhit2->at(j),Vhit2->at(j)));
      std::pair<double,double> xu2(findIntersectionXU(Xhit2->at(j),Uhit2->at(j)));
      std::pair<double,double> uv2(findIntersectionUV(Uhit2->at(j),Vhit2->at(j)));
      
      // cout<<xv.first<<" "<< xv.second<<endl;
      // cout<<xu.first<<" "<< xu.second<<endl;
      // cout<<uv.first<<" "<< uv.second<<endl;

      double first_hit_x = (xv.first+xu.first+uv.first)/3.0;
      double first_hit_y = (xv.second+xu.second+uv.second)/3.0;
      
      double last_hit_x = (xv2.first+xu2.first+uv2.first)/3.0;
      double last_hit_y = (xv2.second+xu2.second+uv2.second)/3.0;

      double track_length = TMath::Sqrt((last_hit_x-first_hit_x)*(last_hit_x-first_hit_x)+(last_hit_y-first_hit_y)*(last_hit_y-first_hit_y));
      double track_angle = 0;
      track_angle = TMath::ATan((last_hit_y-first_hit_y)/(last_hit_x-first_hit_x));
	
      h_first_hits->Fill(first_hit_x,first_hit_y);
      h_last_hits->Fill(last_hit_x,last_hit_y);      
      h_track_lenght->Fill(track_length*0.3);
      h_track_angle->Fill(track_angle);

    }
  }

  TCanvas *c1 = new TCanvas("c1","",900,600);
  h_track_length->GetXaxis()->SetTitle("Track length [cm]");
  h_track_length->Draw();
  c1->Print("plots/tlength_cosmic.png");
  c1->Print("plots/tlength_cosmic.pdf");
  TCanvas *c2 = new TCanvas("c2","",900,600);
  h_track_angle->Draw();
  c2->Print("plots/tangle_cosmic.png");
  c2->Print("plots/tangle_cosmic.pdf");
  //h_first_hits->Draw("colz");
  
}

std::pair<double,double> findIntersectionXV(double x,double v) {
  TF1 *line = new TF1("line","0.5/(cos((pi)/6.0))*x+[0]",0,335);
  line->SetParameter(0,-94.84+v);
  return std::make_pair(x,line->Eval(x));
}

std::pair<double,double> findIntersectionXU(double x,double u) {
  TF1 *line = new TF1("line","-0.5/(cos((pi)/6.0))*x+[0]",0,335);
  line->SetParameter(0,425.84-u);
  return std::make_pair(x,line->Eval(x));
}

std::pair<double,double> findIntersectionUV(double u,double v) {
  TF1 *line1 = new TF1("line1","-0.5/(cos((pi)/6.0))*x+[0]",0,335);
  line1->SetParameter(0,425.84-u);
  TF1 *line2 = new TF1("line2","0.5/(cos((pi)/6.0))*x+[0]",0,335);
  line2->SetParameter(0,-94.84+v);
  double x = ((425.84-u+94.84-v))/(1/(TMath::Cos((TMath::Pi())/6.0)));
  
  return std::make_pair(x,line1->Eval(x));
}
