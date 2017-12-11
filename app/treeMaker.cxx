#include <memory>

// Include to get Event Loop.
#include <eventLoop.hxx>
#include <TReconTrack.hxx>
#include <TReconHit.hxx>
#include <TRealDatum.hxx>

#include <TEvent.hxx>

// Includes for ROOT classes

#include <HEPUnits.hxx>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <TH1F.h>
#include <TPad.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>



std::string toString(int i)
{
    std::ostringstream s;
    s << i;
    return s.str();
}

class TTreeMakerLoop: public CP::TEventLoopFunction {
public:
    TTreeMakerLoop() {

	first_hit_X.clear();
	last_hit_X.clear();

	first_hit_Y.clear();
	last_hit_Y.clear();
	
	first_hit_Z.clear();
	last_hit_Z.clear();

	hfile= new TFile("tracks.root","RECREATE");
	tree = new TTree("tracks","");
	tree->Branch("first_hit_X",&first_hit_X);
	tree->Branch("last_hit_X",&last_hit_X);
	tree->Branch("first_hit_Y",&first_hit_Y);
	tree->Branch("last_hit_Y",&last_hit_Y);
	tree->Branch("first_hit_Z",&first_hit_Z);
	tree->Branch("last_hit_Z",&last_hit_Z);
    }
    virtual ~TTreeMakerLoop() {}
    void Initialize() {
	
	
    }
    bool operator () (CP::TEvent& event) {
	
	std::cout<<event.GetContext()<<std::endl;

	CP::THandle<CP::TReconObjectContainer> tracks = event.Get<CP::TReconObjectContainer>("~/fits/TCaptainRecon/final");
	CP::THandle<CP::TDataVector> dataPMT = event.Get<CP::TDataVector>("~/pmtData");

	TString pdsEvent = "";
	for (int i=0; i<dataPMT->size(); i++) {
	    pdsEvent.Form("~/pmtData/PDSEvent_%d",i);
	    CP::THandle<CP::TEvent> eventPMT = event.Get<CP::TEvent>(pdsEvent);
	    if (!eventPMT) {
		std::cout<<"NO PMT EVENT"<<std::endl;
	    }
	    else {
		double TOF = (eventPMT->Get<CP::TRealDatum>("TOF_ns"))->GetValue();
		if (TOF > 0)
		    std::cout<<"TOF="<<(TOF)<<std::endl;
	    }
	}
	
	if (tracks) {
	    //std::cout<<"TRACKS"<<std::endl;
	    for (CP::TReconObjectContainer::const_iterator t = tracks->begin(); t != tracks->end(); ++t) {
		CP::THandle<CP::THitSelection> hits = (*t)->GetHits();
		//std::cout<<"TRACK"<<std::endl;
		double min_z = 9999;
		double max_z = -9999;
		TVector3 min_hit;
		TVector3 max_hit;		 
		for (CP::THitSelection::const_iterator h = hits->begin(); h != hits->end(); ++h) {
		    TVector3 v = (*h)->GetPosition();
		    //std::cout<<"HIT X="<<v.X()<<" HIT Y="<<v.Y()<<" HIT Z="<<v.Z()<<std::endl;
		    if (v.Z() > max_z) {
			max_z = v.Z();
			max_hit = v;
		    }		    
		    if (v.Z() < min_z) {
			min_z = v.Z();
			min_hit = v;
		    }
		}
		//std::cout<<"MIN HIT X="<<min_hit.X()<<" Y="<<min_hit.Y()<<" Z="<<min_hit.Z()<<std::endl;		 
		//std::cout<<"MAX HIT X="<<max_hit.X()<<" Y="<<max_hit.Y()<<" Z="<<max_hit.Z()<<std::endl;
		first_hit_X.push_back(min_hit.X());
		last_hit_X.push_back(max_hit.X());
		first_hit_Y.push_back(min_hit.Y());
		last_hit_Y.push_back(max_hit.Y());
		first_hit_Z.push_back(min_hit.Z());
		last_hit_Z.push_back(max_hit.Z());
	    }
	}
	else {
	    std::cout<<"NO TRACKS"<<std::endl;
	}

	tree->Fill();		
	first_hit_X.clear();
	last_hit_X.clear();
	first_hit_Y.clear();
	last_hit_Y.clear();
	first_hit_Z.clear();
	last_hit_Z.clear();

	return true;
    }
    // Called at least once.  If multiple file are open, it will be called
    // for each one.   Notice there are two forms...
    void Finalize(CP::TRootOutput * const output) {
	hfile->Write();
	//gSystem->Exec("mv tracks.root "+fName);
    }

private:
  
    TFile* hfile;
    TTree* tree;
  
    std::vector<double> first_hit_X;
    std::vector<double> last_hit_X;

    std::vector<double> first_hit_Y;
    std::vector<double> last_hit_Y;

    std::vector<double> first_hit_Z;
    std::vector<double> last_hit_Z;

    
};

int main(int argc, char **argv) {
    TTreeMakerLoop userCode;
    CP::eventLoop(argc,argv,userCode);
}
