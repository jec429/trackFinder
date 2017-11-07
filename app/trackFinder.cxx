#include <memory>
#include <math.h>
#include <vector>
#include <tuple>
// Include to get Event Loop.
#include <eventLoop.hxx>
//#include "TPlotDigitsHits.hxx"
//#include "TEventDisplay.hxx"
//#include "TGUIManager.hxx"
     
#include <HEPUnits.hxx>
#include <TCaptLog.hxx>
#include <CaptGeomId.hxx>
#include <TEvent.hxx>
#include <TEventFolder.hxx>
#include <TPulseDigit.hxx>
#include <TCalibPulseDigit.hxx>
#include <TMCChannelId.hxx>
#include <TRuntimeParameters.hxx>
#include <TUnitsTable.hxx>
    
#include <TChannelInfo.hxx>
#include <TChannelCalib.hxx>
#include <TGeometryInfo.hxx>
// Includes for ROOT classes
#include <TH1F.h>
#include <TH2F.h>
#include <TPad.h>
#include <TGraph.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>


typedef std::tuple<double,double,double> wTuple;

#ifdef __CINT__
#pragma link C++ class wTuple+;
#endif


class TTrackFinder: public CP::TEventLoopFunction {
public:
  std::size_t GetDigitSampleCount(const CP::TDigit* d) {
    const CP::TPulseDigit* pulse 
      = dynamic_cast<const CP::TPulseDigit*>(d);
    if (pulse) return pulse->GetSampleCount();
    const CP::TCalibPulseDigit* calib 
      = dynamic_cast<const CP::TCalibPulseDigit*>(d);
    if (calib) return calib->GetSampleCount();
    return 0;
  }
  double GetDigitFirstTime(const CP::TDigit* d) {
    const CP::TPulseDigit* pulse 
      = dynamic_cast<const CP::TPulseDigit*>(d);
    if (pulse) return pulse->GetFirstSample();
    const CP::TCalibPulseDigit* calib 
      = dynamic_cast<const CP::TCalibPulseDigit*>(d);
    if (calib) return calib->GetFirstSample()/unit::microsecond;
    return 0.0;
  }
  // This will be 1 for raw digits and 500ns for calibrated digits.
  double GetDigitSampleStep(const CP::TDigit* d) {
    double diff = GetDigitLastTime(d) - GetDigitFirstTime(d);
    return diff/GetDigitSampleCount(d);
  }
  double GetDigitLastTime(const CP::TDigit* d) {
    const CP::TPulseDigit* pulse 
      = dynamic_cast<const CP::TPulseDigit*>(d);
    if (pulse) return pulse->GetFirstSample()+pulse->GetSampleCount();
    const CP::TCalibPulseDigit* calib 
      = dynamic_cast<const CP::TCalibPulseDigit*>(d);
    if (calib) return calib->GetLastSample()/unit::microsecond;
    return 0.0;
  }
  double GetDigitSample(const CP::TDigit* d, int i) {
    const CP::TPulseDigit* pulse 
      = dynamic_cast<const CP::TPulseDigit*>(d);
    if (pulse) return pulse->GetSample(i);
    const CP::TCalibPulseDigit* calib 
      = dynamic_cast<const CP::TCalibPulseDigit*>(d);
    if (calib) return calib->GetSample(i);
    return 0;
  }
  // This will be 3200 for raw digits and 0 for calibrated digits.
  double GetDigitTriggerOffset(const CP::TDigit* d) {
    const CP::TPulseDigit* pulse 
      = dynamic_cast<const CP::TPulseDigit*>(d);
    if (!pulse) return 0.0;
    CP::TChannelCalib chanCalib;
    double off = chanCalib.GetTimeConstant(pulse->GetChannelId(),0);
    double tim = chanCalib.GetTimeConstant(pulse->GetChannelId(),1);
    return - off/tim;
  }
  float calculateSlope(std::vector<std::tuple<double,double,CP::THandle<CP::THit>>> trackSeed){
    float slope = 0.;
    for (auto t1:trackSeed) {
      for (auto t2:trackSeed) {
	if (std::get<0>(t1) == std::get<0>(t2)) continue;
	slope += (std::get<1>(t1)-std::get<1>(t2))/(std::get<0>(t1)-std::get<0>(t2));
      }
    }
    slope = slope/(trackSeed.size()*(trackSeed.size()-1));
    return slope;
  }
  TTrackFinder():
    first_hit(new wTuple),
    last_hit(new wTuple)
  {
    fBeam = false;
    
    hfile= new TFile("tracks.root","RECREATE");
    tree = new TTree("tracks","");
    tree->Branch("first_hit",first_hit);
    tree->Branch("last_hit",last_hit);
    tree->Branch("wplane",&wplane);
    tree->Branch("Event",&Event);

  }
  virtual ~TTrackFinder() {}
    
  void Initialize() {
    tracks3D = new TH1F("tracks3D", "",10,0,10);
    
  }

  virtual bool SetOption(std::string option,std::string value="") {
    if (value != "") return false;
    if (option == "beam") fBeam = true;
    
    return true;
  }

  bool operator () (CP::TEvent& event) {
    // Get the list of hits from the event.  The handle is essentially
    // a pointer to the hit selection.
    //CP::THandle<CP::THitSelection> drift(event.GetHits("drift"));
    
    CP::TChannelInfo::Get().SetContext(event.GetContext());
    CP::THandle<CP::TDigitContainer> drift = event.Get<CP::TDigitContainer>("~/digits/drift");

    std::cout<<"=============================================="<<std::endl;
    std::cout<< "Event " << event.GetContext().GetRun()
	     << "." << event.GetContext().GetEvent() << std::endl;
    std::cout<<"=============================================="<<std::endl;

    //int plane = 0;
    std::map<int,std::string> planes;
    planes[0] = "X";
    planes[1] = "V";
    planes[2] = "U";

    std::vector<std::tuple<wTuple,wTuple,int>> trackEdges[3];
    
    for(int plane=0;plane<3;plane++){
      double signalStart = 1E+6;
      double signalEnd = -1E+6;
      double wireTimeStep = -1.0;
      double digitSampleOffset = 0;
      double digitSampleStep = -1;
      double triggerOffset = digitSampleOffset;
      CP::TChannelCalib chanCalib;

      // Find the Z axis range for the histogram.  The median samples the
      // middle.  The max and min are the "biggest" distance from the median.
      std::vector<double> samples;
      double medianSample = 0.0;
      
      if (drift) {
	for (CP::TDigitContainer::const_iterator d = drift->begin(); d != drift->end(); ++d) {
	  // Figure out if this is in the right plane, and get the wire
	  // number.
	  const CP::TDigit* digit = dynamic_cast<const CP::TDigit*>(*d);
	  if (!digit) continue;
	  CP::TGeometryId id = CP::TChannelInfo::Get().GetGeometry(digit->GetChannelId());
	  if (CP::GeomId::Captain::GetWirePlane(id) != plane) continue;
	  // Save the sample to find the median.
	  for (std::size_t i = 0; i < GetDigitSampleCount(*d); ++i) {
	    double s = GetDigitSample(*d,i);
	    if (!std::isfinite(s)) continue;
	    samples.push_back(s);
	  }
	  if (digitSampleStep < 0) {
	    // Find the time range.
	    digitSampleStep = GetDigitSampleStep(*d);
	    digitSampleOffset = GetDigitTriggerOffset(*d);
	  }
	  if (wireTimeStep < 0.0) {
	    wireTimeStep=chanCalib.GetTimeConstant((*d)->GetChannelId(),1);
	  }
	}      
	     
	// Crash prevention.  It shouldn't be possible to have digits without
	// any samples, but...
	if (samples.empty()) return 0;
           
	std::sort(samples.begin(),samples.end());
	medianSample = samples[0.5*samples.size()];
   
	double maxSample = std::abs(samples[0.99*samples.size()]-medianSample);
	double s = std::abs(samples[0.01*samples.size()]-medianSample);
	maxSample = std::max(maxSample,s);
           
	// Find the time axis range based on the times of the bins with a
	// signal.  In most files, the result will be the full range of the
	// digitizer 0-9595, or the full range of the calibrated times (-1600
	// usec to 3197 usec), but things get a little more complicated for the
	// MC.
	std::vector<double> times;
	for (CP::TDigitContainer::const_iterator d = drift->begin();
	     d != drift->end(); ++d) {
	  // Figure out if this is in the right plane, and get the wire
	  // number.
	  const CP::TDigit* digit 
	    = dynamic_cast<const CP::TDigit*>(*d);
	  if (!digit) continue;
	  CP::TGeometryId id 
	    = CP::TChannelInfo::Get().GetGeometry(digit->GetChannelId());
	  if (CP::GeomId::Captain::GetWirePlane(id) != plane) continue;
	  double maxSignal = 0.0;
	  for (std::size_t i = 0; i < GetDigitSampleCount(*d); ++i) {
	    double s = std::abs(GetDigitSample(*d,i) - medianSample);
	    if (!std::isfinite(s)) continue;
	    if (maxSignal < s) maxSignal = s;
	  }
	  if (maxSignal < 0.25*maxSample) continue;
	  times.push_back(GetDigitFirstTime(*d));
	  times.push_back(GetDigitLastTime(*d));
	}
	std::sort(times.begin(),times.end());
       
	signalStart = times[0.01*times.size()];
	signalEnd = times[0.99*times.size()];
      }
      else if (event.Get<CP::THitSelection>("~/hits/drift")) {
	CP::THandle<CP::THitSelection> hits = event.Get<CP::THitSelection>("~/hits/drift");
	for (CP::THitSelection::iterator h = hits->begin(); h != hits->end(); ++h) {
	  if (signalStart>(*h)->GetTime()) signalStart = (*h)->GetTime();
	  if (signalEnd<(*h)->GetTime()) signalEnd = (*h)->GetTime();
	  CP::TChannelId cid = (*h)->GetChannelId();
	  if (wireTimeStep < 0) {
	    wireTimeStep = chanCalib.GetTimeConstant(cid,1);
	  }
	}
	digitSampleStep = 0.5;
	digitSampleOffset = 0.0;
	signalStart /= unit::microsecond;
	signalEnd /= unit::microsecond;
      }    
    
      CP::THandle<CP::THitSelection> hits = event.Get<CP::THitSelection>("~/hits/drift");
    
      double timeUnit =  wireTimeStep/digitSampleStep;

      std::vector<wTuple> wiresTimesCharges;
    
      for (CP::THitSelection::iterator h = hits->begin(); h != hits->end(); ++h) {
	CP::TGeometryId id = (*h)->GetGeomId();
	CP::TChannelId cid = (*h)->GetChannelId();
	if (CP::GeomId::Captain::GetWirePlane(id) != plane) continue;
	// The wire number (offset for the middle of the bin).
	double wire = CP::GeomId::Captain::GetWireNumber(id) + 0.5;
	// The hit charge
	double charge = (*h)->GetCharge();
	// The hit time.
	double hTime = (*h)->GetTime();
	hTime = hTime;
	// The digitized hit time.
	double dTime = hTime/timeUnit + triggerOffset;
	// The digitized hit start time.
	double dStartTime = ((*h)->GetTimeStart()-(*h)->GetTime());
	dStartTime /= timeUnit;
	dStartTime += dTime;
	// The digitized hit start time.
	double dStopTime = ((*h)->GetTimeStop()-(*h)->GetTime());
	dStopTime /= timeUnit;
	dStopTime += dTime;
	// The hit RMS.
	//double rms = (*h)->GetTimeRMS();
	// The digitized RMS
	//double dRMS = rms/timeUnit;

	wiresTimesCharges.push_back(std::make_tuple(wire,dTime,charge));	            
      }

      std::sort(wiresTimesCharges.begin(), wiresTimesCharges.end(), [] (wTuple const& a, wTuple const& b) { return std::get<0>(a) < std::get<0>(b); });
    
      std::vector<std::vector<wTuple>> trackSeeds;
      std::vector<std::vector<wTuple>> tracks;
      std::vector<double> hotWires;
      double hotWire = 0;
      int nhWires = 0;
    
      for (auto w:wiresTimesCharges) {
	std::vector<wTuple> track_tmp;
	track_tmp.push_back(w);
	trackSeeds.push_back(track_tmp);
	if(std::get<0>(w) == hotWire) {
	  nhWires++;
	}
	else{
	  if (nhWires > 20) hotWires.push_back(hotWire);
	  nhWires = 0;
	}     
	if (w == wiresTimesCharges.back()) {
	  if (nhWires > 20) hotWires.push_back(hotWire);
	}
	hotWire = std::get<0>(w);
      }
      int nTracks = 0;
      std::cout<<"Plane "<<planes[plane]<<std::endl;
      std::cout<<"wires="<<wiresTimesCharges.size()<<std::endl;   
      // Veto hot wires
      std::vector<int> hitsToVeto;
      for (auto hw:hotWires) {
	for (unsigned int i = 0; i < trackSeeds.size(); i++) {
	  if (hw == std::get<0>(trackSeeds[i][0])) {
	    hitsToVeto.push_back(i);
	  }
	}
      }

      reverse(hitsToVeto.begin(),hitsToVeto.end());
    
      for (auto i:hitsToVeto) {
	trackSeeds.erase(trackSeeds.begin() + i);
	wiresTimesCharges.erase(wiresTimesCharges.begin() + i);
      }
      std::cout<<"track seeds="<<trackSeeds.size()<<std::endl;
        
      for (unsigned int i = 0; i < trackSeeds.size(); i++) {
	for (unsigned int k = 0; k < trackSeeds[i].size(); k++) {	
	  for (int j = wiresTimesCharges.size() -1; j  >= 0; j--) {	  
	    double deltaW = fabs(std::get<0>(trackSeeds[i][k])-std::get<0>(wiresTimesCharges[j]));
	    double deltaT = fabs(std::get<1>(trackSeeds[i][k])-std::get<1>(wiresTimesCharges[j]));
	    //std::cout<<"Deltas="<<deltaW<<" "<<deltaT<<std::endl;
	    if ((deltaW !=0 && deltaT != 0) && deltaW < 6 && deltaT < 50) {
	      trackSeeds[i].push_back(wiresTimesCharges[j]);
	      wiresTimesCharges.erase(wiresTimesCharges.begin() + j);
	    }	  
	  }	
	}
      }
      std::cout<<"pruned seeds="<<trackSeeds.size()<<std::endl;
      for (auto ts:trackSeeds) {
	// std::cout<<"track seed size="<<ts.size()<<std::endl;
	// if (ts.size() > 1) {
	//   for(auto t:ts)
	//     std::cout<<"wires="<<std::get<0>(t)<<std::endl;
	// }
	if (ts.size() > 10) {
	  tracks.push_back(ts);
	  nTracks++;
	}
      }
        
      std::cout<<"ntracks = "<<nTracks<<std::endl;
      int ntr = 0;
      for (auto ts: tracks) {
	std::cout<<"New track"<<std::endl;
	ntr++;
	bool split = false;
	bool shower = false;

	double wMax = 0.0; 
	double tMax = -9999.0; 
	double wMin = 9999.0; 
	double tMin = 9999.0; 

	wTuple firstHit;
	wTuple lastHit;
      
	for (auto t1:ts) {
	  for (auto t2:ts) {
	    if (fabs(std::get<0>(t1) - std::get<0>(t2)) < 3 && fabs(std::get<1>(t1) - std::get<1>(t2)) > 50) {
	      split = true;
	    }	  
	  }	
	  //std::cout<<"Hits="<<std::get<0>(t1)<<" "<<std::get<1>(t1)<<std::endl;
	  if (std::get<0>(t1) > wMax) wMax = std::get<0>(t1);
	  if (std::get<1>(t1) > tMax) {tMax = std::get<1>(t1); lastHit = t1;};
	  if (std::get<0>(t1) < wMin) wMin = std::get<0>(t1);
	  if (std::get<1>(t1) < tMin) {tMin = std::get<1>(t1); firstHit = t1;};		
	}
	std::cout<<"EDGES="<<std::get<1>(firstHit)<<","<<std::get<0>(firstHit)<<" "<<std::get<1>(lastHit)<<","<<std::get<0>(lastHit)<<std::endl;
	trackEdges[plane].push_back(make_tuple(firstHit,lastHit,plane));
      }        
      std::cout<<"=============================================="<<std::endl;
    }


    std::sort(trackEdges[0].begin(), trackEdges[0].end(), [] (std::tuple<wTuple,wTuple,int> const& a, std::tuple<wTuple,wTuple,int> const& b) { return std::get<1>(std::get<0>(a)) < std::get<1>(std::get<0>(b)); });
    std::sort(trackEdges[1].begin(), trackEdges[1].end(), [] (std::tuple<wTuple,wTuple,int> const& a, std::tuple<wTuple,wTuple,int> const& b) { return std::get<1>(std::get<0>(a)) < std::get<1>(std::get<0>(b)); });
    std::sort(trackEdges[2].begin(), trackEdges[2].end(), [] (std::tuple<wTuple,wTuple,int> const& a, std::tuple<wTuple,wTuple,int> const& b) { return std::get<1>(std::get<0>(a)) < std::get<1>(std::get<0>(b)); });
    
    for (auto ts:trackEdges[0]) {
      std::cout<<"T="<<std::get<1>(std::get<0>(ts))<<" "<<planes[std::get<2>(ts)]<<std::endl;
      std::cout<<"W1="<<std::get<0>(std::get<0>(ts))<<" "<<std::get<0>(std::get<1>(ts))<<std::endl;
    }
    for (auto ts:trackEdges[1]) {
      std::cout<<"T="<<std::get<1>(std::get<0>(ts))<<" "<<planes[std::get<2>(ts)]<<std::endl;
      std::cout<<"W1="<<std::get<0>(std::get<0>(ts))<<" "<<std::get<0>(std::get<1>(ts))<<std::endl;
    }
    for (auto ts:trackEdges[2]) {
      std::cout<<"T="<<std::get<1>(std::get<0>(ts))<<" "<<planes[std::get<2>(ts)]<<std::endl;
      std::cout<<"W1="<<std::get<0>(std::get<0>(ts))<<" "<<std::get<0>(std::get<1>(ts))<<std::endl;
    }

    int nMatches = 0;
    for (auto tsx:trackEdges[0]) {
      int match = 1;
      for (auto tsv:trackEdges[1]) {
	if ( fabs(std::get<1>(std::get<0>(tsx)) - std::get<1>(std::get<0>(tsv))) < 50 )
	  match++;
      }
      for (auto tsu:trackEdges[2]) {
	if ( fabs(std::get<1>(std::get<0>(tsx)) - std::get<1>(std::get<0>(tsu))) < 50 )
	  match++;	
      }

      if (match==3){
	std::cout<<"match time="<<std::get<1>(std::get<0>(tsx))<<std::endl;
	nMatches++;
      }
    }

    tracks3D->Fill(nMatches);
    std::cout<<"nMatches="<<nMatches<<std::endl;

    std::cout<<"=============================================="<<std::endl;
    std::cout<<"=============================================="<<std::endl;

    return true;
  }

  // Called at least once.  If multiple file are open, it will be called
  // for each one.   Notice there are two forms...
  void Finalize(CP::TRootOutput * const output) {
    tracks3D->Draw();
    gPad->Print("XUV_tracks.root");
    //hfile->Write();
  }
private:
  bool fBeam;
  
  TH1F* tracks3D;

  TFile* hfile;
  TTree* tree;
  Int_t wplane;
  Int_t Event;
  wTuple * first_hit;
  wTuple * last_hit;
  
};
int main(int argc, char **argv) {
  TTrackFinder userCode;
  CP::eventLoop(argc,argv,userCode);
}
