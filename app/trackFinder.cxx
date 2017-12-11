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
#include <TSystem.h>


typedef std::tuple<double,double,double> wTuple;

#ifdef __CINT__
#pragma link C++ class vector+;
#endif

//#define VERBOSE
//#define VVERBOSE

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
    float calculateSlope(std::vector<std::tuple<double,double,double>> trackSeed) {
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

    float minHitDist(std::vector<std::tuple<double,double,double>> t1,std::vector<std::tuple<double,double,double>> t2) {
	float minDist = 999;
	for (auto h1:t1) {
	    for (auto h2:t2) {
		//std::cout<<"dist="<<(std::get<1>(h1)-std::get<1>(h2))*(std::get<1>(h1)-std::get<1>(h2))+(std::get<0>(h1)-std::get<0>(h2))*(std::get<0>(h1)-std::get<0>(h2))<<std::endl;
		if ((std::get<1>(h1)-std::get<1>(h2))*(std::get<1>(h1)-std::get<1>(h2))+(std::get<0>(h1)-std::get<0>(h2))*(std::get<0>(h1)-std::get<0>(h2)) < minDist)
		    minDist = (std::get<1>(h1)-std::get<1>(h2))*(std::get<1>(h1)-std::get<1>(h2))+(std::get<0>(h1)-std::get<0>(h2))*(std::get<0>(h1)-std::get<0>(h2));
	    }	    
	}
	return minDist;
    }

    std::vector<wTuple> instersection(std::vector<wTuple> &v1, std::vector<wTuple> &v2) {
	std::vector<wTuple> v3;
	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());
	set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));
	return v3;
    }


    std::vector<wTuple> combineTuples(std::vector<wTuple> &t1, std::vector<wTuple> &t2) {
	std::vector<wTuple> t12;
	t12.reserve( t1.size() + t2.size() ); // preallocate memory
	t12.insert( t12.end(), t1.begin(), t1.end() );
	t12.insert( t12.end(), t2.begin(), t2.end() );
	return t12;
    }
	    
    TTrackFinder()
    {
	fBeam = false;
	fName = "tracks.root";
    
	first_hit_X.clear();
	first_hit_V.clear();
	first_hit_U.clear();
	last_hit_X.clear();
	last_hit_V.clear();
	last_hit_U.clear();

	beam.clear();
    
	hfile= new TFile("tracks.root","RECREATE");
	tree = new TTree("tracks","");
	tree->Branch("first_hit_X",&first_hit_X);
	tree->Branch("last_hit_X",&last_hit_X);
	tree->Branch("first_hit_V",&first_hit_V);
	tree->Branch("last_hit_V",&last_hit_V);
	tree->Branch("first_hit_U",&first_hit_U);
	tree->Branch("last_hit_U",&last_hit_U);
	tree->Branch("nmatches",&nmatches);
	tree->Branch("nmatches2D",&nmatches2D);
	tree->Branch("Event",&Event);
	tree->Branch("beam",&beam);

    }
    virtual ~TTrackFinder() {}
    
    void Initialize() {
	tracks3D = new TH1F("tracks3D", "",10,0,10);
    
    }

    virtual bool SetOption(std::string option,std::string value="") {
	//if (value != "") return false;
	if (option == "beam") fBeam = true;
	if (option == "name") fName = value;
    
	return true;
    }

    bool operator () (CP::TEvent& event) {
	// Get the list of hits from the event.  The handle is essentially
	// a pointer to the hit selection.
	//CP::THandle<CP::THitSelection> drift(event.GetHits("drift"));
    
	CP::TChannelInfo::Get().SetContext(event.GetContext());
	CP::THandle<CP::TDigitContainer> drift = event.Get<CP::TDigitContainer>("~/digits/drift");

#ifdef VERBOSE
	std::cout<<"=============================================="<<std::endl;
	std::cout<< "Event " << event.GetContext().GetRun()
		 << "." << event.GetContext().GetEvent() << std::endl;
	std::cout<<"=============================================="<<std::endl;
#endif
	Event = event.GetContext().GetEvent();
    
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
		    if (nhWires > 10) hotWires.push_back(hotWire);
		    nhWires = 0;
		}     
		if (w == wiresTimesCharges.back()) {
		    if (nhWires > 10) hotWires.push_back(hotWire);
		}
		hotWire = std::get<0>(w);
	    }
	    int nTracks = 0;
#ifdef VERBOSE
	    std::cout<<"Plane "<<planes[plane]<<std::endl;
	    std::cout<<"wires="<<wiresTimesCharges.size()<<std::endl;
#endif
	    // Veto hot wires
	    std::vector<int> hitsToVeto;
	    for (auto hw:hotWires) {
#ifdef VERBOSE
		std::cout<<"HW="<<hw<<std::endl;      
#endif
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
#ifdef VERBOSE
	    std::cout<<"Track seeds="<<trackSeeds.size()<<std::endl;
#endif
	    for (unsigned int i = 0; i < trackSeeds.size(); i++) {
		// std::cout<<"tseed="<<std::get<0>(trackSeeds[i][0])<<std::endl;
		for (unsigned int k = 0; k < trackSeeds[i].size(); k++) {	
		    // std::cout<<"tseed head="<<std::get<0>(trackSeeds[i][k])<<std::endl;
		    for (int j = wiresTimesCharges.size() -1; j  >= 0; j--) {	    
			double deltaW = fabs(std::get<0>(trackSeeds[i][k])-std::get<0>(wiresTimesCharges[j]));
			double deltaT = fabs(std::get<1>(trackSeeds[i][k])-std::get<1>(wiresTimesCharges[j]));
			// std::cout<<"wire test="<<std::get<0>(wiresTimesCharges[j])<<std::endl;
			// std::cout<<"Deltas="<<deltaW<<" "<<deltaT<<std::endl;
			if ((deltaW !=0 && deltaT != 0) && deltaW < 10 && deltaT < 20) {
			    trackSeeds[i].push_back(wiresTimesCharges[j]);
			    wiresTimesCharges.erase(wiresTimesCharges.begin() + j);
			}	  
		    }	
		}
	    }
	    std::vector<std::vector<wTuple>> pretracks;
#ifdef VERBOSE
	    std::cout<<"pruned seeds="<<trackSeeds.size()<<std::endl;
#endif
	    for (auto ts:trackSeeds) {
#ifdef VVERBOSE
		std::cout<<"track seed size="<<ts.size()<<std::endl;
		if (ts.size() > 1) {
		    for(auto t:ts)
			std::cout<<"wires="<<std::get<0>(t)<<" "<< std::get<1>(t)<<std::endl;
		}
#endif	
		if (ts.size() > 10) {
		    pretracks.push_back(ts);
		    nTracks++;
		}
	    }

	    // Combine tracks by matching slopes and are close by
	    std::vector<std::pair<int,int>> matched_tracks;
	    for (unsigned int i = 0; i< pretracks.size(); i++) {
		std::vector<wTuple> t1 = pretracks[i];
		float slope1 = calculateSlope(t1);
		for (unsigned int j = i; j < pretracks.size(); j++) {
		    std::vector<wTuple> t2 = pretracks[j];
		    float slope2 = calculateSlope(t2);
		    if ( fabs(slope1-slope2)<0.1){
			float minDist = minHitDist(t1,t2);
#ifdef VERBOSE
			std::cout<<"Min dist="<<minDist<<std::endl;
#endif
			if (minDist < 100) {
			    matched_tracks.push_back(std::make_pair(i,j));
#ifdef VERBOSE
			std::cout<<"Match"<<std::endl;
#endif
			}
		    }
		}
	    }

	    //      for (unsigned int i=0; i<pretrack.size(); i++) {

	    std::vector<int> matched_remove;
	    for ( auto p: matched_tracks) {
	
		if (p.first != p.second) {
#ifdef VERBOSE
		    std::cout<<"Pairs="<<p.first<<","<<p.second<<std::endl;
#endif
		    matched_remove.push_back(p.first);
		    matched_remove.push_back(p.second);
		}
	    }

	    for ( auto p: matched_remove) {
#ifdef VERBOSE
		std::cout<<"remove="<<p<<std::endl;
#endif
		matched_tracks.erase(std::remove(matched_tracks.begin(), matched_tracks.end(), std::make_pair(p,p)), matched_tracks.end());
	    }
	    for ( auto p: matched_tracks) {
#ifdef VERBOSE
		std::cout<<"finalw="<<p.first<<" "<<p.second<<std::endl;
#endif
		if (p.first == p.second) tracks.push_back(pretracks[p.first]);
		else tracks.push_back(combineTuples(pretracks[p.first],pretracks[p.second]));
	
	    }      
#ifdef VERBOSE
	    std::cout<<"npretracks = "<<pretracks.size()<<std::endl;
	    std::cout<<"ntracks = "<<tracks.size()<<std::endl;
#endif
	    int ntr = 0;
	    for (auto ts: tracks) {
#ifdef VERBOSE
		std::cout<<"New track"<<std::endl;
#endif
		ntr++;

		double wMax = 0.0; 
		double tMax = -9999.0; 
		double wMin = 9999.0; 
		double tMin = 9999.0; 

		wTuple firstHit;
		wTuple lastHit;

#ifdef VERBOSE
		std::cout<<"Slope="<<calculateSlope(ts)<<std::endl;
#endif      
		for (auto t1:ts) {
		    //for (auto t2:ts) {
		    //if (fabs(std::get<0>(t1) - std::get<0>(t2)) < 3 && fabs(std::get<1>(t1) - std::get<1>(t2)) > 50) {
		    //split = true;
		    //}	  
		    //}	
		    //std::cout<<"Hits="<<std::get<0>(t1)<<" "<<std::get<1>(t1)<<std::endl;
		    if (std::get<0>(t1) > wMax) {wMax = std::get<0>(t1);}
		    if (std::get<1>(t1) > tMax) {tMax = std::get<1>(t1);lastHit = t1; };
		    if (std::get<0>(t1) < wMin) {wMin = std::get<0>(t1);}
		    if (std::get<1>(t1) < tMin) {tMin = std::get<1>(t1); firstHit = t1;};		
		}
#ifdef VERBOSE
		std::cout<<"EDGES="<<std::get<1>(firstHit)<<","<<std::get<0>(firstHit)<<" "<<std::get<1>(lastHit)<<","<<std::get<0>(lastHit)<<std::endl;
#endif	
		trackEdges[plane].push_back(make_tuple(firstHit,lastHit,plane));
	    }
#ifdef VERBOSE
	    std::cout<<"=============================================="<<std::endl;
#endif
	}


	std::sort(trackEdges[0].begin(), trackEdges[0].end(), [] (std::tuple<wTuple,wTuple,int> const& a, std::tuple<wTuple,wTuple,int> const& b) { return std::get<1>(std::get<0>(a)) < std::get<1>(std::get<0>(b)); });
	
	std::sort(trackEdges[1].begin(), trackEdges[1].end(), [] (std::tuple<wTuple,wTuple,int> const& a, std::tuple<wTuple,wTuple,int> const& b) { return std::get<1>(std::get<0>(a)) < std::get<1>(std::get<0>(b)); });
	
	std::sort(trackEdges[2].begin(), trackEdges[2].end(), [] (std::tuple<wTuple,wTuple,int> const& a, std::tuple<wTuple,wTuple,int> const& b) { return std::get<1>(std::get<0>(a)) < std::get<1>(std::get<0>(b)); });

#ifdef VERBOSE
	for (auto ts:trackEdges[0]) {
	    std::cout<<"T="<<std::get<1>(std::get<0>(ts))<<" "<<planes[std::get<2>(ts)]<<std::endl;
	    std::cout<<"T="<<std::get<1>(std::get<1>(ts))<<" "<<planes[std::get<2>(ts)]<<std::endl;
	    std::cout<<"W1="<<std::get<0>(std::get<0>(ts))<<" "<<std::get<0>(std::get<1>(ts))<<std::endl;
	}
	for (auto ts:trackEdges[1]) {
	    std::cout<<"T="<<std::get<1>(std::get<0>(ts))<<" "<<planes[std::get<2>(ts)]<<std::endl;
	    std::cout<<"T="<<std::get<1>(std::get<1>(ts))<<" "<<planes[std::get<2>(ts)]<<std::endl;
	    std::cout<<"W1="<<std::get<0>(std::get<0>(ts))<<" "<<std::get<0>(std::get<1>(ts))<<std::endl;
	}
	for (auto ts:trackEdges[2]) {
	    std::cout<<"T="<<std::get<1>(std::get<0>(ts))<<" "<<planes[std::get<2>(ts)]<<std::endl;
	    std::cout<<"T="<<std::get<1>(std::get<1>(ts))<<" "<<planes[std::get<2>(ts)]<<std::endl;
	    std::cout<<"W1="<<std::get<0>(std::get<0>(ts))<<" "<<std::get<0>(std::get<1>(ts))<<std::endl;
	}
#endif

	int nMatches = 0;
	double fhx,fhv,fhu;
	double lhx,lhv,lhu;

	std::vector<std::tuple<wTuple,wTuple,int>> used_tracks_X;
	std::vector<std::tuple<wTuple,wTuple,int>> used_tracks_V;
	std::vector<std::tuple<wTuple,wTuple,int>> used_tracks_U;

#ifdef VERBOSE
	std::cout<<"tracks X="<<trackEdges[0].size() <<std::endl;
	std::cout<<"tracks V="<<trackEdges[1].size() <<std::endl;
	std::cout<<"tracks U="<<trackEdges[2].size() <<std::endl;
#endif	

	for (auto tsx:trackEdges[0]) {
	    int matchV = 0;
	    int matchU = 0;

	    fhx = std::get<0>(std::get<0>(tsx)); // first hit X plane
	    lhx = std::get<0>(std::get<1>(tsx)); // last hit X plane

	    double delta_xv_f = 100;
	    double delta_xv_l = 100;
	    double delta_xu_f = 100;
	    double delta_xu_l = 100;

	    std::tuple<wTuple,wTuple,int> matched_track_V;
	    std::tuple<wTuple,wTuple,int> matched_track_U;

#ifdef VERBOSE
	    std::cout<<"tracks V 1.5="<<trackEdges[1].size() <<std::endl;
	    std::cout<<"tracks U 1.5="<<trackEdges[2].size() <<std::endl;
#endif
	    
	    for (auto tsv:trackEdges[1]) {	
#ifdef VVERBOSE
		std::cout<<"timesV="<<std::get<1>(std::get<0>(tsx))<<" "<<std::get<1>(std::get<0>(tsv))<<
		    " "<<std::get<1>(std::get<1>(tsx))<<" "<<std::get<1>(std::get<1>(tsv))<<std::endl;

		std::cout<<"deltasV="<<delta_xv_f<<" "<<delta_xv_l<<std::endl;
#endif	    
		if ( fabs(std::get<1>(std::get<0>(tsx)) - std::get<1>(std::get<0>(tsv))) < delta_xv_f ||
		     fabs(std::get<1>(std::get<1>(tsx)) - std::get<1>(std::get<1>(tsv))) < delta_xv_l		     
		     ) {
		    delta_xv_f = fabs(std::get<1>(std::get<0>(tsx)) - std::get<1>(std::get<0>(tsv)));
		    delta_xv_l = fabs(std::get<1>(std::get<1>(tsx)) - std::get<1>(std::get<1>(tsv)));
		    matchV++;
		    fhv = std::get<0>(std::get<0>(tsv)); // first hit V plane
		    lhv = std::get<0>(std::get<1>(tsv)); // last hit V plane
		    matched_track_V = tsv;
		}
	    }
	    for (auto tsu:trackEdges[2]) {
#ifdef VVERBOSE
		std::cout<<"timesU="<<std::get<1>(std::get<0>(tsx))<<" "<<std::get<1>(std::get<0>(tsu))<<
		    " "<<std::get<1>(std::get<1>(tsx))<<" "<<std::get<1>(std::get<1>(tsu))<<std::endl;

		std::cout<<"deltasU="<<delta_xu_f<<" "<<delta_xu_l<<std::endl;
#endif	    
		if ( fabs(std::get<1>(std::get<0>(tsx)) - std::get<1>(std::get<0>(tsu))) < delta_xu_f ||
		     fabs(std::get<1>(std::get<1>(tsx)) - std::get<1>(std::get<1>(tsu))) < delta_xu_l
		     ) {
		    delta_xu_f = fabs(std::get<1>(std::get<0>(tsx)) - std::get<1>(std::get<0>(tsu)));
		    delta_xu_l = fabs(std::get<1>(std::get<1>(tsx)) - std::get<1>(std::get<1>(tsu)));
		    matchU++;
		    fhu = std::get<0>(std::get<0>(tsu)); // first hit U plane
		    lhu = std::get<0>(std::get<1>(tsu)); // last hit U plane
		    matched_track_U = tsu;
		}
	    }
      
	    if (matchV > 0 && matchU > 0){	
#ifdef VERBOSE
		std::cout<<"match time="<<std::get<1>(std::get<0>(tsx)) <<" "<< std::get<1>(std::get<0>(matched_track_V)) << " " << std::get<1>(std::get<0>(matched_track_U)) << std::endl;
#endif
		nMatches++;
		first_hit_X.push_back(fhx);
		first_hit_V.push_back(fhv);
		first_hit_U.push_back(fhu);
		last_hit_X.push_back(lhx);
		last_hit_V.push_back(lhv);
		last_hit_U.push_back(lhu);

		used_tracks_X.push_back(tsx);

		trackEdges[1].erase(std::remove(trackEdges[1].begin(), trackEdges[1].end(), matched_track_V), trackEdges[1].end());
		trackEdges[2].erase(std::remove(trackEdges[2].begin(), trackEdges[2].end(), matched_track_U), trackEdges[2].end());

		if (std::get<1>(std::get<0>(tsx)) > 0 && std::get<1>(std::get<0>(tsx)) < 2000)
		    beam.push_back(true);
		else
		    beam.push_back(false);
	
	    }
	}

	for (auto utx: used_tracks_X) {
	    trackEdges[0].erase(std::remove(trackEdges[0].begin(), trackEdges[0].end(), utx), trackEdges[0].end());
	}
#ifdef VERBOSE
	std::cout<<"tracks X 2="<<trackEdges[0].size() <<std::endl;
	std::cout<<"tracks V 2="<<trackEdges[1].size() <<std::endl;
	std::cout<<"tracks U 2="<<trackEdges[2].size() <<std::endl;
#endif	

	int nMatches_2D = 0;
	for (auto tsx:trackEdges[0]) {
	    fhx = std::get<0>(std::get<0>(tsx)); // first hit X plane
	    lhx = std::get<0>(std::get<1>(tsx)); // last hit X plane	    
	    double delta_xv_f = 100;
	    double delta_xv_l = 100;
      	    std::tuple<wTuple,wTuple,int> matched_track_V;

	    for (auto tsv:trackEdges[1]) {	
		if ( fabs(std::get<1>(std::get<0>(tsx)) - std::get<1>(std::get<0>(tsv))) < delta_xv_f ||
		     fabs(std::get<1>(std::get<1>(tsx)) - std::get<1>(std::get<1>(tsv))) < delta_xv_l		     
		     ) {
		    delta_xv_f = fabs(std::get<1>(std::get<0>(tsx)) - std::get<1>(std::get<0>(tsv)));
		    delta_xv_l = fabs(std::get<1>(std::get<1>(tsx)) - std::get<1>(std::get<1>(tsv)));
		    fhv = std::get<0>(std::get<0>(tsv)); // first hit V plane
		    lhv = std::get<0>(std::get<1>(tsv)); // last hit V plane
		    matched_track_V = tsv;
		    used_tracks_X.push_back(tsx);
		    nMatches_2D++;
		    first_hit_X.push_back(fhx);
		    first_hit_V.push_back(fhv);
		    first_hit_U.push_back(-999);
		    last_hit_X.push_back(lhx);
		    last_hit_V.push_back(lhv);
		    last_hit_U.push_back(-999);

		    if (std::get<1>(std::get<0>(tsx)) > 0 && std::get<1>(std::get<0>(tsx)) < 2000)
			beam.push_back(true);
		    else
			beam.push_back(false);		    
		}
	    }	    
	    trackEdges[1].erase(std::remove(trackEdges[1].begin(), trackEdges[1].end(), matched_track_V), trackEdges[1].end());
	}

	for (auto utx: used_tracks_X) {
	    trackEdges[0].erase(std::remove(trackEdges[0].begin(), trackEdges[0].end(), utx), trackEdges[0].end());
	}
#ifdef VERBOSE
	std::cout<<"tracks X 3="<<trackEdges[0].size() <<std::endl;
	std::cout<<"tracks V 3="<<trackEdges[1].size() <<std::endl;
#endif	
	
	for (auto tsx:trackEdges[0]) {
	    fhx = std::get<0>(std::get<0>(tsx)); // first hit X plane
	    lhx = std::get<0>(std::get<1>(tsx)); // last hit X plane	    
	    double delta_xu_f = 100;
	    double delta_xu_l = 100;
      	    std::tuple<wTuple,wTuple,int> matched_track_U;

	    for (auto tsu:trackEdges[2]) {	
		if ( fabs(std::get<1>(std::get<0>(tsx)) - std::get<1>(std::get<0>(tsu))) < delta_xu_f ||
		     fabs(std::get<1>(std::get<1>(tsx)) - std::get<1>(std::get<1>(tsu))) < delta_xu_l		     
		     ) {
		    delta_xu_f = fabs(std::get<1>(std::get<0>(tsx)) - std::get<1>(std::get<0>(tsu)));
		    delta_xu_l = fabs(std::get<1>(std::get<1>(tsx)) - std::get<1>(std::get<1>(tsu)));
		    fhu = std::get<0>(std::get<0>(tsu)); // first hit U plane
		    lhu = std::get<0>(std::get<1>(tsu)); // last hit U plane
		    matched_track_U = tsu;
		    used_tracks_X.push_back(tsx);
		    nMatches_2D++;
		    first_hit_X.push_back(fhx);
		    first_hit_V.push_back(-999);
		    first_hit_U.push_back(fhu);
		    last_hit_X.push_back(lhx);
		    last_hit_V.push_back(-999);
		    last_hit_U.push_back(lhu);
		    if (std::get<1>(std::get<0>(tsx)) > 0 && std::get<1>(std::get<0>(tsx)) < 2000)
			beam.push_back(true);
		    else
			beam.push_back(false);	    
		}
	    }
	    trackEdges[2].erase(std::remove(trackEdges[2].begin(), trackEdges[2].end(), matched_track_U), trackEdges[2].end());
	}

	for (auto utx: used_tracks_X) {
	    trackEdges[0].erase(std::remove(trackEdges[0].begin(), trackEdges[0].end(), utx), trackEdges[0].end());
	}
#ifdef VERBOSE
	std::cout<<"tracks X 4="<<trackEdges[0].size() <<std::endl;
	std::cout<<"tracks V 4="<<trackEdges[2].size() <<std::endl;
#endif
	for (auto tsu:trackEdges[2]) {
	    fhu = std::get<0>(std::get<0>(tsu)); // first hit U plane
	    lhu = std::get<0>(std::get<1>(tsu)); // last hit U plane	    
	    double delta_uv_f = 100;
	    double delta_uv_l = 100;
      	    std::tuple<wTuple,wTuple,int> matched_track_V;

	    for (auto tsv:trackEdges[1]) {	
		if ( fabs(std::get<1>(std::get<0>(tsu)) - std::get<1>(std::get<0>(tsv))) < delta_uv_f ||
		     fabs(std::get<1>(std::get<1>(tsu)) - std::get<1>(std::get<1>(tsv))) < delta_uv_l		     
		     ) {
		    delta_uv_f = fabs(std::get<1>(std::get<0>(tsu)) - std::get<1>(std::get<0>(tsv)));
		    delta_uv_l = fabs(std::get<1>(std::get<1>(tsu)) - std::get<1>(std::get<1>(tsv)));
		    fhv = std::get<0>(std::get<0>(tsv)); // first hit V plane
		    lhv = std::get<0>(std::get<1>(tsv)); // last hit V plane
		    matched_track_V = tsv;
		    used_tracks_U.push_back(tsu);
		    nMatches_2D++;
		    first_hit_X.push_back(-999);
		    first_hit_V.push_back(fhv);
		    first_hit_U.push_back(fhu);
		    last_hit_X.push_back(-999);
		    last_hit_V.push_back(lhv);
		    last_hit_U.push_back(lhu);
		    if (std::get<1>(std::get<0>(tsu)) > 0 && std::get<1>(std::get<0>(tsu)) < 2000)
			beam.push_back(true);
		    else
			beam.push_back(false);
		}
	    }
	    trackEdges[1].erase(std::remove(trackEdges[1].begin(), trackEdges[1].end(), matched_track_V), trackEdges[1].end());	    
	}

	for (auto utu: used_tracks_U) {
	    trackEdges[2].erase(std::remove(trackEdges[2].begin(), trackEdges[2].end(), utu), trackEdges[2].end());
	}
#ifdef VERBOSE
	std::cout<<"tracks U 5="<<trackEdges[2].size() <<std::endl;
	std::cout<<"tracks V 5="<<trackEdges[1].size() <<std::endl;
#endif	
	

	nmatches = nMatches;
	nmatches2D = nMatches_2D;
	tree->Fill();
	tracks3D->Fill(nMatches);
#ifdef VERBOSE
	std::cout<<"nMatches="<<nMatches<<std::endl;
	std::cout<<"nMatches 2D="<<nMatches_2D<<std::endl;

	std::cout<<"=============================================="<<std::endl;
	std::cout<<"=============================================="<<std::endl;
#endif
	first_hit_X.clear();
	first_hit_V.clear();
	first_hit_U.clear();
	last_hit_X.clear();
	last_hit_V.clear();
	last_hit_U.clear();
	beam.clear();
    
	return true;
    }

    // Called at least once.  If multiple file are open, it will be called
    // for each one.   Notice there are two forms...
    void Finalize(CP::TRootOutput * const output) {
	tracks3D->Draw();
	gPad->Print("XUV_tracks.root");
	std::cout<<"FNAME="<<fName<<std::endl;
	hfile->Write();
	gSystem->Exec("mv tracks.root "+fName);
    }
private:
    bool fBeam;
    TString fName;
    
    TH1F* tracks3D;

    TFile* hfile;
    TTree* tree;
    Int_t nmatches;
    Int_t nmatches2D;
    Int_t Event;
  
    std::vector<bool> beam;
    std::vector<double> first_hit_X;
    std::vector<double> last_hit_X;
    std::vector<double> first_hit_V;
    std::vector<double> last_hit_V;
    std::vector<double> first_hit_U;
    std::vector<double> last_hit_U;
  
};
int main(int argc, char **argv) {
    TTrackFinder userCode;
    CP::eventLoop(argc,argv,userCode);
}
