#include <assert.h>
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include <DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h>
#include <DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h>
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

using namespace std;
using namespace edm;


struct MuonData
{
  void init();
  TTree* book(TTree *t);

  bool has_prop_GE11_L1;
  bool has_prop_GE11_L2;
  int prop_chamber_GE11_L1;
  int prop_chamber_GE11_L2;
  int prop_roll_GE11_L1;
  int prop_roll_GE11_L2;
  float prop_x_GE11_L1;
  float prop_x_GE11_L2;
  float prop_y_GE11_L1;
  float prop_y_GE11_L2;
  float prop_r_GE11_L1;
  float prop_r_GE11_L2;
  float prop_localx_GE11_L1;
  float prop_localx_GE11_L2;
  float prop_localy_GE11_L1;
  float prop_localy_GE11_L2;
  float prop_y_adjusted_GE11_L1;
  float prop_y_adjusted_GE11_L2;
  float prop_localphi_rad_GE11_L1;
  float prop_localphi_rad_GE11_L2; 
  float prop_localphi_deg_GE11_L1;
  float prop_localphi_deg_GE11_L2;

  bool has_rechit_GE11_L1;
  bool has_rechit_GE11_L2;
  int rechit_chamber_GE11_L1;
  int rechit_chamber_GE11_L2;
  int rechit_roll_GE11_L1;
  int rechit_roll_GE11_L2;
  float rechit_x_GE11_L1;
  float rechit_x_GE11_L2;
  float rechit_y_GE11_L1;
  float rechit_y_GE11_L2;
  float rechit_r_GE11_L1;
  float rechit_r_GE11_L2;
  float rechit_localx_GE11_L1;
  float rechit_localx_GE11_L2;
  float rechit_localy_GE11_L1;
  float rechit_localy_GE11_L2;
  float rechit_y_adjusted_GE11_L1;
  float rechit_y_adjusted_GE11_L2;
  float rechit_localphi_rad_GE11_L1;
  float rechit_localphi_rad_GE11_L2;
  float rechit_localphi_deg_GE11_L1;
  float rechit_localphi_deg_GE11_L2;

  float RdPhi_GE11_L1;
  float RdPhi_GE11_L2;

  bool has_fidcut_GE11_L1;
  bool has_fidcut_GE11_L2;


};

void MuonData::init()
{
  has_prop_GE11_L1 = false;
  has_prop_GE11_L2 = false;
  prop_chamber_GE11_L1 = 99999;
  prop_chamber_GE11_L2 = 99999;
  prop_roll_GE11_L1 = 99999;
  prop_roll_GE11_L2 = 99999;
  prop_x_GE11_L1 = 99999;
  prop_x_GE11_L2 = 99999;
  prop_y_GE11_L1 = 99999;
  prop_y_GE11_L2 = 99999;
  prop_r_GE11_L1 = 99999;
  prop_r_GE11_L2 = 99999;
  prop_localx_GE11_L1 = 99999;
  prop_localx_GE11_L2 = 99999;
  prop_localy_GE11_L1 = 99999;
  prop_localy_GE11_L2 = 99999;
  prop_y_adjusted_GE11_L1 = 99999;
  prop_y_adjusted_GE11_L2 = 99999;
  prop_localphi_rad_GE11_L1 = 99999;
  prop_localphi_rad_GE11_L2 = 99999;
  prop_localphi_deg_GE11_L1 = 99999;
  prop_localphi_deg_GE11_L2 = 99999;

  has_rechit_GE11_L1 = false;
  has_rechit_GE11_L2 = false;
  rechit_chamber_GE11_L1 = 999999;
  rechit_chamber_GE11_L2 = 999999;
  rechit_roll_GE11_L1 = 999999;
  rechit_roll_GE11_L2 = 999999;
  rechit_x_GE11_L1 = 999999;
  rechit_x_GE11_L2 = 999999;
  rechit_y_GE11_L1 = 999999;
  rechit_y_GE11_L2 = 999999;
  rechit_r_GE11_L1 = 999999;
  rechit_r_GE11_L2 = 999999;
  rechit_localx_GE11_L1 = 999999;
  rechit_localx_GE11_L2 = 999999;
  rechit_localy_GE11_L1 = 999999;
  rechit_localy_GE11_L2 = 999999;
  rechit_y_adjusted_GE11_L1 = 999999;
  rechit_y_adjusted_GE11_L2 = 999999;
  rechit_localphi_rad_GE11_L1 = 999999;
  rechit_localphi_rad_GE11_L2 = 999999;
  rechit_localphi_deg_GE11_L1 = 999999;
  rechit_localphi_deg_GE11_L2 = 999999; 

  RdPhi_GE11_L1 = 999999;
  RdPhi_GE11_L2 = 999999;

  has_fidcut_GE11_L1 = false;
  has_fidcut_GE11_L2 = false;

}

TTree* MuonData::book(TTree *t){
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("MuonData", "MuonData");


//Propogated Data
  t->Branch("has_prop_GE11_L1", &has_prop_GE11_L1);//, "has_prop_GE11_L1/B");
  t->Branch("has_prop_GE11_L2", &has_prop_GE11_L2);//, "has_prop_GE11_L2/B");
  t->Branch("prop_chamber_GE11_L1", &prop_chamber_GE11_L1);//, "prop_chamber_GE11_L1/I");
  t->Branch("prop_chamber_GE11_L2", &prop_chamber_GE11_L2);//, "prop_chamber_GE11_L2/I");
  t->Branch("prop_roll_GE11_L1", &prop_roll_GE11_L1);//, "prop_roll_GE11_L1/I");
  t->Branch("prop_roll_GE11_L2", &prop_roll_GE11_L2);//, "prop_roll_GE11_L2/I");
  t->Branch("prop_x_GE11_L1", &prop_x_GE11_L1);//, "prop_x_GE11_L1/F");
  t->Branch("prop_x_GE11_L2", &prop_x_GE11_L2);//, "prop_x_GE11_L2/F");
  t->Branch("prop_y_GE11_L1", &prop_y_GE11_L1);//, "prop_y_GE11_L1/F");
  t->Branch("prop_y_GE11_L2", &prop_y_GE11_L2);//, "prop_y_GE11_L2/F");
  t->Branch("prop_r_GE11_L1", &prop_r_GE11_L1);//, "prop_r_GE11_L1/F");
  t->Branch("prop_r_GE11_L2", &prop_r_GE11_L2);//, "prop_r_GE11_L2/F");
  t->Branch("prop_localx_GE11_L1", &prop_localx_GE11_L1);//, "prop_localx_GE11_L1/F");
  t->Branch("prop_localx_GE11_L2", &prop_localx_GE11_L2);//, "prop_localx_GE11_L2/F");
  t->Branch("prop_localy_GE11_L1", &prop_localy_GE11_L1);//, "prop_localy_GE11_L1/F");
  t->Branch("prop_localy_GE11_L2", &prop_localy_GE11_L2);//, "prop_localy_GE11_L2/F");
  t->Branch("prop_y_adjusted_GE11_L1", &prop_y_adjusted_GE11_L1);//, "prop_y_adjusted_GE11_L1/F");
  t->Branch("prop_y_adjusted_GE11_L2", &prop_y_adjusted_GE11_L2);//, "prop_y_adjusted_GE11_L2/F");
  t->Branch("prop_localphi_rad_GE11_L1", &prop_localphi_rad_GE11_L1);//, "prop_loalphi_rad_GE11_L1/F");
  t->Branch("prop_localphi_rad_GE11_L2", &prop_localphi_rad_GE11_L2);//, "prop_loalphi_rad_GE11_L2/F");
  t->Branch("prop_localphi_deg_GE11_L1", &prop_localphi_deg_GE11_L1);//, "prop_loalphi_deg_GE11_L1/F");
  t->Branch("prop_localphi_deg_GE11_L2", &prop_localphi_deg_GE11_L2);//, "prop_loalphi_deg_GE11_L2/F");
//Reconstructed Data
  t->Branch("has_rechit_GE11_L1", &has_rechit_GE11_L1);//, "has_rechit_GE11_L1/B");
  t->Branch("has_rechit_GE11_L2", &has_rechit_GE11_L2);//, "has_rechit_GE11_L2/B");
  t->Branch("rechit_chamber_GE11_L1", &rechit_chamber_GE11_L1);//, "rechit_chamber_GE11_L1/I");
  t->Branch("rechit_chamber_GE11_L2", &rechit_chamber_GE11_L2);//, "rechit_chamber_GE11_L2/I");
  t->Branch("rechit_roll_GE11_L1", &rechit_roll_GE11_L1);//, "rechit_roll_GE11_L1/I");
  t->Branch("rechit_roll_GE11_L2", &rechit_roll_GE11_L2);//, "rechit_roll_GE11_L2/I");
  t->Branch("rechit_x_GE11_L1", &rechit_x_GE11_L1);//, "rechit_x_GE11_L1/F");
  t->Branch("rechit_x_GE11_L2", &rechit_x_GE11_L2);//, "rechit_x_GE11_L2/F");
  t->Branch("rechit_y_GE11_L1", &rechit_y_GE11_L1);//, "rechit_y_GE11_L1/F");
  t->Branch("rechit_y_GE11_L2", &rechit_y_GE11_L2);//, "rechit_y_GE11_L2/F");
  t->Branch("rechit_r_GE11_L1", &rechit_r_GE11_L1);//, "rechit_r_GE11_L1/F");
  t->Branch("rechit_r_GE11_L2", &rechit_r_GE11_L2);//, "rechit_r_GE11_L2/F");
  t->Branch("rechit_localx_GE11_L1", &rechit_localx_GE11_L1);//, "rechit_localx_GE11_L1/F");
  t->Branch("rechit_localx_GE11_L2", &rechit_localx_GE11_L2);//, "rechit_localx_GE11_L2/F");
  t->Branch("rechit_localy_GE11_L1", &rechit_localy_GE11_L1);//, "rechit_localy_GE11_L1/F");
  t->Branch("rechit_localy_GE11_L2", &rechit_localy_GE11_L2);//, "rechit_localy_GE11_L2/F");
  t->Branch("rechit_y_adjusted_GE11_L1", &rechit_y_adjusted_GE11_L1);//, "rechit_y_adjusted_GE11_L1/F");
  t->Branch("rechit_y_adjusted_GE11_L2", &rechit_y_adjusted_GE11_L2);//, "rechit_y_adjusted_GE11_L2/F");
  t->Branch("rechit_localphi_rad_GE11_L1", &rechit_localphi_rad_GE11_L1);//, "rechit_loalphi_rad_GE11_L1/F");
  t->Branch("rechit_localphi_rad_GE11_L2", &rechit_localphi_rad_GE11_L2);//, "rechit_loalphi_rad_GE11_L2/F");
  t->Branch("rechit_localphi_deg_GE11_L1", &rechit_localphi_deg_GE11_L1);//, "rechit_loalphi_deg_GE11_L1/F");
  t->Branch("rechit_localphi_deg_GE11_L2", &rechit_localphi_deg_GE11_L2);//, "rechit_loalphi_deg_GE11_L2/F");
//Residual
  t->Branch("RdPhi_GE11_L1", &RdPhi_GE11_L1);//, "RdPhi_GE11_L1/F");
  t->Branch("RdPhi_GE11_L2", &RdPhi_GE11_L2);//, "RdPhi_GE11_L2/F");
//Cut
  t->Branch("has_fidcut_GE11_L1", &has_fidcut_GE11_L1);//, "has_fidcut_GE11_L1/B");
  t->Branch("has_fidcut_GE11_L2", &has_fidcut_GE11_L2);//, "has_fidcut_GE11_L2/B");

  return t;
}

class analyser : public edm::EDAnalyzer {
public:
  explicit analyser(const edm::ParameterSet&);
  ~analyser(){};

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() ;
  virtual void endJob() ;

  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;

  edm::Service<TFileService> fs;

  MuonServiceProxy* theService_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;

  edm::ESHandle<GEMGeometry> GEMGeometry_;

  TTree * tree_data_;
  MuonData data_;

};

analyser::analyser(const edm::ParameterSet& iConfig)
{
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");

  muons_ = consumes<View<reco::Muon> >(iConfig.getParameter<InputTag>("muons"));
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  theService_ = new MuonServiceProxy(serviceParameters);

  tree_data_ = data_.book(tree_data_);
}


void
analyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  iSetup.get<MuonGeometryRecord>().get(GEMGeometry_);

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder_);

  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");

  edm::Handle<GEMRecHitCollection> gemRecHits;
  iEvent.getByToken(gemRecHits_, gemRecHits);

  edm::Handle<View<reco::Muon> > muons;
  if (! iEvent.getByToken(muons_, muons)) return;

  if (muons->size() == 0) return;

  for (size_t i = 0; i < muons->size(); ++i){
    cout << "new muon" << endl;
    edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
    const reco::Muon* mu = muRef.get();

    if (mu->pt() < 2.0) continue;
    if (not mu->standAloneMuon()) continue;

    data_.init();

    if (not mu->innerTrack()) continue;
    const reco::Track* innerTrack = mu->track().get();
    reco::TransientTrack ttTrack = ttrackBuilder_->build(innerTrack); //Propagates from inner tracker to GEM

    float count = 0;

    for (const auto& ch : GEMGeometry_->etaPartitions()) {
      if (ch->id().station() != 1) continue; //Only takes GE1/1
      TrajectoryStateOnSurface tsos = propagator->propagate(ttTrack.innermostMeasurementState(),ch->surface());
      if (!tsos.isValid()) continue;
      GlobalPoint pos_global = tsos.globalPosition();
      LocalPoint pos_local = ch->toLocal(tsos.globalPosition());
      const GlobalPoint pos2D_global(pos_global.x(), pos_global.y(), 0);
      const LocalPoint pos2D_local(pos_local.x(), pos_local.y(), 0);
      const BoundPlane& bps(ch->surface());

      if (pos_global.eta() * mu->eta() < 0.0) continue;

      if (bps.bounds().inside(pos2D_local) and ch->id().station() == 1 and ch->id().ring() == 1){
        const float fidcut_angle = 1.0;
        const float cut_ang = 5.0 - fidcut_angle;
        const float fidcut_y = 5.0;
        const float cut_even_high = 250.0 - fidcut_y;
        const float cut_odd_high = 250.0 - fidcut_y;
        const float cut_low = 130.0 + fidcut_y;

        const auto& etaPart_ch = GEMGeometry_->etaPartition(ch->id());
//        float strip = etaPart_ch->strip(pos_local);

        const float prop_y_to_center = etaPart_ch->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).perp();
        LocalPoint local_to_center(pos_local.x(), prop_y_to_center + pos_local.y(), 0);
        const float prop_localphi_rad = (3.14159265/2.) - local_to_center.phi();
        const float prop_localphi_deg = prop_localphi_rad*180/3.14169265;
//        cout << "new ch passed" << endl;
//        cout << "chamber is " << ch->id().chamber() << " and roll is " << ch->id().roll() << " and layer is " << ch->id().layer() << endl;
//        cout << "count is " << count << endl;
        count++;
//        cout << "New angle is " << prop_localphi_deg << endl;

        if (ch->id().layer() == 1){
          data_.has_prop_GE11_L1 = true;
          data_.prop_chamber_GE11_L1 = ch->id().chamber();
          data_.prop_roll_GE11_L1 = ch->id().roll();
          data_.prop_x_GE11_L1 = pos_global.x();
          data_.prop_y_GE11_L1 = pos_global.y();
          data_.prop_r_GE11_L1 = pos_global.mag();
          data_.prop_localx_GE11_L1 = pos_local.x();
          data_.prop_localy_GE11_L1 = pos_local.y();
          data_.prop_y_adjusted_GE11_L1 = prop_y_to_center + pos_local.y();
          data_.prop_localphi_rad_GE11_L1 = prop_localphi_rad;
          data_.prop_localphi_deg_GE11_L1 = prop_localphi_deg;

          if (ch->id().chamber()%2 == 0){
            if (prop_localphi_deg > -cut_ang && prop_localphi_deg < cut_ang && prop_y_to_center + pos_local.y() > cut_low && prop_y_to_center + pos_local.y() < cut_even_high){
              //if(fabs(prop_localphi_deg) < cut_ang) cout << prop_localphi_deg << " is less than " << cut_ang << endl;
              data_.has_fidcut_GE11_L1 = true;
//              cout << "has fidcut = " << data_.has_fidcut_GE11_L1 << " and angle is " << prop_localphi_deg << endl;
            }
          }
          if (ch->id().chamber()%2 == 1){
            if (prop_localphi_deg > -cut_ang && prop_localphi_deg < cut_ang && prop_y_to_center + pos_local.y() > cut_low && prop_y_to_center + pos_local.y() < cut_odd_high){
              //if(fabs(prop_localphi_deg) < cut_ang) cout << prop_localphi_deg << " is less than " << cut_ang << endl;
              data_.has_fidcut_GE11_L1 = true;
//              cout << "has fidcut = " << data_.has_fidcut_GE11_L1 << " and angle is " << prop_localphi_deg << endl;
            }
          }
        }

        if (ch->id().layer() == 2){
          data_.has_prop_GE11_L2 = true;
          data_.prop_chamber_GE11_L2 = ch->id().chamber();
          data_.prop_roll_GE11_L2 = ch->id().roll();
          data_.prop_x_GE11_L2 = pos_global.x();
          data_.prop_y_GE11_L2 = pos_global.y();
          data_.prop_r_GE11_L2 = pos_global.mag();
          data_.prop_localx_GE11_L2 = pos_local.x();
          data_.prop_localy_GE11_L2 = pos_local.y();
          data_.prop_y_adjusted_GE11_L2 = prop_y_to_center + pos_local.y();
          data_.prop_localphi_rad_GE11_L2 = prop_localphi_rad;
          data_.prop_localphi_deg_GE11_L2 = prop_localphi_deg;
          if (ch->id().chamber()%2 == 0){
            if (prop_localphi_deg > -cut_ang && prop_localphi_deg < cut_ang && prop_y_to_center + pos_local.y() > cut_low && prop_y_to_center + pos_local.y() < cut_even_high){
              data_.has_fidcut_GE11_L2 = true;
            }
          }
          if (ch->id().chamber()%2 == 1){
            if (prop_localphi_deg > -cut_ang && prop_localphi_deg < cut_ang && prop_y_to_center + pos_local.y() > cut_low && prop_y_to_center + pos_local.y() < cut_odd_high){
              data_.has_fidcut_GE11_L2 = true;
            }
          }
        }


        for (auto hit = gemRecHits->begin(); hit != gemRecHits->end(); hit++){
          if ( (hit)->geographicalId().det() == DetId::Detector::Muon && (hit)->geographicalId().subdetId() == MuonSubdetId::GEM){
            GEMDetId gemid((hit)->geographicalId());
            if (gemid.chamber() == ch->id().chamber() and gemid.layer() == ch->id().layer() and abs(gemid.roll() - ch->id().roll()) <= 1){
              cout << "starting rechit" << endl;
              const auto& etaPart = GEMGeometry_->etaPartition(gemid);
              float strip = etaPart->strip(hit->localPosition());
              float stripAngle = etaPart->specificTopology().stripAngle(strip);
              float cosAngle = cos(stripAngle);
              float sinAngle = sin(stripAngle);

              float rechit_y_to_center = etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2)).perp();
              LocalPoint local_to_center((hit)->localPosition().x(), rechit_y_to_center + (hit)->localPosition().y(), 0);
              float rechit_localphi_rad = (3.14159265/2.) - local_to_center.phi();
              float rechit_localphi_deg = rechit_localphi_rad*180/3.14159265;
              float deltay_roll =  etaPart_ch->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).perp() - etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2)).perp();


              if (ch->id().station() == 1 and ch->id().ring() == 1 and fabs((hit)->localPosition().x() - pos_local.x()) < 999.0){
 		cout << "Station 1 Ring 1 and dX = " << fabs((hit)->localPosition().x() - pos_local.x()) << endl;
		cout << "RdPhi = " << cosAngle * (pos_local.x() - (hit)->localPosition().x()) + sinAngle * (pos_local.y() + deltay_roll) << endl;
                if (ch->id().layer() == 1){
                  data_.has_rechit_GE11_L1 = true;
                  data_.rechit_chamber_GE11_L1 = gemid.chamber();
                  data_.rechit_roll_GE11_L1 = gemid.roll();
                  data_.rechit_x_GE11_L1 = etaPart->toGlobal((hit)->localPosition()).x();
                  data_.rechit_y_GE11_L1 = etaPart->toGlobal((hit)->localPosition()).y();
                  data_.rechit_r_GE11_L1 = etaPart->toGlobal((hit)->localPosition()).mag();
                  data_.rechit_localx_GE11_L1 = (hit)->localPosition().x();
                  data_.rechit_localy_GE11_L1 = (hit)->localPosition().y();
                  data_.rechit_y_adjusted_GE11_L1 = rechit_y_to_center + (hit)->localPosition().y();
                  data_.rechit_localphi_rad_GE11_L1 = rechit_localphi_rad;
                  data_.rechit_localphi_deg_GE11_L1 = rechit_localphi_deg;
                  data_.RdPhi_GE11_L1 = cosAngle * (pos_local.x() - (hit)->localPosition().x()) + sinAngle * (pos_local.y() + deltay_roll);
                }

                if (ch->id().layer() == 2){
                  data_.has_rechit_GE11_L2 = true;
                  data_.rechit_chamber_GE11_L2 = gemid.chamber();
                  data_.rechit_roll_GE11_L2 = gemid.roll();
                  data_.rechit_x_GE11_L2 = etaPart->toGlobal((hit)->localPosition()).x();
                  data_.rechit_y_GE11_L2 = etaPart->toGlobal((hit)->localPosition()).y();
                  data_.rechit_r_GE11_L2 = etaPart->toGlobal((hit)->localPosition()).mag();
                  data_.rechit_localx_GE11_L2 = (hit)->localPosition().x();
                  data_.rechit_localy_GE11_L2 = (hit)->localPosition().y();
                  data_.rechit_y_adjusted_GE11_L2 = rechit_y_to_center + (hit)->localPosition().y();
                  data_.rechit_localphi_rad_GE11_L2 = rechit_localphi_rad;
                  data_.rechit_localphi_deg_GE11_L2 = rechit_localphi_deg;
                  data_.RdPhi_GE11_L2 = cosAngle * (pos_local.x() - (hit)->localPosition().x()) + sinAngle * (pos_local.y() + deltay_roll);
                }
              }
            }
          }
        }
      }
    }
//    cout << "FILLING. fidcut = " << data_.has_fidcut_GE11_L1 << " angle is = " << data_.prop_localphi_deg_GE11_L1 << endl;
    tree_data_->Fill();
//    cout << "fill tree" << endl;
  }
}

void analyser::beginJob(){}
void analyser::endJob(){}

DEFINE_FWK_MODULE(analyser);
