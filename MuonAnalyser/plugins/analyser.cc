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

  bool has_rechit_GE11;
  int rechit_region_GE11;
  int rechit_layer_GE11;
  int rechit_chamber_GE11;
  int rechit_roll_GE11;
  float rechit_x_GE11;
  float rechit_y_GE11;
  float rechit_r_GE11;
  float rechit_localx_GE11;
  float rechit_localy_GE11;
  float rechit_y_adjusted_GE11;
  float rechit_localphi_rad_GE11;
  float rechit_localphi_deg_GE11;
  int det_id;



};

void MuonData::init()
{
  has_rechit_GE11 = false;
  rechit_region_GE11 = 999999;
  rechit_layer_GE11 = 999999;
  rechit_chamber_GE11 = 999999;
  rechit_roll_GE11 = 999999;
  rechit_x_GE11 = 999999;
  rechit_y_GE11 = 999999;
  rechit_r_GE11 = 999999;
  rechit_localx_GE11 = 999999;
  rechit_localy_GE11 = 999999;
  rechit_y_adjusted_GE11 = 999999;
  rechit_localphi_rad_GE11 = 999999;
  rechit_localphi_deg_GE11 = 999999;
  det_id = 999999;

}

TTree* MuonData::book(TTree *t){
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("MuonData", "MuonData");


//Propogated Inner
//Propogated CSC
//Reconstructed
  t->Branch("has_rechit_GE11", &has_rechit_GE11);
  t->Branch("rechit_region_GE11", &rechit_region_GE11);
  t->Branch("rechit_layer_GE11", &rechit_layer_GE11);
  t->Branch("rechit_chamber_GE11", &rechit_chamber_GE11);
  t->Branch("rechit_roll_GE11", &rechit_roll_GE11);
  t->Branch("rechit_x_GE11", &rechit_x_GE11);
  t->Branch("rechit_y_GE11", &rechit_y_GE11);
  t->Branch("rechit_r_GE11", &rechit_r_GE11);
  t->Branch("rechit_localx_GE11", &rechit_localx_GE11);
  t->Branch("rechit_localy_GE11", &rechit_localy_GE11);
  t->Branch("rechit_y_adjusted_GE11", &rechit_y_adjusted_GE11);
  t->Branch("rechit_localphi_rad_GE11", &rechit_localphi_rad_GE11);
  t->Branch("rechit_localphi_deg_GE11", &rechit_localphi_deg_GE11);
  t->Branch("det_id", &det_id);

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
    const reco::Track* muonTrack = 0;
    if ( mu->globalTrack().isNonnull() ) muonTrack = mu->globalTrack().get();
    else if ( mu->outerTrack().isNonnull()  ) muonTrack = mu->outerTrack().get();
    else 
      continue;
    reco::TransientTrack ttTrack_tracker = ttrackBuilder_->build(innerTrack); //tracker to GEM
    reco::TransientTrack ttTrack_CSC = ttrackBuilder_->build(muonTrack); //CSC to GEM

    for (auto hit = gemRecHits->begin(); hit != gemRecHits->end(); hit++){
      if ( (hit)->geographicalId().det() == DetId::Detector::Muon && (hit)->geographicalId().subdetId() == MuonSubdetId::GEM){
        GEMDetId gemid((hit)->geographicalId());
        cout << "starting rechit" << endl;
        const auto& etaPart = GEMGeometry_->etaPartition(gemid);
        float strip = etaPart->strip(hit->localPosition());
        float rechit_y_to_center = etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2)).perp();
        LocalPoint local_to_center((hit)->localPosition().x(), rechit_y_to_center + (hit)->localPosition().y(), 0);
        float rechit_localphi_rad = (3.14159265/2.) - local_to_center.phi();
        float rechit_localphi_deg = rechit_localphi_rad*180/3.14159265;



        data_.has_rechit_GE11 = true;
        data_.rechit_region_GE11 = gemid.region();
        data_.rechit_layer_GE11 = gemid.layer();
        data_.rechit_chamber_GE11 = gemid.chamber();
        data_.rechit_roll_GE11 = gemid.roll();
        data_.rechit_x_GE11 = etaPart->toGlobal((hit)->localPosition()).x();
        data_.rechit_y_GE11 = etaPart->toGlobal((hit)->localPosition()).y();
        data_.rechit_r_GE11 = etaPart->toGlobal((hit)->localPosition()).mag();
        data_.rechit_localx_GE11 = (hit)->localPosition().x();
        data_.rechit_localy_GE11 = (hit)->localPosition().y();
        data_.rechit_y_adjusted_GE11 = rechit_y_to_center + (hit)->localPosition().y();
        data_.rechit_localphi_rad_GE11 = rechit_localphi_rad;
        data_.rechit_localphi_deg_GE11 = rechit_localphi_deg;
        data_.det_id = gemid.region()*(gemid.station()*100 + gemid.chamber());
      }
    }
    tree_data_->Fill();
  }
}

void analyser::beginJob(){}
void analyser::endJob(){}

DEFINE_FWK_MODULE(analyser);
