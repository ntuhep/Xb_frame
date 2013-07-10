// vim:set ts=4 sw=4 fdm=marker et:
//Update:
//  2012Mar28   pchen       Add GenInfo.
//  2012Mar28   pchen       Fix TrgResult.
//  2012Aug23   pchen       Fix pion mass typo. Fix xxxInfo initial bug.
//  2012Sep30   pchen       Remove opposite charge dipion examination.
//  2012Nov06   pchen       Use official soft muon cut.
//  2013Jan22   pchen       Fix GenInfo signal selection, add PU informations.
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"

        //message logger
#include "FWCore/MessageLogger/interface/MessageLogger.h"
        //magnetic field
#include "MagneticField/Engine/interface/MagneticField.h"
        //vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/Ref.h"
        //For Kalman vertex fitters
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
        //KalmanTrimmedVertexFinder
#include "RecoVertex/TrimmedKalmanVertexFinder/interface/KalmanTrimmedVertexFinder.h"
        //ROOT
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"
        //others
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include <iostream>
        //Xb frame
//#include "XbFrame/Xb_frame/interface/Xb_use.h"
#include "XbFrame/Xb_frame/interface/format.h"
#include "XbFrame/Xb_frame/interface/TriggerBooking.h"

#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TBranch.h"
#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>
#include <math.h>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <list>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#define MUON_MASS   0.10565837
#define PION_MASS   0.13957018
#define JPSI_MASS   3.096916
#define UPS_MASS    9.46030

//
// class declaration
//

class Xb_frame : public edm::EDAnalyzer
{//{{{
    public:
        explicit Xb_frame(const edm::ParameterSet&);
        ~Xb_frame();
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 
    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;
 
        virtual void beginRun(edm::Run const&, edm::EventSetup const&);
        virtual void endRun(edm::Run const&, edm::EventSetup const&);
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
 
        // ----------member data ---------------------------

        edm::Service<TFileService> fs;
        TTree *root;
        EvtInfoBranches     EvtInfo;
        MuonInfoBranches    MuonInfo;
        TrackInfoBranches   TrackInfo;
        XbInfoBranches      XbInfo;
        GenInfoBranches     GenInfo;

        edm::ParameterSet theConfig;
//      std::vector<std::string> TriggersForMatching_;
        edm::InputTag hltLabel_;
        edm::InputTag genLabel_;
        edm::InputTag muonLabel_;
        edm::InputTag trackLabel_;
        edm::InputTag puInfoLabel_;
        std::string   ntupleType_;//'upsilon','jpsi','all','no'

        //histograms
        TH1F *MuonCutLevel;
        TH1F *TrackCutLevel;
        TH1F *XbujCutLevel;
        TH1F *XbMassCutLevel;
};//}}}

void Xb_frame::beginJob()
{//{{{
    root = fs->make<TTree>("root","root");
    EvtInfo.regTree(root);
    MuonInfo.regTree(root);
    TrackInfo.regTree(root);
    XbInfo.regTree(root);
    GenInfo.regTree(root);
}//}}}

Xb_frame::Xb_frame(const edm::ParameterSet& iConfig):theConfig(iConfig)
{//{{{
    //now do what ever initialization is needed
//  TriggersForMatching_= iConfig.getUntrackedParameter<std::vector<std::string> >("TriggersForMatching");
    genLabel_           = iConfig.getParameter<edm::InputTag>("GenLabel");
    trackLabel_         = iConfig.getParameter<edm::InputTag>("TrackLabel");
    muonLabel_          = iConfig.getParameter<edm::InputTag>("MuonLabel");
    hltLabel_           = iConfig.getParameter<edm::InputTag>("HLTLabel");
    puInfoLabel_        = iConfig.getParameter<edm::InputTag>("PUInfoLabel");
    ntupleType_         = iConfig.getUntrackedParameter<std::string>("NtupleType","no_c");

    MuonCutLevel        = fs->make<TH1F>("MuonCutLevel"     , "MuonCutLevel"    , 10, 0, 10);
    TrackCutLevel       = fs->make<TH1F>("TrackCutLevel"    , "TrackCutLevel"   , 10, 0, 10);
    XbujCutLevel        = fs->make<TH1F>("XbujCutLevel"     , "XbujCutLevel"    , 10, 0, 10);
    XbMassCutLevel      = fs->make<TH1F>("XbMassCutLevel"   , "XbMassCutLevel"  , 10, 0, 10);
}//}}}

Xb_frame::~Xb_frame()
{//{{{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}//}}}

//
// member functions
//

// ------------ method called for each event  ------------
void Xb_frame::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //std::cout << "*************************\nReconstructing event number: " << iEvent.id() << "\n";
    using namespace edm;
    using namespace reco;
    ESHandle<MagneticField> bField;
    iSetup.get<IdealMagneticFieldRecord>().get(bField);

    // Change used muon and track collections
    edm::Handle< std::vector<pat::Muon> > muons;
    iEvent.getByLabel(muonLabel_,muons);
    edm::Handle< std::vector<pat::GenericParticle> > tks;
    iEvent.getByLabel(trackLabel_, tks);

    //CLEAN all memory
    memset(&EvtInfo     ,0x00,sizeof(EvtInfo)   );
    memset(&MuonInfo    ,0x00,sizeof(MuonInfo)  );
    memset(&TrackInfo   ,0x00,sizeof(TrackInfo) );
    memset(&XbInfo      ,0x00,sizeof(XbInfo)    );
    memset(&GenInfo     ,0x00,sizeof(GenInfo)   );
    GenInfo.mhmu1_index = -1;
    GenInfo.mhmu2_index = -1;
    GenInfo.mhtk1_index = -1;
    GenInfo.mhtk2_index = -1;
    GenInfo.mhujMass = -1;
    GenInfo.mhxbMass = -1;

    // EvtInfo section{{{
    //EvtInfo.hltnames->clear();
    EvtInfo.RunNo   = iEvent.id().run();
    EvtInfo.EvtNo   = iEvent.id().event();
    EvtInfo.BxNo    = iEvent.bunchCrossing();
    EvtInfo.LumiNo  = iEvent.luminosityBlock();
    EvtInfo.Orbit   = iEvent.orbitNumber();
    EvtInfo.McFlag  = !iEvent.isRealData();
    EvtInfo.nTrgBook= N_TRIGGER_BOOKINGS;
    
    //HLT
    edm::Handle<TriggerResults> TrgResultsHandle; //catch triggerresults
    bool with_TriggerResults = iEvent.getByLabel(hltLabel_,TrgResultsHandle);
    if(!with_TriggerResults){//{{{
        std::cout << "Sorry there is no TriggerResult in the file" << std::endl;
    }else{
        //get the names of the triggers
        const edm::TriggerNames &TrgNames = iEvent.triggerNames(*TrgResultsHandle);
        EvtInfo.trgCount = 0;
        for(int i=0; i< N_TRIGGER_BOOKINGS; i++){
            unsigned int TrgIndex = TrgNames.triggerIndex(TriggerBooking[i]);
            if (TrgIndex == TrgNames.size()) {
                EvtInfo.trgBook[i] = -4; // The trigger path is not known in this event.
            }else if ( !TrgResultsHandle->wasrun( TrgIndex ) ) {
                EvtInfo.trgBook[i] = -3; // The trigger path was not included in this event.
            }else if ( !TrgResultsHandle->accept( TrgIndex ) ) {
                EvtInfo.trgBook[i] = -2; // The trigger path was not accepted in this event.
            }else if (  TrgResultsHandle->error ( TrgIndex ) ) {
                EvtInfo.trgBook[i] = -1; // The trigger path has an error in this event.
            }else {
                EvtInfo.trgBook[i] = +1; // It's triggered.
                EvtInfo.trgCount++; 
            }
        }
        EvtInfo.nHLT = TrgNames.size();
        for(unsigned int i=0; i<TrgNames.size(); i++){
            EvtInfo.hltBits[i] = (TrgResultsHandle->accept(i) == true) ? 1:0;
        }
    }//end(!with_TriggerResults)}}}

    // Handle primary vertex properties
    Vertex thePrimaryV;
    math::XYZPoint RefVtx;
        //get beamspot information
    Vertex theBeamSpotV;
    reco::BeamSpot beamSpot;
    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
    if (beamSpotHandle.isValid()){
        beamSpot = *beamSpotHandle;
        theBeamSpotV = Vertex(beamSpot.position(), beamSpot.covariance3D());
    }else{
        std::cout<< "No beam spot available from EventSetup \n";
    }

        //get vertex informationa
    edm::Handle<reco::VertexCollection> VertexHandle;
    iEvent.getByLabel("offlinePrimaryVertexHandle", VertexHandle);
    if (!VertexHandle.failedToGet() && VertexHandle->size()>0){
        //int nVtxTrks = 0;//outdated PV definition
        double max_tkSt = 0;
        for(std::vector<reco::Vertex>::const_iterator it_vtx = VertexHandle->begin();
            it_vtx != VertexHandle->end(); it_vtx++){
            if (!it_vtx->isValid()) continue;
            ////find primary vertex with most number of associated tracks, outdated PV defination
            //int tksize = it_vtx->tracksSize();
            //if(nVtxTrks < tksize){
            //    nVtxTrks = tksize;
            //    thePrimaryV = Vertex(*it_vtx);
            //}

            //find primary primary vertex with largest St
            double tkSt = 0;
            for(std::vector<reco::TrackBaseRef>::const_iterator it_tk = it_vtx->tracks_begin();
                it_tk != it_vtx->tracks_end(); it_tk++){
                tkSt += it_tk->get()->pt();
            }
            if (tkSt > max_tkSt){
                max_tkSt = tkSt;
                thePrimaryV = Vertex(*it_vtx);
            }
        }
    }else{ 
        thePrimaryV = Vertex(beamSpot.position(), beamSpot.covariance3D());
    }
    RefVtx = thePrimaryV.position();

    EvtInfo.PVx     = thePrimaryV.position().x();
    EvtInfo.PVy     = thePrimaryV.position().y();
    EvtInfo.PVz     = thePrimaryV.position().z();
    EvtInfo.PVxE    = thePrimaryV.xError();
    EvtInfo.PVyE    = thePrimaryV.yError();
    EvtInfo.PVzE    = thePrimaryV.zError();
    EvtInfo.PVnchi2 = thePrimaryV.normalizedChi2();
    EvtInfo.PVchi2  = thePrimaryV.chi2();

        // get pile-up information
    std::cout << "puInfoLabel_=" << puInfoLabel_ << std::endl;
    if (!iEvent.isRealData()){
        edm::Handle<std::vector< PileupSummaryInfo > >  PUHandle;
        iEvent.getByLabel(puInfoLabel_, PUHandle);
        std::vector<PileupSummaryInfo>::const_iterator PVI;
        for(PVI = PUHandle->begin(); PVI != PUHandle->end(); ++PVI) {
            EvtInfo.nPU[EvtInfo.nBX]   = PVI->getPU_NumInteractions();
            EvtInfo.BXPU[EvtInfo.nBX]  = PVI->getBunchCrossing();
            EvtInfo.trueIT[EvtInfo.nBX]= PVI->getTrueNumInteractions();
            EvtInfo.nBX += 1;
        }
    }else{
        EvtInfo.nBX = 0;
    }

    //}}}
    //printf("-----*****DEBUG:End of EvtInfo.\n");

    // Double check size=0.
    MuonInfo.size   = 0;
    TrackInfo.size  = 0;
    XbInfo.uj_size  = 0;
    XbInfo.size     = 0;
    GenInfo.size    = 0;

    std::vector<pat::Muon>              input_muons;
    std::vector<pat::GenericParticle>   input_tracks;
    input_muons = *muons;
    input_tracks = *tks;
    try{
        //standard check for validity of input data
        if (input_muons.size() == 0){
            std::cout << "There's no muon : " << iEvent.id() << std::endl;
        }else{
            std::cout << "Got " << input_muons.size() << " muons" << std::endl;
            if (input_tracks.size() == 0){
                std::cout << "There's no track: " << iEvent.id() << std::endl;
            }else{
                std::cout << "Got " << input_tracks.size() << " tracks" << std::endl;
                if (input_tracks.size() > 1 && input_muons.size() > 1){

                    //MuonInfo section{{{
                    int mu_hindex = -1;
                    const reco::GenParticle* genMuonPtr[MAX_MUON];
                    memset(genMuonPtr,0x00,MAX_MUON);
                    for(std::vector<pat::Muon>::const_iterator mu_it=input_muons.begin();
                        mu_it != input_muons.end() ; mu_it++){
                        if(mu_hindex >= MAX_MUON){
                            fprintf(stderr,"ERROR: number of muons exceeds the size of array.\n");
                            break;//exit(0);
                        }
                        mu_hindex = int(mu_it - input_muons.begin());
                        MuonCutLevel->Fill(0);
                        //if(!(mu_it->innerTrack().isNonnull()*mu_it->globalTrack().isNonnull())) {continue;}
                        //MuonCutLevel->Fill(1);
                        if (!(mu_it->isTrackerMuon() || mu_it->isGlobalMuon())) continue;
                        MuonCutLevel->Fill(2);
                        //if (!(mu_it->isGlobalMuon()*mu_it->track().isAvailable()*mu_it->globalTrack().isAvailable())) continue;
                        //    MuonCutLevel->Fill(3);
                        if (mu_it->p()>200 || mu_it->pt()>200)                  continue;
                        MuonCutLevel->Fill(4);
                        if (!muon::isGoodMuon(*mu_it,muon::TMOneStationTight))  continue;
                        MuonCutLevel->Fill(5);
                        if (fabs(mu_it->innerTrack()->dxy(thePrimaryV.position())) >= 3.        || 
                            fabs(mu_it->innerTrack()->dz(thePrimaryV.position()))  >= 30.       
                           ) continue;
                        if (mu_it->innerTrack()->hitPattern().pixelLayersWithMeasurement()<2    ||
                            mu_it->innerTrack()->normalizedChi2()>1.8                           //||
                           ) continue;
                        if (mu_it->innerTrack()->hitPattern().trackerLayersWithMeasurement()<6  &&
                            mu_it->innerTrack()->hitPattern().numberOfValidStripHits()<11
                           ) continue;
                        MuonCutLevel->Fill(6);

                        MuonInfo.index          [MuonInfo.size] = MuonInfo.size;
                        MuonInfo.handle_index   [MuonInfo.size] = mu_hindex;
//                      MuonInfo.position       [MuonInfo.size] = distance(input_muons.begin(),mu_it);
//                      std::cout<<"distance: "<<std::distance(input_muons.begin(),mu_it)<<endl;
                        MuonInfo.charge         [MuonInfo.size] = mu_it->charge();
                        MuonInfo.pt             [MuonInfo.size] = mu_it->pt();
                        MuonInfo.eta            [MuonInfo.size] = mu_it->eta();
                        MuonInfo.phi            [MuonInfo.size] = mu_it->phi();
//                      MuonInfo.p              [MuonInfo.size] = mu_it->p();
                        MuonInfo.i_striphit     [MuonInfo.size] = mu_it->innerTrack()->hitPattern().numberOfValidStripHits();
                        MuonInfo.i_pixelhit     [MuonInfo.size] = mu_it->innerTrack()->hitPattern().numberOfValidPixelHits();
                        MuonInfo.i_nStripLayer  [MuonInfo.size] = mu_it->innerTrack()->hitPattern().stripLayersWithMeasurement();
                        MuonInfo.i_nPixelLayer  [MuonInfo.size] = mu_it->innerTrack()->hitPattern().pixelLayersWithMeasurement();
                        MuonInfo.i_chi2         [MuonInfo.size] = mu_it->innerTrack()->chi2();
                        MuonInfo.i_ndf          [MuonInfo.size] = mu_it->innerTrack()->ndof();
                        MuonInfo.fpbarrelhit    [MuonInfo.size] = mu_it->innerTrack()->hitPattern().hasValidHitInFirstPixelBarrel();
                        MuonInfo.fpendcaphit    [MuonInfo.size] = mu_it->innerTrack()->hitPattern().hasValidHitInFirstPixelEndcap();
                        MuonInfo.d0             [MuonInfo.size] = mu_it->track()->d0();
                        MuonInfo.dz             [MuonInfo.size] = mu_it->track()->dz();
                        MuonInfo.dzPV           [MuonInfo.size] = mu_it->track()->dz(RefVtx);
                        MuonInfo.dxyPV          [MuonInfo.size] = mu_it->track()->dxy(RefVtx);
                        MuonInfo.iso_trk        [MuonInfo.size] = mu_it->trackIso();//R<0.3
                        //MuonInfo.iso_calo       [MuonInfo.size] = mu_it->caloIso();//sum of two iso
                        MuonInfo.iso_ecal       [MuonInfo.size] = mu_it->ecalIso();
                        MuonInfo.iso_hcal       [MuonInfo.size] = mu_it->hcalIso();
                        MuonInfo.n_matches      [MuonInfo.size] = mu_it->numberOfMatches();//only in chamber

                        if (!iEvent.isRealData())
                            genMuonPtr [MuonInfo.size] = mu_it->genParticle();

                        if(mu_it->isGlobalMuon()){
                            MuonInfo.g_striphit [MuonInfo.size] = mu_it->globalTrack()->hitPattern().numberOfValidStripHits();
                            MuonInfo.g_pixelhit [MuonInfo.size] = mu_it->globalTrack()->hitPattern().numberOfValidPixelHits();
                            MuonInfo.g_chi2     [MuonInfo.size] = mu_it->globalTrack()->chi2();
                            MuonInfo.g_ndf      [MuonInfo.size] = mu_it->globalTrack()->ndof();
                            MuonInfo.nmuhit     [MuonInfo.size] = mu_it->globalTrack()->hitPattern().numberOfValidMuonHits();
                        }else{
                            MuonInfo.g_striphit [MuonInfo.size] = -1;
                            MuonInfo.g_pixelhit [MuonInfo.size] = -1;
                            MuonInfo.g_chi2     [MuonInfo.size] = -1;
                            MuonInfo.g_ndf      [MuonInfo.size] = -1;
                            MuonInfo.nmuhit     [MuonInfo.size] = -1;
                        }

                        int qm = 0;
                        for(int qi=1; qi!= 24; ++qi){
                            if (muon::isGoodMuon(*mu_it, muon::SelectionType(qi))){
                                qm += 1 << qi;
                            }
                        }
                        MuonInfo.muqual         [MuonInfo.size] = qm;   

                        //2.tracker hits requirement for MC/2011/2012 soft pion
                        //1.make sure non-zero for binary test
                        //0.pass all cuts
                        if (!iEvent.isRealData() || 
                            (iEvent.id().run() < 180297 && mu_it->innerTrack()->hitPattern().numberOfValidStripHits() > 10) || 
                            (iEvent.id().run() > 180296 && mu_it->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 )
                           ) MuonInfo.isGoodCand[MuonInfo.size] += 1 << 2;
                        MuonInfo.isGoodCand[MuonInfo.size] += 1 << 1;
                        if (((MuonInfo.isGoodCand[MuonInfo.size]>>1)&((MuonInfo.isGoodCand[MuonInfo.size]>>1)+1)) == 0)
                            MuonInfo.isGoodCand[MuonInfo.size] += 1;

                        MuonInfo.size++;
                    }//end of MuonInfo}}}
                    //printf("-----*****DEBUG:End of MuonInfo.\n");

                    //Preselect tracks{{{
                    bool isNeededTrack[MAX_TRACK];// Are the tracks redundant?
                    memset(isNeededTrack,false,MAX_TRACK);
                    for(std::vector<pat::GenericParticle>::const_iterator tk_it=input_tracks.begin();
                        tk_it != input_tracks.end(); tk_it++){
                        TrackCutLevel->Fill(0);//number of all tracks
                        bool isMuonTrack = false; //remove muon track
                        for(std::vector<pat::Muon>::iterator it=input_muons.begin() ; it != input_muons.end() ; it++){
                            if (!it->track().isNonnull())                   continue;
                            if (fabs(tk_it->pt() -it->track()->pt() )<0.00001 &&
                                fabs(tk_it->eta()-it->track()->eta())<0.00001 &&
                                fabs(tk_it->phi()-it->track()->phi())<0.00001 ){
                                    isMuonTrack = true;
                                    break;
                            }
                        }
                        if (isMuonTrack)                                    continue;
                        TrackCutLevel->Fill(1);//number of non muon tracks
                        if (tk_it->track()->normalizedChi2()>5)             continue;
                        TrackCutLevel->Fill(2);
                        if (tk_it->pt()<0.4)                                continue;
                        TrackCutLevel->Fill(3);
                        if (tk_it->p()>200 || tk_it->pt()>200)              continue;
                        TrackCutLevel->Fill(4);
                        if (fabs(tk_it->eta()) > 2.5)                       continue;
                        TrackCutLevel->Fill(5);
                        if (tk_it->track()->hitPattern().numberOfValidStripHits()<10)continue;
                        TrackCutLevel->Fill(6);
                        if (tk_it->track()->hitPattern().numberOfValidPixelHits()<2) continue;
                        TrackCutLevel->Fill(7);
                        isNeededTrack[tk_it-input_tracks.begin()] = true;
                    }//end of track preselection}}}
                    //printf("-----*****DEBUG:End of track preselection.\n");

                    // XbInfo section{{{
                    int mu1_index = -1;
                    int mu1_hindex = -1;
                    bool gogogo = false;
                    for(std::vector<pat::Muon>::const_iterator mu_it1=input_muons.begin();
                        mu_it1 != input_muons.end() ; mu_it1++){
                        //Check if it in MuonInfo
                        mu1_hindex = int(mu_it1 - input_muons.begin());
                        gogogo = false;
                        for(int i=0; i < MuonInfo.size; i++){
                            if (mu1_hindex == MuonInfo.handle_index[i]){
                                gogogo = true;
                                break;
                            }
                        }
                        if (!gogogo) continue;
                        //Get the corrisponding index in MuonInfo
                        mu1_index ++;
                        int mu2_index = mu1_index;
                        int mu2_hindex = -1; 
                        for(std::vector<pat::Muon>::const_iterator mu_it2=mu_it1+1;
                            mu_it2 != input_muons.end() ; mu_it2++){
                            mu2_hindex = int(mu_it2 - input_muons.begin()); 
                            gogogo = false;
                            for(int j=0; j < MuonInfo.size; j++){
                                if(mu2_hindex == MuonInfo.handle_index[j]){
                                    gogogo = true;
                                    break;
                                }
                            }
                            if (!gogogo) continue;
                            mu2_index ++;   
                            if (mu_it1->charge() * mu_it2->charge()>0) continue;
                            XbujCutLevel->Fill(0);
    
                            //Fit 2 muon
                            reco::TransientTrack muonPTT(mu_it1->track(), &(*bField) );
                            reco::TransientTrack muonMTT(mu_it2->track(), &(*bField) );
                            if(!muonPTT.isValid()) continue;
                            if(!muonMTT.isValid()) continue;
                            XbujCutLevel->Fill(1);
    
                            const reco::Muon* rmu1 = dynamic_cast<const reco::Muon * >(mu_it1->originalObject());
                            const reco::Muon* rmu2 = dynamic_cast<const reco::Muon * >(mu_it2->originalObject());
                            if(muon::overlap(*rmu1, *rmu2)) continue;
                            XbujCutLevel->Fill(2);
    
                            KinematicParticleFactoryFromTransientTrack pFactory;
                            ParticleMass muon_mass = MUON_MASS; //pdg mass
                            float muon_sigma = muon_mass*1.e-6;
                            float chi = 0.;
                            float ndf = 0.;
                            std::vector<RefCountedKinematicParticle> muonParticles;
                            muonParticles.push_back(pFactory.particle(muonPTT,muon_mass,chi,ndf,muon_sigma));
                            muonParticles.push_back(pFactory.particle(muonMTT,muon_mass,chi,ndf,muon_sigma));
    
                            KinematicParticleVertexFitter   fitter;   
                            RefCountedKinematicTree         ujVFT;
                            ujVFT = fitter.fit(muonParticles); 
                            if (!ujVFT->isValid()) continue;
                            XbujCutLevel->Fill(3); 
                            ujVFT->movePointerToTheTop();
    
                            RefCountedKinematicParticle ujVFP       = ujVFT->currentParticle();
                            RefCountedKinematicVertex   ujVFPvtx    = ujVFT->currentDecayVertex();
                            ujVFT->movePointerToTheFirstChild();
                            KinematicParameters         ujmu1KP     = ujVFT->currentParticle()->currentState().kinematicParameters();
                            ujVFT->movePointerToTheNextChild();
                            KinematicParameters         ujmu2KP     = ujVFT->currentParticle()->currentState().kinematicParameters();
                            double chi2_prob_uj = TMath::Prob(ujVFPvtx->chiSquared(), ujVFPvtx->degreesOfFreedom());
                            if(chi2_prob_uj < 0.01) continue;
                            XbujCutLevel->Fill(4);

                            double dimuon_mass = ujVFP->currentState().mass();
                            if ((dimuon_mass < 11.5 && dimuon_mass > 8.5) && ntupleType_.compare("jpsi") != 0){
                                XbujCutLevel->Fill(5);
                            }else if ((dimuon_mass < 6. && dimuon_mass > 3.) && ntupleType_.compare("upsilon") != 0){
                                XbujCutLevel->Fill(6);
                            }else{
                                continue;
                            }


                            TLorentzVector uj_4vec,uj_mu1_4vec,uj_mu2_4vec;
                            uj_4vec.SetPxPyPzE(ujVFP->currentState().kinematicParameters().momentum().x(),
                                               ujVFP->currentState().kinematicParameters().momentum().y(),
                                               ujVFP->currentState().kinematicParameters().momentum().z(),
                                               ujVFP->currentState().kinematicParameters().energy());
                            uj_mu1_4vec.SetPxPyPzE( ujmu1KP.momentum().x(),
                                                    ujmu1KP.momentum().y(),
                                                    ujmu1KP.momentum().z(),
                                                    ujmu1KP.energy());
                            uj_mu2_4vec.SetPxPyPzE( ujmu2KP.momentum().x(),
                                                    ujmu2KP.momentum().y(),
                                                    ujmu2KP.momentum().z(),
                                                    ujmu2KP.energy());
                            //uj_4vec.Print();
                            //uj_mu1_4vec.Print();
                            //uj_mu2_4vec.Print();

                            XbInfo.uj_index         [XbInfo.uj_size]= XbInfo.uj_size;
                            XbInfo.uj_mass          [XbInfo.uj_size]= uj_4vec.Mag();
                            XbInfo.uj_px            [XbInfo.uj_size]= uj_4vec.Px();
                            XbInfo.uj_py            [XbInfo.uj_size]= uj_4vec.Py();
                            XbInfo.uj_pz            [XbInfo.uj_size]= uj_4vec.Pz();
                            XbInfo.uj_vtxX          [XbInfo.uj_size]= ujVFPvtx->position().x();
                            XbInfo.uj_vtxY          [XbInfo.uj_size]= ujVFPvtx->position().y();
                            XbInfo.uj_vtxZ          [XbInfo.uj_size]= ujVFPvtx->position().z();
                            XbInfo.uj_vtxXE         [XbInfo.uj_size]= sqrt(ujVFPvtx->error().cxx());
                            XbInfo.uj_vtxYE         [XbInfo.uj_size]= sqrt(ujVFPvtx->error().cyy());
                            XbInfo.uj_vtxZE         [XbInfo.uj_size]= sqrt(ujVFPvtx->error().czz());
                            XbInfo.uj_vtxdof        [XbInfo.uj_size]= ujVFPvtx->degreesOfFreedom();
                            XbInfo.uj_vtxchi2       [XbInfo.uj_size]= ujVFPvtx->chiSquared();
                            XbInfo.uj_rfmu1_index   [XbInfo.uj_size]= mu1_index;
                            XbInfo.uj_rfmu2_index   [XbInfo.uj_size]= mu2_index;
                            XbInfo.uj_rfmu1_px      [XbInfo.uj_size]= uj_mu1_4vec.Px();
                            XbInfo.uj_rfmu1_py      [XbInfo.uj_size]= uj_mu1_4vec.Py();
                            XbInfo.uj_rfmu1_pz      [XbInfo.uj_size]= uj_mu1_4vec.Pz();
                            //XbInfo.uj_rfmu1_e       [XbInfo.uj_size]= uj_mu1_4vec.E()
                            XbInfo.uj_rfmu2_px      [XbInfo.uj_size]= uj_mu2_4vec.Px();
                            XbInfo.uj_rfmu2_py      [XbInfo.uj_size]= uj_mu2_4vec.Py();
                            XbInfo.uj_rfmu2_pz      [XbInfo.uj_size]= uj_mu2_4vec.Pz();
                            //XbInfo.uj_rfmu2_e       [XbInfo.uj_size]= uj_mu2_4vec.E();

                            //4.|uj_Eta| < 1.2?
                            //3.good muon candidates?
                            //2.make sure non-zero for binary test
                            //1.Xb or X3872
                            //0.passed all selections
                            if (fabs(uj_4vec.Eta()) < 1.2)//Eta < 1.2
                                XbInfo.uj_isGoodCand[XbInfo.uj_size] += 1 << 4;
                            if (MuonInfo.isGoodCand[mu1_index]%2 == 1 && MuonInfo.isGoodCand[mu2_index]%2 == 1)
                                XbInfo.uj_isGoodCand[XbInfo.uj_size] += 1 << 3;
                            XbInfo.uj_isGoodCand[XbInfo.uj_size] += 1 << 2;
                            if (uj_4vec.Mag() > 7)//Xb or X3872
                                XbInfo.uj_isGoodCand[XbInfo.uj_size] += 1 << 1;
                            if (((XbInfo.uj_isGoodCand[XbInfo.uj_size]>>2)&((XbInfo.uj_isGoodCand[XbInfo.uj_size]>>2)+1))==0)
                                XbInfo.uj_isGoodCand[XbInfo.uj_size] += 1;

                            XbInfo.uj_size++;
                            muonParticles.clear();

                            int tk1_hindex = -1;
                            for(std::vector<pat::GenericParticle>::const_iterator tk_it1=input_tracks.begin();
                                tk_it1 != input_tracks.end() ; tk_it1++){
                                tk1_hindex = int(tk_it1 - input_tracks.begin());
                                if (!isNeededTrack[tk1_hindex]) continue;
                                if (abs(tk_it1->charge()) != 1) continue;
    
                                int tk2_hindex = -1;
                                for(std::vector<pat::GenericParticle>::const_iterator tk_it2=tk_it1+1;
                                    tk_it2 != input_tracks.end() ; tk_it2++){
                                    tk2_hindex = int(tk_it2 - input_tracks.begin()); 
                                    if (!isNeededTrack[tk2_hindex]) continue;
            
                                    if (abs(tk_it2->charge()) != 1) continue;
                                    //if (tk_it1->charge() * tk_it2->charge()>0) continue;// ignore charge test
            
                                    reco::TransientTrack pionPTT(tk_it1->track(), &(*bField) );
                                    reco::TransientTrack pionMTT(tk_it2->track(), &(*bField) );
                                    if (!pionPTT.isValid()) continue;
                                    if (!pionMTT.isValid()) continue;
            
                                    ParticleMass pion_mass = PION_MASS;
                                    float pion_sigma = pion_mass*1.e-6;
            
                                    if (XbInfo.size >= MAX_XB) continue;
                                    //doing tktk fit
                                    std::vector<RefCountedKinematicParticle> tktk_candidate;
                                    tktk_candidate.push_back(pFactory.particle(pionPTT,pion_mass,chi,ndf,pion_sigma));
                                    tktk_candidate.push_back(pFactory.particle(pionMTT,pion_mass,chi,ndf,pion_sigma));
        
                                    XbMassCutLevel->Fill(0);
                                    KinematicParticleVertexFitter   tktk_fitter;
                                    RefCountedKinematicTree         tktk_VFT;
                                    tktk_VFT = tktk_fitter.fit(tktk_candidate);
                                    if(!tktk_VFT->isValid()) continue;
                                    tktk_VFT->movePointerToTheTop();
                                    XbMassCutLevel->Fill(1);
                                    RefCountedKinematicVertex tktk_VFPvtx = tktk_VFT->currentDecayVertex();
                                    double chi2_prob_tktk = TMath::Prob(tktk_VFPvtx->chiSquared(),
                                                                        tktk_VFPvtx->degreesOfFreedom());
                                    if(chi2_prob_tktk < 0.1) continue;
                                    XbMassCutLevel->Fill(2);
            
                                    std::vector<RefCountedKinematicParticle> Xb_candidate;
                                    Xb_candidate.push_back(pFactory.particle(muonPTT,muon_mass,chi,ndf,muon_sigma));
                                    Xb_candidate.push_back(pFactory.particle(muonMTT,muon_mass,chi,ndf,muon_sigma));
                                    Xb_candidate.push_back(pFactory.particle(pionPTT,pion_mass,chi,ndf,pion_sigma));
                                    Xb_candidate.push_back(pFactory.particle(pionMTT,pion_mass,chi,ndf,pion_sigma));
                                    RefCountedKinematicTree xbVFT;
        
                                    if (ntupleType_.compare("no_c") != 0){ //refit on uj mass
                                        ParticleMass uj_mass;
                                        if (dimuon_mass > 7 ){
                                            uj_mass = UPS_MASS;
                                        }else{
                                            uj_mass = JPSI_MASS;
                                        }
                                        MultiTrackKinematicConstraint *uj_c = new  TwoTrackMassKinematicConstraint(uj_mass);
                                        KinematicConstrainedVertexFitter kcvFitter;
                                        xbVFT = kcvFitter.fit(Xb_candidate, uj_c);
                                    }else{   //only 4 track refit
                                        KinematicParticleVertexFitter kcvFitter;
                                        xbVFT = kcvFitter.fit(Xb_candidate);
                                    }
            
                                    if (!xbVFT->isValid()) continue;
                                    xbVFT->movePointerToTheTop();
                                    RefCountedKinematicParticle     xbVFP       = xbVFT->currentParticle();
                                    RefCountedKinematicVertex       xbVFPvtx    = xbVFT->currentDecayVertex();
                                    if (!xbVFPvtx->vertexIsValid()) continue;
                                    std::vector<RefCountedKinematicParticle> xCands  = xbVFT->finalStateParticles();
            
                                    double chi2_prob = TMath::Prob(xbVFPvtx->chiSquared(),xbVFPvtx->degreesOfFreedom());
                                    if (chi2_prob < 0.1) continue;
                                    XbMassCutLevel->Fill(3);

                                    //Cut out a mass window
                                    if (!(xbVFP->currentState().mass() < 12. && xbVFP->currentState().mass() > 9.) &&
                                        !(xbVFP->currentState().mass() <  6. && xbVFP->currentState().mass() > 3.) ){
                                        continue;
                                    }
                                    XbMassCutLevel->Fill(4);

                                    TLorentzVector xb_4vec,xb_mu1_4vec,xb_mu2_4vec,xb_tk1_4vec,xb_tk2_4vec;
                                    xb_4vec.SetPxPyPzE(xbVFP->currentState().kinematicParameters().momentum().x(),
                                                       xbVFP->currentState().kinematicParameters().momentum().y(),
                                                       xbVFP->currentState().kinematicParameters().momentum().z(),
                                                       xbVFP->currentState().kinematicParameters().energy());
                                    xb_mu1_4vec.SetPxPyPzE(xCands[0]->currentState().kinematicParameters().momentum().x(),
                                                           xCands[0]->currentState().kinematicParameters().momentum().y(),
                                                           xCands[0]->currentState().kinematicParameters().momentum().z(),
                                                           xCands[0]->currentState().kinematicParameters().energy());
                                    xb_mu2_4vec.SetPxPyPzE(xCands[1]->currentState().kinematicParameters().momentum().x(),
                                                           xCands[1]->currentState().kinematicParameters().momentum().y(),
                                                           xCands[1]->currentState().kinematicParameters().momentum().z(),
                                                           xCands[1]->currentState().kinematicParameters().energy());
                                    xb_tk1_4vec.SetPxPyPzE(xCands[2]->currentState().kinematicParameters().momentum().x(),
                                                           xCands[2]->currentState().kinematicParameters().momentum().y(),
                                                           xCands[2]->currentState().kinematicParameters().momentum().z(),
                                                           xCands[2]->currentState().kinematicParameters().energy());
                                    xb_tk2_4vec.SetPxPyPzE(xCands[3]->currentState().kinematicParameters().momentum().x(),
                                                           xCands[3]->currentState().kinematicParameters().momentum().y(),
                                                           xCands[3]->currentState().kinematicParameters().momentum().z(),
                                                           xCands[3]->currentState().kinematicParameters().energy());
                                    //xb_4vec.Print();
                                    //xb_mu1_4vec.Print();
                                    //xb_mu2_4vec.Print();
                                    //xb_tk1_4vec.Print();
                                    //xb_tk2_4vec.Print();
            
                                    XbInfo.index[XbInfo.size]   = XbInfo.size;
                                    XbInfo.mass[XbInfo.size]    = xb_4vec.Mag();
                                    //std::cout << "xbCandMass[" << XbInfo.size << "]=" << xb_4vec.Mag() << std::endl;
                                    XbInfo.px[XbInfo.size]      = xb_4vec.Px();
                                    XbInfo.py[XbInfo.size]      = xb_4vec.Py();
                                    XbInfo.pz[XbInfo.size]      = xb_4vec.Pz();
                                    XbInfo.pxE[XbInfo.size]     = sqrt(xbVFP->currentState().kinematicParametersError().matrix()(3,3));
                                    XbInfo.pyE[XbInfo.size]     = sqrt(xbVFP->currentState().kinematicParametersError().matrix()(4,4));
                                    XbInfo.pzE[XbInfo.size]     = sqrt(xbVFP->currentState().kinematicParametersError().matrix()(5,5));
                                    XbInfo.vtxX[XbInfo.size]    = xbVFPvtx->position().x();
                                    XbInfo.vtxY[XbInfo.size]    = xbVFPvtx->position().y();
                                    XbInfo.vtxZ[XbInfo.size]    = xbVFPvtx->position().z();
                                    XbInfo.vtxXE[XbInfo.size]   = sqrt(xbVFPvtx->error().cxx());
                                    XbInfo.vtxYE[XbInfo.size]   = sqrt(xbVFPvtx->error().cyy());
                                    XbInfo.vtxZE[XbInfo.size]   = sqrt(xbVFPvtx->error().czz());
                                    XbInfo.vtxdof[XbInfo.size]  = xbVFPvtx->degreesOfFreedom();
                                    XbInfo.vtxchi2[XbInfo.size] = xbVFPvtx->chiSquared();
    
                                    XbInfo.rfuj_index[XbInfo.size]  = XbInfo.uj_size-1;
                                    //XbInfo.rfmu1_index[XbInfo.size] = mu1_index; //Not used since format_1_7
                                    //XbInfo.rfmu2_index[XbInfo.size] = mu2_index; 
                                    XbInfo.rftk1_index[XbInfo.size] = -tk1_hindex-1;//For later matching in track section.
                                    XbInfo.rftk2_index[XbInfo.size] = -tk2_hindex-1;

                                    XbInfo.rfmu1_px[XbInfo.size]=xb_mu1_4vec.Px();
                                    XbInfo.rfmu1_py[XbInfo.size]=xb_mu1_4vec.Py();
                                    XbInfo.rfmu1_pz[XbInfo.size]=xb_mu1_4vec.Pz(); 
                                    //XbInfo.rfmu1_e [XbInfo.size]=xb_mu1_4vec.E(); 
                                    XbInfo.rfmu2_px[XbInfo.size]=xb_mu2_4vec.Px(); 
                                    XbInfo.rfmu2_py[XbInfo.size]=xb_mu2_4vec.Py(); 
                                    XbInfo.rfmu2_pz[XbInfo.size]=xb_mu2_4vec.Pz(); 
                                    //XbInfo.rfmu2_e [XbInfo.size]=xb_mu2_4vec.E();  
                                    XbInfo.rftk1_px[XbInfo.size]=xb_tk1_4vec.Px(); 
                                    XbInfo.rftk1_py[XbInfo.size]=xb_tk1_4vec.Py(); 
                                    XbInfo.rftk1_pz[XbInfo.size]=xb_tk1_4vec.Pz(); 
                                    //XbInfo.rftk1_e [XbInfo.size]=xb_tk1_4vec.E();  
                                    XbInfo.rftk2_px[XbInfo.size]=xb_tk2_4vec.Px(); 
                                    XbInfo.rftk2_py[XbInfo.size]=xb_tk2_4vec.Py(); 
                                    XbInfo.rftk2_pz[XbInfo.size]=xb_tk2_4vec.Pz(); 
                                    //XbInfo.rftk2_e [XbInfo.size]=xb_tk2_4vec.E();  
                                    
                                    //9.dR(ujPi2) < 0.8
                                    //8.dR(ujPi1) < 0.8
                                    //7.uj prob > 2%
                                    //6.uj mass window.
                                    //5.opposite charge sign?
                                    //4.good uj candidates?
                                    //3.good tracks?(Will be assigned in track section)
                                    //2.make sure non-zero for further binary test of all standard cut.
                                    //1.jpsi or upsilon?
                                    //0.Passed all standard cut?
                                    if (reco::deltaR(uj_4vec.Eta(),uj_4vec.Phi(),xb_tk2_4vec.Eta(),xb_tk2_4vec.Phi())< 0.8)
                                        XbInfo.isGoodCand[XbInfo.size] += 1 << 9;
                                    if (reco::deltaR(uj_4vec.Eta(),uj_4vec.Phi(),xb_tk1_4vec.Eta(),xb_tk1_4vec.Phi())< 0.8 )
                                        XbInfo.isGoodCand[XbInfo.size] += 1 << 8;
                                    if (chi2_prob_uj > 0.02 )
                                        XbInfo.isGoodCand[XbInfo.size] += 1 << 7;
                                    if ((uj_4vec.Mag() > 9.300 && uj_4vec.Mag()< 9.620) ||
                                        (uj_4vec.Mag() > 3.000 && uj_4vec.Mag()< 3.194) )
                                        XbInfo.isGoodCand[XbInfo.size] += 1 << 6;
                                    if (tk_it1->charge()*tk_it2->charge()<0)
                                        XbInfo.isGoodCand[XbInfo.size] += 1 << 5;
                                    if ((XbInfo.uj_isGoodCand[XbInfo.uj_size]&1) == 1)
                                        XbInfo.isGoodCand[XbInfo.size] += 1 << 4;
                                    ////if ((TrackInfo.isGoodCand[XbInfo.rftk1_index[XbInfo.size]]&1)*
                                    ////    (TrackInfo.isGoodCand[XbInfo.rftk2_index[XbInfo.size]]&1)==1)
                                    ////    XbInfo.isGoodCand[XbInfo.size] += 1 << 3;
                                    XbInfo.isGoodCand[XbInfo.size] += 1 << 2;
                                    if (xbVFP->currentState().mass() > 7)
                                        XbInfo.isGoodCand[XbInfo.size] += 1 << 1;
                                    if (((XbInfo.isGoodCand[XbInfo.size]>>2)&((XbInfo.isGoodCand[XbInfo.size]>>2)+1))==0)
                                        XbInfo.isGoodCand[XbInfo.size] += 1;
        
                                    Xb_candidate.clear();
                                    xCands.clear();
                                    XbInfo.size++;
                                }//Tk2
                            }//Tk1              
                        }//Mu2
                    }//Mu1}}}
                    //printf("-----*****DEBUG:End of XbInfo.\n");

                    // TrackInfo section {{{
                    const reco::GenParticle* genTrackPtr[MAX_GEN];
                    memset(genTrackPtr,0x00,MAX_GEN);
                    for(std::vector<pat::GenericParticle>::const_iterator tk_it=input_tracks.begin();
                        tk_it != input_tracks.end() ; tk_it++){
                        if(TrackInfo.size >= MAX_TRACK){
                            fprintf(stderr,"ERROR: number of tracks exceeds the size of array.\n");
                            break;;
                        }
                        int tk_hindex = int(tk_it - input_tracks.begin());
                        if (isNeededTrack[tk_hindex]==false) continue;

                        //Create list of relative xb candidates for later filling
                        std::vector<int> listOfRelativeXbCands;//1~nXb
                        for(int iXb=0; iXb < XbInfo.size; iXb++){
                            if(XbInfo.rftk1_index[iXb] == -tk_hindex-1){
                                listOfRelativeXbCands.push_back(iXb+1);
                            }else if(XbInfo.rftk2_index[iXb] == -tk_hindex-1){
                                listOfRelativeXbCands.push_back(-iXb-1);
                            }
                        }
                        if (listOfRelativeXbCands.size() == 0) continue;//drop unused tracks

                        TrackInfo.index          [TrackInfo.size] = TrackInfo.size;
                        TrackInfo.handle_index   [TrackInfo.size] = tk_hindex;
                        TrackInfo.charge         [TrackInfo.size] = tk_it->charge();
                        TrackInfo.pt             [TrackInfo.size] = tk_it->pt();
                        TrackInfo.eta            [TrackInfo.size] = tk_it->eta();
                        TrackInfo.phi            [TrackInfo.size] = tk_it->phi();
                        //TrackInfo.p              [TrackInfo.size] = tk_it->p();
                        TrackInfo.striphit       [TrackInfo.size] = tk_it->track()->hitPattern().numberOfValidStripHits();
                        TrackInfo.pixelhit       [TrackInfo.size] = tk_it->track()->hitPattern().numberOfValidPixelHits();
                        TrackInfo.fpbarrelhit    [TrackInfo.size] = tk_it->track()->hitPattern().hasValidHitInFirstPixelBarrel();
                        TrackInfo.fpendcaphit    [TrackInfo.size] = tk_it->track()->hitPattern().hasValidHitInFirstPixelEndcap();
                        TrackInfo.chi2           [TrackInfo.size] = tk_it->track()->chi2();
                        TrackInfo.ndf            [TrackInfo.size] = tk_it->track()->ndof();
                        TrackInfo.d0             [TrackInfo.size] = tk_it->track()->d0();
                        TrackInfo.d0error        [TrackInfo.size] = tk_it->track()->d0Error();
                        TrackInfo.dzPV           [TrackInfo.size] = tk_it->track()->dz(RefVtx);
                        TrackInfo.dxyPV          [TrackInfo.size] = tk_it->track()->dxy(RefVtx);
                        if (!iEvent.isRealData())
                            genTrackPtr [TrackInfo.size] = tk_it->genParticle();

                        //1. make sure non-zero for binary test
                        //0. isGoodCand
                        TrackInfo.isGoodCand     [TrackInfo.size]+= 1 << 1;
                        if (((TrackInfo.isGoodCand[TrackInfo.size]>>1)&((TrackInfo.isGoodCand[TrackInfo.size]>>1)+1))==0)
                            TrackInfo.isGoodCand     [TrackInfo.size]+= 1;

                        //Fill correct track index and track quality to correspond Xb candidate
                        for(unsigned int iCands=0; iCands < listOfRelativeXbCands.size(); iCands++){
                            if (listOfRelativeXbCands[iCands]>0){
                                XbInfo.rftk1_index[listOfRelativeXbCands[iCands]-1] = TrackInfo.size;
                                if ((TrackInfo.isGoodCand[TrackInfo.size]&1) == 1){
                                    if (XbInfo.rftk2_index[listOfRelativeXbCands[iCands]-1]>=0){
                                        XbInfo.isGoodCand[listOfRelativeXbCands[iCands]-1] += ((XbInfo.isGoodCand[listOfRelativeXbCands[iCands]-1]>>2)&1)<<2;
                                    }else{
                                        XbInfo.isGoodCand[listOfRelativeXbCands[iCands]-1] += 1 << 2;
                                    }
                                }
                            }else{
                                XbInfo.rftk2_index[-listOfRelativeXbCands[iCands]-1] = TrackInfo.size;
                                if ((TrackInfo.isGoodCand[TrackInfo.size]&1) == 1){
                                    if (XbInfo.rftk1_index[-listOfRelativeXbCands[iCands]-1]>=0){
                                        XbInfo.isGoodCand[-listOfRelativeXbCands[iCands]-1] += ((XbInfo.isGoodCand[-listOfRelativeXbCands[iCands]-1]>>2)&1)<<2;
                                    }else{
                                        XbInfo.isGoodCand[-listOfRelativeXbCands[iCands]-1] += 1 << 2;
                                    }
                                }
                            }
                        }
                        TrackInfo.size++;
                    }//end of TrackInfo}}}
                    //printf("-----*****DEBUG:End of TrackInfo.\n");

                    // GenInfo section{{{
                    if (!iEvent.isRealData()){
                        edm::Handle< std::vector<reco::GenParticle> > gens;
                        iEvent.getByLabel(genLabel_, gens);
    
                        std::vector<const reco::Candidate *> cands;
                        for(std::vector<reco::GenParticle>::const_iterator it_gen = gens->begin();
                            it_gen != gens->end(); it_gen++ ){
                            cands.push_back(&*it_gen);
                        }
    
                        for(std::vector<reco::GenParticle>::const_iterator it_gen=gens->begin();
                            it_gen != gens->end(); it_gen++){
                            if (it_gen->status() > 2)                           continue;
                            //if (it_gen->pdgId() == 92)                          continue;
                            //if (it_gen->pdgId() == 21)                          continue;
                            //if (it_gen->pdgId() == 22 && it_gen->status() != 1) continue;
                            if (abs(it_gen->pdgId()) != 211 &&
                                abs(it_gen->pdgId()) != 13  &&
                                it_gen->pdgId() != 443      &&
                                it_gen->pdgId() != 100443   &&
                                it_gen->pdgId() != 553      &&
                                it_gen->pdgId() != 100553
                               ) continue;

                            bool isGenSignal = false;
                            if (abs(it_gen->pdgId()) == 13                              &&
                                it_gen->numberOfMothers() == 1                          &&
                                (it_gen->mother()->pdgId() == 443 ||
                                 it_gen->mother()->pdgId() == 553 )                     &&
                                it_gen->mother()->numberOfDaughters() >= 2              &&
                                it_gen->mother()->numberOfMothers() == 1                &&
                                (it_gen->mother()->mother()->pdgId() == 100553 ||
                                 it_gen->mother()->mother()->pdgId() == 100443 )        &&
                                it_gen->mother()->mother()->numberOfDaughters() == 3    &&
                                abs(it_gen->mother()->mother()->daughter(1)->pdgId()) == 211 
                               ) isGenSignal = true;//signal mu

                            if (abs(it_gen->pdgId()) == 211                             &&
                                it_gen->numberOfMothers() == 1                          &&
                                (it_gen->mother()->mother()->pdgId() == 100553 ||
                                 it_gen->mother()->mother()->pdgId() == 100443 )        &&
                                it_gen->mother()->numberOfDaughters()== 3               &&
                                (it_gen->mother()->daughter(0)->pdgId() == 553 ||
                                 it_gen->mother()->daughter(0)->pdgId() == 443 )        &&
                                it_gen->mother()->daughter(0)->numberOfDaughters()>=2   &&
                                abs(it_gen->mother()->daughter(0)->daughter(0)->pdgId()) == 13
                               ) isGenSignal = true;//signal pi
                            
                            if ((it_gen->pdgId() == 443 || it_gen->pdgId() == 553)      &&
                                (it_gen->mother()->pdgId() == 100443 ||
                                 it_gen->mother()->pdgId() == 100553 )                  &&
                                it_gen->mother()->numberOfDaughters() == 3              &&
                                abs(it_gen->mother()->daughter(1)->pdgId()) == 211      &&
                                it_gen->numberOfDaughters() >= 2                        &&
                                abs(it_gen->daughter(0)->pdgId()) == 13    
                               ) isGenSignal = true;//signal uj

                            if ((it_gen->pdgId()==100443 || it_gen->pdgId()==100553)    &&
                                it_gen->numberOfDaughters() == 3                        &&
                                it_gen->daughter(0)->numberOfDaughters() >= 2           &&
                                abs(it_gen->daughter(0)->daughter(0)->pdgId()) == 13    &&
                                abs(it_gen->daughter(1)->pdgId()) == 211 
                               ) isGenSignal = true;//signal xb

                            if (!isGenSignal) continue;

                            int iMo1 = -1,  iMo2 = -1,  iDa1 = -1,  iDa2 = -1;
                            for(std::vector<const reco::Candidate *>::iterator iCands = cands.begin();
                                iCands != cands.end(); iCands++){
                                if (it_gen->numberOfMothers() >= 2){
                                    if (it_gen->mother(0) == *iCands)
                                        iMo1 = iCands - cands.begin();
                                    if (it_gen->mother(1) == *iCands)
                                        iMo2 = iCands - cands.begin();
                                }else if(it_gen->numberOfMothers() == 1){
                                    if (it_gen->mother(0) == *iCands)
                                        iMo1 = iCands - cands.begin();
                                }
                                if (it_gen->numberOfDaughters() >= 2){
                                    if (it_gen->daughter(0) == *iCands)
                                        iDa1 = iCands - cands.begin();
                                    else if (it_gen->daughter(1) == *iCands)
                                        iDa2 = iCands - cands.begin();
                                }else if(it_gen->numberOfDaughters() == 1){
                                    if (it_gen->daughter(0) == *iCands)
                                        iDa1 = iCands - cands.begin();
                                }
                            }
    
                            //depend on TrackInfo.handle_index and MuonInfo.handle_index
                            //printf("-----*****DEBUG:Start of matching.\n");
                            if (abs(it_gen->pdgId()) == 13){
                                //printf("-----*****DEBUG:Entered muon matching block.\n");
                                int best_mu_index = -1;
                                for(int muonIdx = 0; muonIdx < MuonInfo.size; muonIdx++){
                                    ////Private matching
                                    //if (fabs(MuonInfo.pt [muonIdx]/it_gen->pt ()-1) > 0.1) continue;
                                    //if (fabs(MuonInfo.eta[muonIdx]/it_gen->eta()-1) > 0.1) continue;
                                    //if (fabs(MuonInfo.phi[muonIdx]/it_gen->phi()-1) > 0.1) continue;
                                    //if (best_mu_index == -1){
                                    //    printf("-----*****DEBUG:[Mu]Tar.Pt /Ref.Pt = %9f/%9f\n",
                                    //        MuonInfo.pt [muonIdx],it_gen->pt ());
                                    //    printf("-----*****DEBUG:[Mu]Tar.Eta/Ref.Eta= %9f/%9f\n",
                                    //        MuonInfo.eta[muonIdx],it_gen->eta());
                                    //    printf("-----*****DEBUG:[Mu]Tar.Phi/Ref.Phi= %9f/%9f\n",
                                    //        MuonInfo.phi[muonIdx],it_gen->phi());
                                    //    best_mu_index = muonIdx;
                                    //}else if (sqrt(pow(MuonInfo.pt [muonIdx]/it_gen->pt()-1,2)+
                                    //          pow(MuonInfo.eta[muonIdx]/it_gen->eta()-1,2)+
                                    //          pow(MuonInfo.phi[muonIdx]/it_gen->phi()-1,2))-
                                    //          sqrt(pow(MuonInfo.pt[best_mu_index]/it_gen->pt()-1,2)-
                                    //          pow(MuonInfo.eta[best_mu_index]/it_gen->eta()-1,2)-
                                    //          pow(MuonInfo.phi[best_mu_index]/it_gen->phi()-1,2))<=0){
                                    //    printf("-----*****DEBUG:[Mu]Tar.Pt /Ref.Pt = %9f/%9f\n",
                                    //        MuonInfo.pt [muonIdx],it_gen->pt ());
                                    //    printf("-----*****DEBUG:[Mu]Tar.Eta/Ref.Eta= %9f/%9f\n",
                                    //        MuonInfo.eta[muonIdx],it_gen->eta());
                                    //    printf("-----*****DEBUG:[Mu]Tar.Phi/Ref.Phi= %9f/%9f\n",
                                    //        MuonInfo.phi[muonIdx],it_gen->phi());
                                    //    best_mu_index = muonIdx;
                                    //}
                                    
                                    // match by pat::Muon
                                    if (genMuonPtr[muonIdx] == 0) continue;
                                    if (it_gen->p4() == genMuonPtr[muonIdx]->p4()){
                                        best_mu_index = muonIdx;
                                        //printf("-----*****DEBUG:[Mu]Tar.Pt /Ref.Pt = %9f/%9f\n",
                                        //    MuonInfo.pt [muonIdx],it_gen->pt ());
                                        //printf("-----*****DEBUG:[Mu]Tar.Eta/Ref.Eta= %9f/%9f\n",
                                        //    MuonInfo.eta[muonIdx],it_gen->eta());
                                        //printf("-----*****DEBUG:[Mu]Tar.Phi/Ref.Phi= %9f/%9f\n",
                                        //    MuonInfo.phi[muonIdx],it_gen->phi());
                                        break;
                                    }
                                }
                                if (GenInfo.mhmu1_index == -1){
                                    GenInfo.mhmu1_index = best_mu_index;
                                    //printf("-----*****DEBUG:mhmu1_index=%d\n",GenInfo.mhmu1_index);
                                }else{
                                    GenInfo.mhmu2_index = best_mu_index;
                                    //printf("-----*****DEBUG:mhmu2_index=%d\n",GenInfo.mhmu2_index);
                                }
                            }else if(abs(it_gen->pdgId()) == 211){
                                //printf("-----*****DEBUG:Entered pion matching block.\n");
                                int best_tk_index = -1;
                                for(int trackIdx = 0; trackIdx < TrackInfo.size; trackIdx++){
                                    // Private matching
                                    //if (fabs(TrackInfo.pt [trackIdx]/it_gen->pt ()-1) > 0.1) continue;
                                    //if (fabs(TrackInfo.eta[trackIdx]/it_gen->eta()-1) > 0.1) continue;
                                    //if (fabs(TrackInfo.phi[trackIdx]/it_gen->phi()-1) > 0.1) continue;
                                    //if (best_tk_index == -1){
                                    //    printf("-----*****DEBUG:[Tk]Tar.Pt /Ref.Pt = %9f/%9f\n",
                                    //        TrackInfo.pt [trackIdx],it_gen->pt ());
                                    //    printf("-----*****DEBUG:[Tk]Tar.Eta/Ref.Eta= %9f/%9f\n",
                                    //        TrackInfo.eta[trackIdx],it_gen->eta());
                                    //    printf("-----*****DEBUG:[Tk]Tar.Phi/Ref.Phi= %9f/%9f\n",
                                    //        TrackInfo.phi[trackIdx],it_gen->phi());
                                    //    best_tk_index = trackIdx;
                                    //}else if(sqrt(pow(TrackInfo.pt [trackIdx]/it_gen->pt ()-1,2)+
                                    //         pow(TrackInfo.eta[trackIdx]/it_gen->eta()-1,2)+
                                    //         pow(TrackInfo.phi[trackIdx]/it_gen->phi()-1,2))-
                                    //         sqrt(pow(TrackInfo.pt[best_tk_index]/it_gen->pt()-1,2)-
                                    //         pow(TrackInfo.eta[best_tk_index]/it_gen->eta()-1,2)-
                                    //         pow(TrackInfo.phi[best_tk_index]/it_gen->phi()-1,2))<=0){
                                    //    printf("-----*****DEBUG:[Tk]Tar.Pt /Ref.Pt = %9f/%9f\n",
                                    //        TrackInfo.pt [trackIdx],it_gen->pt ());
                                    //    printf("-----*****DEBUG:[Tk]Tar.Eta/Ref.Eta= %9f/%9f\n",
                                    //        TrackInfo.eta[trackIdx],it_gen->eta());
                                    //    printf("-----*****DEBUG:[Tk]Tar.Phi/Ref.Phi= %9f/%9f\n",
                                    //        TrackInfo.phi[trackIdx],it_gen->phi());
                                    //    best_tk_index = trackIdx;
                                    //}
                                    
                                    //match by pat::GenericParticle
                                    if (genTrackPtr[trackIdx] == 0 ) continue;
                                    if (it_gen->p4() == genTrackPtr[trackIdx]->p4()){
                                        best_tk_index = trackIdx;
                                        //printf("-----*****DEBUG:[Tk]Tar.Pt /Ref.Pt = %9f/%9f\n",
                                        //    TrackInfo.pt [trackIdx],it_gen->pt ());
                                        //printf("-----*****DEBUG:[Tk]Tar.Eta/Ref.Eta= %9f/%9f\n",
                                        //    TrackInfo.eta[trackIdx],it_gen->eta());
                                        //printf("-----*****DEBUG:[Tk]Tar.Phi/Ref.Phi= %9f/%9f\n",
                                        //    TrackInfo.phi[trackIdx],it_gen->phi());
                                        break;
                                    }
                                }
                                if (GenInfo.mhtk1_index == -1){
                                    GenInfo.mhtk1_index = best_tk_index;
                                    //printf("-----*****DEBUG:mhtk1_index=%d\n",GenInfo.mhtk1_index);
                                }else{
                                    GenInfo.mhtk2_index = best_tk_index;
                                    //printf("-----*****DEBUG:mhtk2_index=%d\n",GenInfo.mhtk2_index);
                                }
                            }
                            //printf("-----*****DEBUG:mhmu1_index=%d\n",GenInfo.mhmu1_index);
                            //printf("-----*****DEBUG:mhmu2_index=%d\n",GenInfo.mhmu2_index);
                            //printf("-----*****DEBUG:mhtk1_index=%d\n",GenInfo.mhtk1_index);
                            //printf("-----*****DEBUG:mhtk2_index=%d\n",GenInfo.mhtk2_index);
    
                            GenInfo.index[GenInfo.size]         = GenInfo.size;
                            GenInfo.handle_index[GenInfo.size]  = it_gen-gens->begin();
                            GenInfo.pt[GenInfo.size]            = it_gen->pt();
                            GenInfo.eta[GenInfo.size]           = it_gen->eta();
                            GenInfo.phi[GenInfo.size]           = it_gen->phi();
                            GenInfo.mass[GenInfo.size]          = it_gen->mass();
                            GenInfo.pdgId[GenInfo.size]         = it_gen->pdgId();
                            GenInfo.status[GenInfo.size]        = it_gen->status();
                            GenInfo.nMo[GenInfo.size]           = it_gen->numberOfMothers();
                            GenInfo.nDa[GenInfo.size]           = it_gen->numberOfDaughters();
                            GenInfo.mo1[GenInfo.size]           = iMo1;//To be matched later.
                            GenInfo.mo2[GenInfo.size]           = iMo2;
                            GenInfo.da1[GenInfo.size]           = iDa1;
                            GenInfo.da2[GenInfo.size]           = iDa2;
                            GenInfo.size++;
                        }
                        //printf("-----*****DEBUG:End of gens loop.\n");
                        
                        if (GenInfo.mhmu1_index != -1     &&
                            GenInfo.mhmu2_index != -1     ){
                            TLorentzVector mhmu1_4vec,mhmu2_4vec;
                            mhmu1_4vec.SetPtEtaPhiM(
                                MuonInfo.pt [GenInfo.mhmu1_index],
                                MuonInfo.eta[GenInfo.mhmu1_index],
                                MuonInfo.phi[GenInfo.mhmu1_index],
                                MUON_MASS);
                            mhmu2_4vec.SetPtEtaPhiM(
                                MuonInfo.pt [GenInfo.mhmu2_index],
                                MuonInfo.eta[GenInfo.mhmu2_index],
                                MuonInfo.phi[GenInfo.mhmu2_index],
                                MUON_MASS);
                            GenInfo.mhujMass = (mhmu1_4vec+mhmu2_4vec).Mag();
                            //printf("-----*****DEBUG:mhujMass=%9f\n",(mhmu1_4vec+mhmu2_4vec).Mag());
                            if (GenInfo.mhtk1_index != -1     &&
                                GenInfo.mhtk2_index != -1     ){
                                TLorentzVector mhtk1_4vec,mhtk2_4vec;
                                mhtk1_4vec.SetPtEtaPhiM(
                                    TrackInfo.pt [GenInfo.mhtk1_index],
                                    TrackInfo.eta[GenInfo.mhtk1_index],
                                    TrackInfo.phi[GenInfo.mhtk1_index],
                                    PION_MASS);
                                mhtk2_4vec.SetPtEtaPhiM(
                                    TrackInfo.pt [GenInfo.mhtk2_index],
                                    TrackInfo.eta[GenInfo.mhtk2_index],
                                    TrackInfo.phi[GenInfo.mhtk2_index],
                                    PION_MASS);
                                GenInfo.mhxbMass = (mhmu1_4vec+mhmu2_4vec+mhtk1_4vec+mhtk2_4vec).Mag();
                                //printf("-----*****DEBUG:mhxbMass=%9f\n",(mhmu1_4vec+mhmu2_4vec+mhtk1_4vec+mhtk2_4vec).Mag());
                            }
                        }
    
                        //Pass handle_index to igen
                        for(int igen = 0; igen < GenInfo.size; igen++){
                            int iMo1 = GenInfo.mo1[igen];
                            int iMo2 = GenInfo.mo2[igen];
                            int iDa1 = GenInfo.da1[igen];
                            int iDa2 = GenInfo.da2[igen];
                            for(int k = 0; k < GenInfo.size; k++){
                                if (iMo1 == GenInfo.handle_index[k])
                                    GenInfo.mo1[igen] = k;
                                else if (iMo2 == GenInfo.handle_index[k])
                                    GenInfo.mo2[igen] = k;
                                else if (iDa1 == GenInfo.handle_index[k])
                                    GenInfo.da1[igen] = k;
                                else if (iDa2 == GenInfo.handle_index[k])
                                    GenInfo.da2[igen] = k;
                            }
                            //In case that GEN particles are omitted from GenInfo
                            //handle_index couldn't be the same as igen
                            //since the very first proton pair has status 3.
                            if (iMo1 == GenInfo.mo1[igen])
                                GenInfo.mo1[igen] = -1;
                            if (iMo2 == GenInfo.mo2[igen])
                                GenInfo.mo2[igen] = -1;
                            if (iDa1 == GenInfo.da1[igen])
                                GenInfo.da1[igen] = -1;
                            if (iDa2 == GenInfo.da2[igen])
                               GenInfo.da2[igen] = -1;
                        }
                        //printf("-----*****DEBUG:End of IndexToIgen\n");
                    }//isRealData}}}
                    //printf("-----*****DEBUG:End of GenInfo.\n");
                    //std::cout<<"Start to fill!\n";
                }//has nMuons>1 nTracks>1
            }//Tracks
        }//Muonss
    }//try
    catch (std::exception & err){
            std::cout  << "Exception during event number: " << iEvent.id()
                << "\n" << err.what() << "\n";
    }//catch 
    std::cout << "XbInfo.size=" << XbInfo.size << std::endl;
    root->Fill();
    //std::cout<<"filled!\n";
}


// ------------ method called once each job just after ending the event loop  ------------
void Xb_frame::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void Xb_frame::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void Xb_frame::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void Xb_frame::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void Xb_frame::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Xb_frame::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Xb_frame);
