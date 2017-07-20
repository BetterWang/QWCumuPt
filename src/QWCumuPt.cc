// -*- C++ -*-
//
// Package:    QWCumuPt
// Class:      QWCumuPt
// 
/**\class QWCumuPt QWCumuPt.cc QWAna/QWCumuPt/src/QWCumuPt.cc

Description: [Cumulant differential package]

Implementation:
[Quan Wang is cool]
*/
//
// Original Author:  Quan Wang
//         Created:  05/23/2014
// $Id: QWCumuPt.cc,v 1.0 2017/02/11 15:56:58 qwang Exp $
//
//


// system include files
#include <memory>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "TComplex.h"
#include <complex>


#include "QWAna/QWCumuPt/interface/QWCumuPt.h"


using namespace std;

//#ifdef QW_DEBUG
//
// constructors and destructor
//
QWCumuPt::QWCumuPt(const edm::ParameterSet& iConfig):
	vertexZ_( iConfig.getUntrackedParameter<edm::InputTag>("vertexZ") ),
	centralityTag_( iConfig.getUntrackedParameter<edm::InputTag>("centrality") )
{
	const edm::ParameterSet& track = iConfig.getUntrackedParameter<edm::ParameterSet>("trackSet");
	trackEta_ = track.getUntrackedParameter<edm::InputTag>("Eta");
	trackPhi_ = track.getUntrackedParameter<edm::InputTag>("Phi");
	trackPt_  = track.getUntrackedParameter<edm::InputTag>("Pt");
	trackWeight_  = track.getUntrackedParameter<edm::InputTag>("Weight");

	//now do what ever initialization is needed
	minvz_ = iConfig.getUntrackedParameter<double>("minvz", -15.);
	maxvz_ = iConfig.getUntrackedParameter<double>("maxvz", 15.);

	ptBin_ = iConfig.getUntrackedParameter< std::vector<double> >("ptBin");
	Npt_ = ptBin_.size();

	rfpmineta_ = iConfig.getUntrackedParameter<double>("rfpmineta", -2.4);
	rfpmaxeta_ = iConfig.getUntrackedParameter<double>("rfpmaxeta", 2.4);
	rfpminpt_ = iConfig.getUntrackedParameter<double>("rfpminpt", 0.3);
	rfpmaxpt_ = iConfig.getUntrackedParameter<double>("rfpmaxpt", 3.0);

	dEtaGap_ = iConfig.getUntrackedParameter<double>("EtaGap", 2.0);

	cmode_ = iConfig.getUntrackedParameter<int>("cmode", 1);

        consumes<int>(centralityTag_);
        consumes<std::vector<double> >(vertexZ_);

        consumes<std::vector<double> >(trackEta_);
        consumes<std::vector<double> >(trackPhi_);
        consumes<std::vector<double> >(trackPt_);
        consumes<std::vector<double> >(trackWeight_);

	for ( int n = 1; n < 7; n++ ) {
		q[n] = correlations::QVector(0, 0, true);
	}

	edm::Service<TFileService> fs;

	trV = fs->make<TTree>("trV", "trV");
	trV->Branch("Noff", &gNoff, "Noff/I");
	trV->Branch("Mult", &gMult, "Mult/I");


	for ( int np = 0; np < 4; np++ ) {
		for ( int n = 2; n < 7; n++ ) {
			trV->Branch(Form("rQ%i%i", n, 2+2*np), &rQ[n][np], Form("rQ%i%i/D", n, 2+2*np));

			trV->Branch(Form("rVQp%i%i", n, 2+2*np), &rVQp[n][np], Form("rVQp%i%i[24]/D", n, 2+2*np));
		}

		int n = 2;
		trV->Branch(Form("wQ%i%i", n, 2+2*np), &wQ[n][np], Form("wQ%i%i/D", n, 2+2*np));
		trV->Branch(Form("wVQp%i%i", n, 2+2*np), wVQp[n][np], Form("wVQp%i%i[24]/D", n, 2+2*np));
	}

	int n = 2;
	trV->Branch(Form("wQGap%i", n), &wQGap[n], Form("wQGap%i/D", n));
	trV->Branch(Form("wQpGap%i", n), &wQpGap[n], Form("wQpGap%i[24]/D", n));

	for ( int n = 2; n < 7; n++ ) {
		trV->Branch(Form("rQGap%i", n), &rQGap[n], Form("rQGap%i/D", n));

		trV->Branch(Form("rQpGap%i", n), rQpGap[n], Form("rQpGap%i[24]/D", n));
	}

	cout << " cmode_ = " << cmode_ << endl;

	initQ();
}


QWCumuPt::~QWCumuPt()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

void
QWCumuPt::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	Handle<std::vector<double> >	hVz;

	Handle<std::vector<double> >	hEta;
	Handle<std::vector<double> >	hPhi;
	Handle<std::vector<double> >	hPt;
	Handle<std::vector<double> >	hWeight;


	iEvent.getByLabel(trackEta_,	hEta);
	iEvent.getByLabel(trackPhi_,	hPhi);
	iEvent.getByLabel(trackPt_,	hPt);
	iEvent.getByLabel(trackWeight_, hWeight);

	iEvent.getByLabel(vertexZ_, 	hVz);

	if ( hVz->size() < 1 ) return;
	if ( fabs((*hVz)[0]) > maxvz_ or fabs((*hVz)[0]) < minvz_ ) return;
	int sz = int(hEta->size());
	if ( sz == 0 ) return;

	std::vector<int>	RFP;
	RFP.reserve(sz);
	int rfp_sz = 0;
	for ( int i = 0; i < sz; i++ ) {
		if ( (*hEta)[i] < rfpmaxeta_ and (*hEta)[i] > rfpmineta_
		and (*hPt)[i] < rfpmaxpt_ and (*hPt)[i] > rfpminpt_ ) {
			RFP[i] = 1;
			rfp_sz++;
		} else {
			RFP[i] = 0;
		}
	}

	for ( int n = 0; n < 7; n++ ) {
		for ( int np = 0; np < 4; np++ ) {
			rQ[n][np] = 0;
			iQ[n][np] = 0;
			wQ[n][np] = 0;

			for ( int i = 0; i < 24; i++ ) {
				rVQp[n][np][i] = 0;
				iVQp[n][np][i] = 0;
				wVQp[n][np][i] = 0;
			}
		}

		rQGap[n] = 0;
		wQGap[n] = 0;

		for ( int i = 0; i < 24; i++ ) {
			rQpGap[n][i] = 0;
			wQpGap[n][i] = 0;
		}
	}

	for ( int i = 0; i < sz; i++ ) {
		if ( RFP[i] != 1 ) continue;
		for ( int n = 1; n < 7; n++ ) {
			q[n].fill((*hPhi)[i], (*hWeight)[i]);
		}
	}

	correlations::Result r[7][4];
	for ( int n = 1; n < 7; n++ ) {
		for ( int np = 0; np < 4; np++ ) {
			r[n][np] = cq[n]->calculate(2+2*np, hc[n]);
		}
	}

	// RFP
	for ( int n = 2; n < 7; n++ ) {
		for ( int np = 0; np < 4; np++ ) {
			rQ[n][np] = r[n][np].sum().real();
			iQ[n][np] = r[n][np].sum().imag();
			wQ[n][np] = r[n][np].weight();
		}
	}

	// 2part RFP
	for ( int i = 0; i < sz; i++ ) {
		if ( RFP[i] != 1 ) continue;
		for ( int j = i+1; j < sz; j++ ) {
			if ( RFP[j] != 1 ) continue;
			if ( fabs( (*hEta)[i] - (*hEta)[j] ) < dEtaGap_ ) continue;
			for ( int n = 2; n < 7; n++ ) {
				rQGap[n] += cos( n*( (*hPhi)[i] - (*hPhi)[j] ) ) * (*hWeight)[i] * (*hWeight)[j];
				wQGap[n] += (*hWeight)[i] * (*hWeight)[j];
			}
		}
	}

	// 2part pt
	for ( int i = 0; i < sz; i++ ) {
		int ipt = 0;
		while ( (*hPt)[i] > ptBin_[ipt+1] ) ipt++;
		for ( int j = i+1; j < sz; j++ ) {
			if ( RFP[j] != 1 ) continue;
			if ( fabs( (*hEta)[i] - (*hEta)[j] ) < dEtaGap_ ) continue;
			for ( int n = 2; n < 7; n++ ) {
				rQpGap[n][ipt] += cos( n*( (*hPhi)[i] - (*hPhi)[j] ) ) * (*hWeight)[i] * (*hWeight)[j];
				wQpGap[n][ipt] += (*hWeight)[i] * (*hWeight)[j];
			}
		}
	}

	for ( int n = 1; n < 7; n++ ) {
		for ( int np = 0; np < 4; np++ ) {
			correlations::Complex qp = 0;
			double wt = 0;
			for ( int ipt = 0; ipt < Npt_; ipt++ ) {
				qp = 0;
				wt = 0;
				for ( int i = 0; i < sz; i++ ) {
					if ( (*hPt)[i] < ptBin_[ipt] or (*hPt)[i] > ptBin_[ipt+1] ) continue;
					correlations::QVector tq = q[n];
					for ( int j = i+1; j < sz; j++ ) {
						if ( RFP[j] == 1 ) {
							tq.unfill( (*hPhi)[j], (*hWeight)[j] );
						}
					}
					correlations::FromQVector *cq = 0;
					switch ( cmode_ ) {
						case 1:
							cq = new correlations::closed::FromQVector(tq);
							break;
						case 2:
							cq = new correlations::recurrence::FromQVector(tq);
							break;
						case 3:
							cq = new correlations::recursive::FromQVector(tq);
							break;
					}
					correlations::Result r = cq->calculate(np*2+1, hc[n]);
					qp += (*hWeight)[i] * correlations::Complex( TMath::Cos((*hPhi)[i] * n) , TMath::Sin((*hPhi)[i] * n) ) * r.sum();
					wt += (*hWeight)[i] * r.weight();
					delete cq;
				}
				rVQp[n][np][ipt] = qp.real();
				iVQp[n][np][ipt] = qp.imag();
				wVQp[n][np][ipt] = wt;
			}
		}
	}

	edm::Handle<int> ch;
	iEvent.getByLabel(centralityTag_,ch);
	gNoff = *ch;
	gMult = rfp_sz;

//	t->RunId = iEvent.id().run();
//	t->EventId = iEvent.id().event();
	trV->Fill();
	doneQ();
}


void
QWCumuPt::initQ()
{
	hc[1] = correlations::HarmonicVector(8);
	hc[1][0] = -1;
	hc[1][1] =  1;
	hc[1][2] = -1;
	hc[1][3] =  1;
	hc[1][4] = -1;
	hc[1][5] =  1;
	hc[1][6] = -1;
	hc[1][7] =  1;

	hc[2] = correlations::HarmonicVector(8);
	hc[2][0] = -2;
	hc[2][1] =  2;
	hc[2][2] = -2;
	hc[2][3] =  2;
	hc[2][4] = -2;
	hc[2][5] =  2;
	hc[2][6] = -2;
	hc[2][7] =  2;

	hc[3] = correlations::HarmonicVector(8);
	hc[3][0] = -3;
	hc[3][1] =  3;
	hc[3][2] = -3;
	hc[3][3] =  3;
	hc[3][4] = -3;
	hc[3][5] =  3;
	hc[3][6] = -3;
	hc[3][7] =  3;

	hc[4] = correlations::HarmonicVector(8);
	hc[4][0] = -4;
	hc[4][1] =  4;
	hc[4][2] = -4;
	hc[4][3] =  4;
	hc[4][4] = -4;
	hc[4][5] =  4;
	hc[4][6] = -4;
	hc[4][7] =  4;

	hc[5] = correlations::HarmonicVector(8);
	hc[5][0] = -5;
	hc[5][1] =  5;
	hc[5][2] = -5;
	hc[5][3] =  5;
	hc[5][4] = -5;
	hc[5][5] =  5;
	hc[5][6] = -5;
	hc[5][7] =  5;

	hc[6] = correlations::HarmonicVector(8);
	hc[6][0] = -6;
	hc[6][1] =  6;
	hc[6][2] = -6;
	hc[6][3] =  6;
	hc[6][4] = -6;
	hc[6][5] =  6;
	hc[6][6] = -6;
	hc[6][7] =  6;



	q[1].resize(hc[1]);
	q[2].resize(hc[2]);
	q[3].resize(hc[3]);
	q[4].resize(hc[4]);
	q[5].resize(hc[5]);
	q[6].resize(hc[6]);
	switch ( cmode_ ) {
		case 1:
			cq[1] = new correlations::closed::FromQVector(q[1]);
			cq[2] = new correlations::closed::FromQVector(q[2]);
			cq[3] = new correlations::closed::FromQVector(q[3]);
			cq[4] = new correlations::closed::FromQVector(q[4]);
			cq[5] = new correlations::closed::FromQVector(q[5]);
			cq[6] = new correlations::closed::FromQVector(q[6]);
			break;
		case 2:
			cq[1] = new correlations::recurrence::FromQVector(q[1]);
			cq[2] = new correlations::recurrence::FromQVector(q[2]);
			cq[3] = new correlations::recurrence::FromQVector(q[3]);
			cq[4] = new correlations::recurrence::FromQVector(q[4]);
			cq[5] = new correlations::recurrence::FromQVector(q[5]);
			cq[6] = new correlations::recurrence::FromQVector(q[6]);
			break;
		case 3:
			cq[1] = new correlations::recursive::FromQVector(q[1]);
			cq[2] = new correlations::recursive::FromQVector(q[2]);
			cq[3] = new correlations::recursive::FromQVector(q[3]);
			cq[4] = new correlations::recursive::FromQVector(q[4]);
			cq[5] = new correlations::recursive::FromQVector(q[5]);
			cq[6] = new correlations::recursive::FromQVector(q[6]);
			break;
	}
}

void
QWCumuPt::doneQ()
{
	q[1].reset();
	q[2].reset();
	q[3].reset();
	q[4].reset();
	q[5].reset();
	q[6].reset();
}

// ------------ method called once each job just before starting event loop  ------------
	void 
QWCumuPt::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
QWCumuPt::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
	void 
QWCumuPt::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
	void 
QWCumuPt::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
	void 
QWCumuPt::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
	void 
QWCumuPt::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
QWCumuPt::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);

	//Specify that only 'tracks' is allowed
	//To use, remove the default given above and uncomment below
	//ParameterSetDescription desc;
	//desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
	//descriptions.addDefault(desc);
}

//////////////////////////////////////////


//define this as a plug-in
DEFINE_FWK_MODULE(QWCumuPt);
