#include <correlations/Types.hh>
#include <correlations/Result.hh>
#include <correlations/QVector.hh>
#include <correlations/recursive/FromQVector.hh>
#include <correlations/recurrence/FromQVector.hh>
#include <correlations/closed/FromQVector.hh>
#include <TComplex.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TNtupleD.h>
#include <TRandom3.h>
#include <TFile.h>
//#include "QWConstV3.h"
#include <RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h>
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//
// constants, enums and typedefs
//

//#define QW_DEBUG 1
//#define QW_PEREVENT 1

#define PRD(x) cout << "!!QW!! " << __LINE__ << " DEBUG OUTPUT " << (#x) << " = " << (x) << endl;
#define PR(x) cout << "!!QW!! " << __LINE__ << " DEBUG OUTPUT " << (#x) << endl;
//
// class declaration
//

///////////////// Class ////////////////////////////

class QWCumuPt : public edm::EDAnalyzer {
	public:
		explicit QWCumuPt(const edm::ParameterSet&);
		~QWCumuPt();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		virtual void beginRun(edm::Run const&, edm::EventSetup const&);
		virtual void endRun(edm::Run const&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

	/////////////////////////////////////////////
		//TRandom3 * gRandom;
		// ----------member data ---------------------------

		edm::InputTag					trackEta_;
		edm::InputTag					trackPhi_;
		edm::InputTag					trackPt_;
		edm::InputTag					trackWeight_;
		edm::InputTag					vertexZ_;

		edm::InputTag					centralityTag_;

		double	minvz_, maxvz_;
		std::vector<double>		ptBin_;
		int				Npt_;
	/////////////////////////////////////////////
		double	rfpmineta_, rfpmaxeta_;
		double	rfpminpt_, rfpmaxpt_;
		double	dEtaGap_;

		int	cmode_;
	/////////////////////////////////////////////
		TTree * trV;

		int gNoff;
		int gMult;

		// RFP
		double rQ[7][4];
		double iQ[7][4];
		double wQ[7][4];

		// RFP Gap
		double rQGap[7];
		double iQGap[7];
		double wQGap[7];

		// 2-part Gap
		double rQpGap[7][24];
		double wQpGap[7][24];

		// m-part pT
		double rVQp[7][4][24];
		double iVQp[7][4][24];
		double wVQp[7][4][24];

		correlations::HarmonicVector	hc[7];
		correlations::QVector		q[7];
		correlations::FromQVector	*cq[7];

		void initQ();
		void doneQ();
};



