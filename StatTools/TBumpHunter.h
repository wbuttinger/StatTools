#ifndef TBUMPHUNTER_H
#define TBUMPHUNTER_H

///
///Implementation of BumpHunter algorithm (http://arxiv.org/pdf/1101.0390v2.pdf by G. Choudalakis)
///Author (of implementation): Will Buttinger (will@cern.ch)
///

#include <TH1D.h>
#include "TMath.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include "TRandom3.h"
#include "Math/GSLRndmEngines.h"
#include "Math/Random.h"
#include <vector>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TList.h>

class TBumpHunter : public TObject {


   public:

      enum TestStatisticType { BUMPHUNTER, DIPHUNTER, DISCREPENCYHUNTER, LOGLIKELIHOOD };

      struct Channel { 
         TH1* _histBack;
         TH1* _histData;
         TH1* _histPseudoData;
         Double_t _searchLowEdge;
         Double_t _searchHighEdge;
         Double_t _bumpWindowLowEdge;
         Double_t _bumpWindowHighEdge;
         Double_t _bumpPvalue;
         std::vector<std::pair<Int_t,Int_t> > _centralWindows;
         Channel() : _histBack(0), _histData(0), _histPseudoData(0), _searchLowEdge(0./0.), _searchHighEdge(0./0.), _bumpWindowLowEdge(0), _bumpWindowHighEdge(0), _bumpPvalue(1) {}
      };


      TBumpHunter(TH1* mc,TH1* data) : TObject() { init();SetChannelDistribution(mc,data); };
      TBumpHunter() : TObject() { init(); };
      ~TBumpHunter();

      ///specify the background distribution (e.g. the poisson means and the uncertainty on these means)
      void SetChannelDistribution(TH1* hist,TH1* data) { //fHistBack = hist; fHistData=data;
         if(fChannels.size()==0) fChannels.resize(1);
         fChannels[0]._histBack = hist; fChannels[0]._histData = data;
         if( fChannels[0]._histPseudoData ) delete fChannels[0]._histPseudoData;
         fChannels[0]._histPseudoData = (TH1*)data->Clone(Form("pseudo_%s",data->GetName())); 
         fChannels[0]._histPseudoData->SetDirectory(0);
         fChannels[0]._histPseudoData->Reset();
         SetChannel(0);
         EvaluateSearchPattern(); 
      }

      void AddChannelDistribution(TH1* hist,TH1* data) { 
         Channel t;
         t._histBack = hist; t._histData = hist; 
         t._histPseudoData = (TH1*)data->Clone(Form("pseudo_%s",data->GetName())); 
         t._histPseudoData->SetDirectory(0);
         t._histPseudoData->Reset();
         fChannels.push_back(t);
         SetChannel(fChannels.size()-1);
         EvaluateSearchPattern();
      }

      void SetChannel(Int_t in) { 
         fCurrentChannel = in;
      }

      void SetSearchRegion(Double_t low, Double_t high) { fChannels[fCurrentChannel]._searchLowEdge=low; fChannels[fCurrentChannel]._searchHighEdge=high; EvaluateSearchPattern(); }
      Double_t GetSearchLowEdge();
      Double_t GetSearchHighEdge();

      Double_t GetMinWindowSize();
      Double_t GetMaxWindowSize();

      Double_t GetWindowSizeStep();

      void SetMinWindowSize(Double_t in) { fMinWindowSize=in; kUseBinsForWindowSize=false; EvaluateSearchPattern();}
      void SetMaxWindowSize(Double_t in) { fMaxWindowSize=in; kUseBinsForMaxWindowSize=false; EvaluateSearchPattern();}

      void SetMinWindowBinSize(Int_t in) { fMinWindowSize=in; kUseBinsForWindowSize=true; EvaluateSearchPattern();}
      void SetMaxWindowBinSize(Int_t in) { fMaxWindowSize=in; kUseBinsForMaxWindowSize=true; EvaluateSearchPattern();}

      void SetWindowSizeStep(Double_t in) { fWindowSizeStep=in; kUseBinsForWindowSize=false; EvaluateSearchPattern();}
      void SetWindowSizeBinStep(Int_t in) { fWindowSizeStep=in; kUseBinsForWindowSize=true; EvaluateSearchPattern();}

      void SetNPseudoExperiments(Int_t n) { fnPseudo=n; }


      void SetBinModel(Int_t model) { if(model!=0&&model!=1&&model!=2&&model!=3) return; fBinModel=model;}
      void SetTestStatisticType(Int_t t) { fTestStatisticType = t; }

      Int_t GetNChannels() { return fChannels.size(); }

      Double_t GetBumpLowEdge() { return fChannels[fCurrentChannel]._bumpWindowLowEdge; }
      Double_t GetBumpHighEdge() { return fChannels[fCurrentChannel]._bumpWindowHighEdge; }

      Double_t GetGlobalPValue() { return fGlobalPValue; }
      Double_t GetLocalPValue() { return fLocalPValue; }

      ///do bumphunting. 0 = success
      Int_t Run();

      ///update the current pseudo-experiment histogram with a new set of values
      TH1* GenerateToyMC();

      TH1* GeneratePseudoData(TH1* mc); 

      TH1D* GetTestStatisticPDF() { return fBumpHunterStatisticPDF; }
      TGraph* GetPValueGraph() { return fBumpHunterStatisticConvergenceGraph; }

      ///evaluates the bumphunter test statistic for the given data
      Double_t EvaluateTestStatistic(TH1* data, Bool_t printOut=false);

      ///evaluate for multichannel
      Double_t EvaluateMultiChannelTestStatistic(Bool_t generatePseudo, Bool_t printOut=false);

      void EvaluateSearchPattern();
      void PrintSearchPattern();

      ///analytical determination of pvalue of nObs under poisson pdf convoluted with gamma distribution for the mean parameter (mean=E,variance=err^2)
      static Double_t GetPoissonConvGammaPValue(Double_t nObs, Double_t E, Double_t err);

      static Double_t GetPoissonPValue(Double_t nObs, Double_t E);
      static Double_t GetGaussianPValue(Double_t nObs, Double_t E, Double_t err);

      static Double_t GetZValue(Double_t pValue, bool& overflow);

      void SetBumpOverlapFactor(Double_t in) { fMultiChannelBumpOverlapFactor = in; }

      ClassDef(TBumpHunter,1);
   private:
      void init();
      Double_t fLocalPValue; //smallest p-value found in the actual data.. i.e. the bump's local pValue
      Double_t fGlobalPValue;

      std::vector<Channel> fChannels;

      TRandom3 pRand;
      ROOT::Math::Random<ROOT::Math::GSLRngMT> r; //used for gamma distribution random sampling

      Int_t fBinModel;
      Int_t fTestStatisticType;


      Double_t fMinWindowSize; Double_t fMaxWindowSize; Bool_t kUseBinsForWindowSize; Bool_t kUseBinsForMaxWindowSize;
      Double_t fWindowSizeStep;


      TH1D* fBumpHunterStatisticPDF;
      TGraphAsymmErrors* fBumpHunterStatisticConvergenceGraph;

      Int_t fnPseudo;
      UInt_t fCurrentChannel;
      Double_t fMultiChannelBumpOverlapFactor; //how much of the bump's regions must overlap. Default is 1 (perfect overlap)

      TCanvas* fBumpHunterResultCanvas;

};

#endif