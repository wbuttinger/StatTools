#include "StatTools/TBumpHunter.h"

#include <iostream>
#include <sstream>
#include <TLatex.h>
#include <TLine.h>
#include <cfloat>

#include "Math/QuantFuncMathCore.h"


ClassImp(TBumpHunter)

TBumpHunter::~TBumpHunter() {
}

void TBumpHunter::init() {
   fLocalPValue = 1;
   fGlobalPValue = 1;

//    fHistBack = 0;
//    fHistData = 0;
//    fHistPseudoData = 0;
// 
//    fHistBacks = new TList; fHistDatas = new TList; 
//    fHistPseudoDatas = new TList; fHistPseudoDatas->SetOwner(true); //owns its objects 

   fBinModel = 0; //poisson convoluted with gamma
   fTestStatisticType = BUMPHUNTER; 

//    fSearchLowEdge = 0./0.;
//    fSearchHighEdge = 0./0.;
// 
//    fBumpWindowLowEdge = 0.;fBumpWindowHighEdge=0.;

   fMinWindowSize = 0./0.;
   fMaxWindowSize = 0./0.;
   fWindowSizeStep = 0./0.;

   kUseBinsForWindowSize=false;

   pRand.SetSeed(123345);
   r.SetSeed(431155);

   fBumpHunterStatisticPDF = new TH1D("bumpHunterPDF","BumpHunter Test Statistic PDF",100,0,10);
   fBumpHunterStatisticPDF->SetDirectory(0);
   fBumpHunterStatisticPDF->GetXaxis()->SetTitle("BumpHunter Statistic");
   fBumpHunterStatisticPDF->GetYaxis()->SetTitle("Probability Density");

   fBumpHunterStatisticConvergenceGraph = 0;
   fBumpHunterResultCanvas = 0;
   fnPseudo = 0;
   fCurrentChannel = 0;
   fMultiChannelBumpOverlapFactor=1.;

}

Int_t TBumpHunter::Run() {
   delete fBumpHunterStatisticPDF;

   Channel& currChannel = fChannels[fCurrentChannel];


   currChannel._bumpWindowLowEdge = 0.;currChannel._bumpWindowHighEdge=0.;


   //evaluate for data 
   Double_t tObs = EvaluateMultiChannelTestStatistic(false);

   if(currChannel._bumpWindowLowEdge==0 && currChannel._bumpWindowHighEdge==0) {
      Error("Run","No bumps found?"); return -1;
   }

   fLocalPValue = exp(-tObs);
   //loop over the channels to show most significant bumps 
   for(unsigned int i=0;i<fChannels.size();i++) {
      bool overflow(false);
      double zVal = GetZValue(fChannels[i]._bumpPvalue,overflow);
      if(overflow) {
         Info("Run","Channel #%d: Observed most significant bump in [%f,%f] (> %f local sigma)",i,fChannels[i]._bumpWindowLowEdge,fChannels[i]._bumpWindowHighEdge,zVal);
      } else {
         Info("Run","Channel #%d: Observed most significant bump in [%f,%f] (%f local sigma)",i,fChannels[i]._bumpWindowLowEdge,fChannels[i]._bumpWindowHighEdge,zVal);
      }
   }
   if(fChannels.size()>1) {
      bool overflow(false);
      double zVal = GetZValue(fLocalPValue,overflow);
      if(overflow) {
         Info("Run","MultiChannel p-value (with OverlapFactor=%f): %f (> %f sigma)",fMultiChannelBumpOverlapFactor,fLocalPValue,zVal);
      } else {
         Info("Run","MultiChannel p-value (with OverlapFactor=%f): %f (%f sigma)",fMultiChannelBumpOverlapFactor,fLocalPValue,zVal);
      }
   }
   Double_t b1 = currChannel._bumpWindowLowEdge; Double_t b2 = currChannel._bumpWindowHighEdge;

   fBumpHunterStatisticPDF = new TH1D("bumpHunterPDF","BumpHunter Test Statistic PDF",100,0,tObs*1.5);
   fBumpHunterStatisticPDF->SetDirectory(0);
   fBumpHunterStatisticPDF->GetXaxis()->SetTitle("BumpHunter Statistic");
   fBumpHunterStatisticPDF->GetYaxis()->SetTitle("Probability Density");

   Int_t nDone = 0;//Int_t nGreater=0;

   fGlobalPValue = 1.;

   std::vector<Double_t> trial;
   std::vector<Double_t> pValues;std::vector<Double_t> pValuesLow; std::vector<Double_t> pValuesHigh;

   TH1D* nGreater = new TH1D("nGreater","",1,0,1);
   TH1D* nTotal = new TH1D("nTotal","",1,0,1);

   //if no nPseudo set, then use the pvalue to estimate the number of pseudo needed
   
   Int_t nPseudo = fnPseudo;
   if(nPseudo==0) {nPseudo = (fLocalPValue<0.0000001) ? 1000000 : 1./fLocalPValue; if(nPseudo<1000) nPseudo=1000;}

   Info("Run","Performing %d Pseudo-experiments....",nPseudo);
   std::cout << "|0%--------------------------------------------100%|" << std::endl;
   std::cout << "|" << std::flush;
   Int_t tickPoint = nPseudo/50;
   Int_t graphPoint = nPseudo/1000 + 1;
   while(/*!converged && */nDone < nPseudo) {
      if((nDone % tickPoint) ==0) std::cout << "*" << std::flush;
      Double_t tPseudo = EvaluateMultiChannelTestStatistic(true);
      fBumpHunterStatisticPDF->Fill(tPseudo);
      if(tPseudo>tObs) nGreater->Fill(0.5);//nGreater++;
      nDone++;nTotal->Fill(0.5);
      fGlobalPValue = (double)nGreater->GetEntries()/(double)nDone;
      if(nDone==nPseudo || nGreater->GetEntries()==1 || (nDone % graphPoint)==0) {
         pValues.push_back(fGlobalPValue);trial.push_back(nDone);
         TGraphAsymmErrors f; f.BayesDivide(nGreater,nTotal); 
         pValuesLow.push_back(f.GetErrorYlow(0));pValuesLow.push_back(f.GetErrorYhigh(0));
      }
   }
   std::cout << "|" << std::endl;


   if(fBumpHunterStatisticConvergenceGraph) delete fBumpHunterStatisticConvergenceGraph;
   fBumpHunterStatisticConvergenceGraph = new TGraphAsymmErrors(trial.size(),&trial[0],&pValues[0],&pValuesLow[0],&pValuesHigh[0]);

   delete nGreater; delete nTotal;


   currChannel._bumpWindowLowEdge=b1;currChannel._bumpWindowHighEdge=b2;


   if(fBumpHunterResultCanvas) delete fBumpHunterResultCanvas;

   std::cout << "blah" << std::endl;

      fBumpHunterResultCanvas = new TCanvas("bhCanvas","BumpHunter Results",500,600);
      fBumpHunterResultCanvas->Divide(1,3);
      fBumpHunterResultCanvas->cd(1);
      currChannel._histBack->Draw();currChannel._histBack->GetXaxis()->SetRangeUser(GetSearchLowEdge(),GetSearchHighEdge());
      currChannel._histData->Draw("same");
      TLine *bumpLineLeft = new TLine(currChannel._bumpWindowLowEdge,0,currChannel._bumpWindowLowEdge,currChannel._histBack->GetMaximum());
      TLine *bumpLineRight = new TLine(currChannel._bumpWindowHighEdge,0,currChannel._bumpWindowHighEdge,currChannel._histBack->GetMaximum());
      bumpLineLeft->SetLineColor(kRed);bumpLineLeft->SetLineStyle(2);bumpLineLeft->SetLineWidth(2);
      bumpLineRight->SetLineColor(kRed);bumpLineRight->SetLineStyle(2);bumpLineRight->SetLineWidth(2);
      bumpLineLeft->Draw("same");
      bumpLineRight->Draw("same");

      fBumpHunterResultCanvas->cd(2);
      fBumpHunterStatisticPDF->Scale(1./fBumpHunterStatisticPDF->GetSumOfWeights());
      fBumpHunterStatisticPDF->Draw();
      TLine *l1 = new TLine(tObs,0,tObs,fBumpHunterStatisticPDF->GetMaximum()/2.);
      l1->SetLineColor(kBlue);
      l1->SetLineWidth(2);
      l1->Draw("same");
      fBumpHunterResultCanvas->cd(3);
      fBumpHunterStatisticConvergenceGraph->Draw("ALP");
      fBumpHunterStatisticConvergenceGraph->GetHistogram()->GetXaxis()->SetTitle("# Pseudo-experiments");
      fBumpHunterStatisticConvergenceGraph->GetHistogram()->GetYaxis()->SetTitle("P-Value");
      TLine *l = new TLine(fBumpHunterStatisticConvergenceGraph->GetHistogram()->GetXaxis()->GetXmin(),pValues[pValues.size()-1],fBumpHunterStatisticConvergenceGraph->GetHistogram()->GetXaxis()->GetXmax(),pValues[pValues.size()-1]);
      l->SetLineColor(kRed);
      l->Draw("same");
      TLatex n;
      n.SetNDC();n.SetTextFont(43);n.SetTextSize(14);n.SetTextColor(kRed);
      std::stringstream m;m.precision(3);
      bool overflow(false);
      if(fGlobalPValue>0.) {
         double zVal = GetZValue(fGlobalPValue,overflow);
         m << zVal << "#sigma (local ";
         overflow=false;
         zVal = GetZValue(fLocalPValue,overflow);
         if(overflow) m << "> ";
         m << zVal << "#sigma)";
      } else {
         double zVal = GetZValue(1./double(nPseudo),overflow);
         m << " > " << zVal << "#sigma (local ";
         overflow=false;
         zVal = GetZValue(fLocalPValue,overflow);
         if(overflow) m << "> ";
         m << zVal << "#sigma)";
      }
      n.DrawLatex(0.5,0.6,m.str().c_str());

   if(fGlobalPValue>0.) {
      Info("Run","global p = %f (%f sigma)",fGlobalPValue,GetZValue(fGlobalPValue,overflow));
   } else {
      //none of the pseudodata was bigger than our pvalue, so estimate the global p-value as 1/nPseudo 
      Info("Run","global p < %f (%f sigma)",1./double(nPseudo),GetZValue(1./double(nPseudo),overflow));
   } 
   return 0;
}




Double_t TBumpHunter::GetSearchLowEdge() {
   Channel& currChannel = fChannels[fCurrentChannel];
   if(!currChannel._histBack) return 0;
   if(currChannel._searchLowEdge!=currChannel._searchLowEdge) {
      //look for first non-zero mc prediction bin
      int i=1;
      while(currChannel._histBack->GetBinContent(i)==0 && currChannel._histBack->GetBinError(i)==0 && i<=currChannel._histBack->GetNbinsX()) i++;
      if(i>currChannel._histBack->GetNbinsX()) {Error("GetSearchLowEdge","All bins are empty!?"); return 0;}
      return currChannel._histBack->GetBinLowEdge(i);
   }
   return currChannel._histBack->GetBinLowEdge(currChannel._histBack->FindFixBin(currChannel._searchLowEdge));
}

Double_t TBumpHunter::GetSearchHighEdge() {
   Channel& currChannel = fChannels[fCurrentChannel];
   if(!currChannel._histBack) return 0;
   if(currChannel._searchHighEdge!=currChannel._searchHighEdge) {
      //look for first zero mc prediction bin after low edge
      int i=currChannel._histBack->FindFixBin(GetSearchLowEdge());
      while(currChannel._histBack->GetBinContent(i)!=0 && i<=currChannel._histBack->GetNbinsX()) i++;
      if(i==currChannel._histBack->FindFixBin(GetSearchLowEdge())) {Error("GetSearchHighEdge","All bins are empty!?"); return 0;}
      return currChannel._histBack->GetBinLowEdge(i);
   }
   return currChannel._histBack->GetBinLowEdge(currChannel._histBack->FindFixBin(currChannel._searchHighEdge)+1);
}

Double_t TBumpHunter::GetMinWindowSize() {
   Channel& currChannel = fChannels[fCurrentChannel];
   if(!currChannel._histBack) return 0;
   //default is twice the average bin width
   if(fMinWindowSize!=fMinWindowSize) return 2*(currChannel._histBack->GetBinLowEdge(currChannel._histBack->GetNbinsX()+1) - currChannel._histBack->GetBinLowEdge(1))/currChannel._histBack->GetNbinsX();
   return fMinWindowSize;
}

Double_t TBumpHunter::GetMaxWindowSize() {
   Channel& currChannel = fChannels[fCurrentChannel];
   if(!currChannel._histBack) return 0;
   //default is half the search window
   if(fMaxWindowSize!=fMaxWindowSize) return (GetSearchHighEdge() - GetSearchLowEdge())/2.;
   return fMaxWindowSize;
}

Double_t TBumpHunter::GetWindowSizeStep() {
   Channel& currChannel = fChannels[fCurrentChannel];
   if(!currChannel._histBack) return 0;
   //default is the average bin width in the search region, or 1 bin if using bins for the window size
   if(fWindowSizeStep!=fWindowSizeStep) {
      if(kUseBinsForWindowSize) return 1; 
      else return (GetSearchHighEdge() - GetSearchLowEdge())/(currChannel._histBack->FindFixBin(GetSearchHighEdge())-currChannel._histBack->FindFixBin(GetSearchLowEdge()));
   }
   return fWindowSizeStep;
}

//decide what sequence to use for the various windows for the current settings
void TBumpHunter::EvaluateSearchPattern() {
   //find current hist in the list
   //Int_t index = fHistBacks->IndexOf(fHistBack);
   //if(index==-1) {Error("EvaluateSearchPattern","Something went badly wrong"); return;}

   Channel& currChannel = fChannels[fCurrentChannel];
   currChannel._centralWindows.clear();

   //fCentralWindows[index].clear();

   TH1* fHistBack = currChannel._histBack;


   Int_t startBin = fHistBack->FindFixBin(GetSearchLowEdge()); //first bin to use in search
   Int_t stopBin = fHistBack->FindFixBin(GetSearchHighEdge())-1; //last bin to use in search

   if(stopBin<startBin) return;

   //Double_t searchWidth = (GetSearchHighEdge() - GetSearchLowEdge());

   Double_t windowSizeStepSize = GetWindowSizeStep();
   Double_t currentWindowSize = GetMinWindowSize();
   Double_t maxWsize = GetMaxWindowSize();
   //Info("EvaluateSearchPattern","MinWindowSize=%f",currentWindowSize);

   bool hasGoodWindows(true);

   while(hasGoodWindows) {
      hasGoodWindows = false ;
      Int_t firstBin = startBin;
      Int_t lastBin = (kUseBinsForWindowSize) ? firstBin+(int)(currentWindowSize-0.5) : fHistBack->FindFixBin(fHistBack->GetBinLowEdge(firstBin)+currentWindowSize-0.00000001);
      while(lastBin <= stopBin) {
         //check the current bin isn't too wide 
         if((kUseBinsForMaxWindowSize && ((lastBin-firstBin+1)<=maxWsize)) || (!kUseBinsForMaxWindowSize && (fHistBack->GetBinLowEdge(lastBin+1)-fHistBack->GetBinLowEdge(firstBin))<=maxWsize)) {
            std::pair<Int_t,Int_t> p; p.first = firstBin; p.second=lastBin; currChannel._centralWindows.push_back(p);
            hasGoodWindows=true;
         }
         //Info("ee","%d to %d",p.first,p.second);
         if(kUseBinsForWindowSize) {
            //shift by 1 or half the current window size 
            if(currentWindowSize>1) firstBin += currentWindowSize/2.;
            else firstBin += 1;
         } else {
            //shift the first bin by 1 or half the current window size, whichever is bigger 
            if((currentWindowSize/2.-0.0000000001) > fHistBack->GetBinWidth(firstBin)) {
               firstBin = fHistBack->FindFixBin(fHistBack->GetBinLowEdge(firstBin)+(currentWindowSize/2.)-0.0000000001);
            } else {
               firstBin++;
            }
         }
         lastBin = (kUseBinsForWindowSize) ? firstBin+(int)(currentWindowSize-0.5) : fHistBack->FindFixBin(fHistBack->GetBinLowEdge(firstBin)+currentWindowSize-0.00000000001);
      }
      currentWindowSize += windowSizeStepSize;
   }
}

void TBumpHunter::PrintSearchPattern() {
   std::vector<std::pair<Int_t,Int_t> >& myWindows = fChannels[fCurrentChannel]._centralWindows;
   for(std::vector<std::pair<Int_t,Int_t> >::iterator it=myWindows.begin();it!=myWindows.end();++it) {
      Info("PrintSearchPattern","%d to %d [%f,%f]",it->first,it->second,fChannels[fCurrentChannel]._histBack->GetBinLowEdge(it->first),fChannels[fCurrentChannel]._histBack->GetBinLowEdge(it->second+1));
   }
}

Double_t TBumpHunter::EvaluateTestStatistic(TH1* data, Bool_t printOut) {
   //start at low edge
   //calculate pvalue in the given window
   //shift the window across max(first-bin-width in window,window/2)
   //recalculate pvalue
   //keep going until right edge of window is beyond high edge
   //increase window size by (highEdge-lowEdge)/(nBins in range from lowEdge to highEdge) - for equal bin widths this is just one bin
   //keep repeating all this until window size > maxWindowSize

   Double_t minPValue = 1.;

   Channel& currChannel = fChannels[fCurrentChannel];
   TH1* fHistBack = currChannel._histBack;

   std::vector<std::pair<Int_t,Int_t> >& myWindows = currChannel._centralWindows;

   for(std::vector<std::pair<Int_t,Int_t> >::iterator it=myWindows.begin();it!=myWindows.end();++it) {
      Double_t localPValue(1.);         //Double_t zVal(0.); //unimportant
      Double_t nObs = data->Integral(it->first,it->second);
      Double_t errExp(0.);
      Double_t nExp(0.);
      if(fBinModel==3) {
         //need to treat errors as correlated, so add up all the errors in the range 
         for(int i=it->first;i<=it->second;i++) {
            nExp += fHistBack->GetBinContent(i); errExp += fHistBack->GetBinError(i);
         }
      } else {
         nExp = fHistBack->IntegralAndError(it->first,it->second,errExp);
      }

      if(fTestStatisticType==BUMPHUNTER && nObs < nExp) localPValue = 1.; //dips are considered insignificant if only bumphunting
      else if(fTestStatisticType==DIPHUNTER && nObs > nExp) localPValue = 1.; //bumps are considered insigificant if diphunting
      else if(fBinModel==0) {
         localPValue = GetPoissonConvGammaPValue(nObs,nExp,errExp);
      } else if(fBinModel==1) {
         localPValue = GetPoissonPValue(nObs,nExp);
      } else if(fBinModel==2) {
         localPValue = GetGaussianPValue(nObs,nExp,errExp);
      } else if(fBinModel==3) {
         localPValue = GetPoissonConvGammaPValue(nObs,nExp,errExp);
      }
      if(printOut) Info("EvaluateTestStatistic","search region: [%f,%f] gave pvalue of %f (mc=%f, data=%f)", fHistBack->GetBinLowEdge(it->first),fHistBack->GetBinLowEdge(it->second+1),localPValue,nExp,nObs);
      if(localPValue < minPValue) {
        if(printOut) Info("EvaluateTestStatistic","smallest so far");
        minPValue = localPValue;
        currChannel._bumpWindowLowEdge = fHistBack->GetBinLowEdge(it->first);
        currChannel._bumpWindowHighEdge = fHistBack->GetBinLowEdge(it->second) + fHistBack->GetBinWidth(it->second);
        currChannel._bumpPvalue = localPValue;
      }
   }

   return -log(minPValue);
}

Double_t TBumpHunter::EvaluateMultiChannelTestStatistic(Bool_t generatePseudo, Bool_t printOut) {
   if(fChannels.size()==1) return EvaluateTestStatistic(((generatePseudo)? GenerateToyMC() : fChannels[0]._histData));
   //loop over channels, requiring the the worst bump overlap to still be better than the overlapfactor 
   std::vector<Double_t> bumpOverlaps;

   Double_t commonWindowLow = -DBL_MAX; Double_t commonWindowHigh = DBL_MAX;

   Double_t out(0.);

   for(unsigned int i=0;i<fChannels.size();i++) {
      SetChannel(i);
      Channel& currChannel = fChannels[fCurrentChannel];
      out += ((generatePseudo) ? EvaluateTestStatistic(GenerateToyMC()) :  EvaluateTestStatistic(currChannel._histData)); //can just add the nll of the pvalues to make the total test statistic
      bumpOverlaps.push_back(1.);
      //update the common window
      if(commonWindowLow<currChannel._bumpWindowLowEdge) commonWindowLow=currChannel._bumpWindowLowEdge;
      if(commonWindowHigh<currChannel._bumpWindowHighEdge) commonWindowHigh=currChannel._bumpWindowHighEdge;
      //calculate overlap factors for all bump windows so far 
      for(unsigned int j=0;j<bumpOverlaps.size();j++) {
         Double_t low = (commonWindowLow<fChannels[j]._bumpWindowLowEdge) ? fChannels[j]._bumpWindowLowEdge : commonWindowLow;
         Double_t high = (commonWindowHigh>fChannels[j]._bumpWindowHighEdge) ? fChannels[j]._bumpWindowHighEdge : commonWindowHigh;
         bumpOverlaps[j] = (high-low)/(fChannels[j]._bumpWindowHighEdge-fChannels[j]._bumpWindowLowEdge);
         if(printOut) Info("EvaluateMultiChannelTestStatistic","bump overlap (%d) = %f",j,bumpOverlaps[j]);
         if(bumpOverlaps[j]<fMultiChannelBumpOverlapFactor) return 0; //one bump not overlapping enough, so just return 0;
      }
   }
   return out;
}

TH1* TBumpHunter::GenerateToyMC() {
   return GeneratePseudoData(fChannels[fCurrentChannel]._histBack);
}

TH1* TBumpHunter::GeneratePseudoData(TH1* bgHist) {

  if(!bgHist) {
      Error("GeneratePseudoData","No background distribution given"); return 0;
  }

   //find in the channel list
//    Int_t index = fHistBacks->IndexOf(bgHist);
//    if(index==-1) {
//       Error("GeneratePseudoData","Background distribution is not in the list"); return 0;
//    }
// 
//    TH1* histData = (TH1*)(fHistPseudoDatas->At(index));

   TH1* histData = fChannels[fCurrentChannel]._histPseudoData;

   histData->Reset();

   //Double_t gausFactor = ROOT::Math::gaussian_cdf(pRand.Gaus(0,1),1); //where in the quantile spectrum to be for correlated uncertainties
   Double_t gausFactor = pRand.Uniform(); //choose quantile to use uniformly so that we fairly sample the possible gamma means

   Int_t startBin = bgHist->FindFixBin(GetSearchLowEdge()); //first bin to use in search
   Int_t stopBin = bgHist->FindFixBin(GetSearchHighEdge())-1; //last bin to use in search

   //loop over bins in the search region and fill with pseudo data
   for(Int_t j=startBin;j<=stopBin;j++) {
      Double_t mean = bgHist->GetBinContent(j);
      Double_t err = bgHist->GetBinError(j);

      if(mean==0 && err==0) continue; //leaves bin as 0 entries

      //define gamma parameters (a,b)
      Double_t b = mean/(err*err); // = E/V
      Double_t a = mean*b; // = E^2/V
      Double_t randomMean(0.); //used in binModel=0

      //which model are using 
      switch(fBinModel) {
         case 0: //poisson convoluted with a gamma distribution for the mean parameter
            //pick a random value for the mean
            randomMean = r.Gamma(a,1./b);
            histData->SetBinContent(j,pRand.PoissonD(randomMean));
            break;
         case 1: //poisson with no uncertainty on the mean (i.e. the bin error is meaningless 
            histData->SetBinContent(j,pRand.PoissonD(mean));
            break;
         case 2: //gaussian with mean and sigma as given 
            histData->SetBinContent(j,pRand.Gaus(mean,err));
         case 3: //poisson convoluted with a gamma on the mean, but all the errors are correlated
            //use the gausfactor to decide how far from the 0.5 (mean) position to stray 
            randomMean = ROOT::Math::gamma_quantile(gausFactor,a,1./b);
            histData->SetBinContent(j,pRand.PoissonD(randomMean));
            break;
      }
   }
   return histData;
}


Double_t TBumpHunter::GetPoissonPValue(Double_t nObs, Double_t E) {
   if(nObs>E) { //excess
      double p = ROOT::Math::inc_gamma_c(nObs,E);
      if(p==1. && nObs>100.) { //excess very excessive, and have a large nObs so no big problem summing one extra term
         return ROOT::Math::inc_gamma(nObs,E);
      } else {
         return 1-p;
      }
   } 
   return ROOT::Math::inc_gamma_c(nObs+1,E); //deficit
}

Double_t TBumpHunter::GetGaussianPValue(Double_t nObs, Double_t E, Double_t sigma) {
   if(nObs>E) return ROOT::Math::gaussian_cdf_c(nObs,sigma,E);
   return ROOT::Math::gaussian_cdf(nObs,sigma,E);
}


//return the p-value (and associated zValue = number of sigma) of seeing nObs "events"
//under the model a poisson distribution with the value of the mean parameter being itself uncertain
//. Assumption is that the mean parameter is distributed as a gamma density, with the expectation (mean) of the gamma function = E
// and the variance of the gamma density equal to err^2
Double_t TBumpHunter::GetPoissonConvGammaPValue(Double_t nObs, Double_t E, Double_t err) {
   //parameters (a,b) of gamma density can be written in terms of the expectation and variance:
   //Gammma(x, a,b);   - note this is not equal to gamma(x) or gamma(a,b), which are different functions
   double b = E/(err*err); // = E/V
   double a = E*b; // = E^2/V

   double pval = 0.;
   //decide if we should ignore systematics or not 
   if (a>100*nObs) {
      //stat error is big enough to ignore the syst error
      //considering only stat error means the p-value is given by:
      // (nData>nMC): pval = sum(x = nData->inf, Poisson(x,nMC)) = 1 - sum(x = 0->nData-1,Poisson(x,nMC))
      // (nData<=nMC): pval = sum(x = 0->nData, Poisson(x,nMC))
      // But sum(x = 0->nData-1,Poisson(x,nMC)) = gamma(nData,nMC)/gamma(nData); <---- see maths websites
      // so we have:
      // (nData>nMC): pval = 1 - gamma(nData,nMC)/gamma(nData);
      // (nData<=nMC): pval = gamma(nData+1,nMC)/gamma(nData);
      // .....And ROOT provides gamma(a,b)/gamma(a) = ROOT::Math::inc_gamma_c(a,b)
      //pval = (nObs<=E) ? ROOT::Math::inc_gamma_c(nObs+1,E) : (1. - ROOT::Math::inc_gamma_c(nObs,E));
      pval = GetPoissonPValue(nObs,E);
   } else {
      //use recursive formula to solve:
      // (nData>nMC): pval = 1 - sum(x=0->nData-1, Integral(y=0->inf, Poisson(x,y)Gamma(y,a,b) dy))
      // (nData<=nMC): pval = sum(x=0->nData, Integral(y=0->inf, Poisson(x,y)Gamma(y,a,b) dy))
      //i.e. integrating out the unknown parameter y
      // Recursive formula: P(n;A,B) = P(n-1;A,B) (A+n-1)/(n*(1+B))
      unsigned stop=nObs;
      if (nObs>E) --stop;
      double sum = 0;
      if(a>100) {
         /// NB: must work in log-scale otherwise troubles!
         double logProb = a*log(b/(1+b));
         sum=exp(logProb); // P(n=0)
         for (unsigned u=1; u<=stop; ++u) {
            logProb += log((a+u-1)/(u*(1+b)));
            sum += exp(logProb);
         }
      } else {
         double p0 = pow(b/(1+b),a); // P(0;A,B)
         double pLast = p0;
         sum = p0;
         for (unsigned k=1; k<=stop; ++k) {
            double p = pLast * (a+k-1) / (k*(1+b));
            sum += p;pLast = p;
         }
      }
      pval = (nObs>E) ?  1-sum : sum; 
   } 

   //bool overflow(false);
   //zValue = GetZValue(pValue,overflow) //large z-values correspond to small p-values.... i.e. significant diff
   //zValue = (nObs<E) ? -1.*zValue : zValue; //flip the z-values of deficits

   return pval;
}

Double_t TBumpHunter::GetZValue(Double_t pValue, bool& overflow) {
   if(pValue<0.000000001) { overflow=true; return 6.; }
   return sqrt(2.)*TMath::ErfInverse(1.-2.*pValue);
}