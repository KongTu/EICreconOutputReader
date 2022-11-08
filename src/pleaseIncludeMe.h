#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>

#include "ROOT/RDataFrame.hxx"
#include <TH1D.h>
#include <TFitResult.h>
#include <TRandom3.h>
#include <TCanvas.h>

#include "TFile.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TVector2.h"
#include "TVector3.h"

#include "fmt/color.h"
#include "fmt/core.h"

// #include "nlohmann/json.hpp"
#include "edm4eic/InclusiveKinematicsData.h"
#include "edm4eic/ReconstructedParticleData.h"
// #include "edm4eic/ReconstructedParticleCollection.h"

#define PI            3.1415926
#define MASS_ELECTRON 0.00051
#define MASS_PROTON   0.93827
#define MASS_PION     0.13957
#define MASS_KAON     0.493667
#define MASS_AU197    183.45406466643374

auto getNtrk(const std::vector<edm4eic::ReconstructedParticleData>& parts)
{
  std::vector<int> mult;
  int n=0;
  for(auto& i1 : parts){
    if(i1.charge!=0) n++;
  }
  mult.push_back( n );
  return mult;
}

auto getEta(const std::vector<edm4eic::ReconstructedParticleData>& parts)
{
	std::vector<double> eta;
	for(auto& i1 : parts){
		TVector3 trk(i1.momentum.x,i1.momentum.y,i1.momentum.z);
		if(i1.charge!=0) eta.push_back(trk.Eta());
	}
	return eta;
}
