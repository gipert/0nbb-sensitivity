/* GROIStatAna.cxx
 *
 * Author: Luigi Pertoldi - pertoldi@pd.infn.it
 * Created: 11 Apr 2018
 *
 */

#include "GROIStatAna.h"

#include "TMath.h"
#include "TF1.h"

#include <BAT/BCMath.h>

GROIStatAna::GROIStatAna(GROIRndExp* initExp, bool hasSignal, std::string name) :
    BCModel(name.c_str()),
    fHasSignal(hasSignal),
    fExp(initExp)
{
    // define parameters
    double B = fExp->GetExpectedBkgCounts();
    // bkg range is kinda arbitrary, in principle should be [0,+∞]
    // given the prior, 4B means 8σ
                   this->AddParameter("B", 0, 4*B);
    // signal range goes from 0 to S_max, given by the current most stringent lower limit
    if (hasSignal) this->AddParameter("S", 0, 4.1615E24*fExp->GetExposure()/currentT12lowerLimit);

    // set priors
    // gaussian if x>0
    // flat     if x<0
    TF1 posgaus("posgaus", "gaus", 0, 4*B);
    posgaus.SetParameters(1, B, B/2);

                   this->SetPrior(0, posgaus);
    // flat prior for signal
    if (hasSignal) this->SetPriorConstant(1);
}

double GROIStatAna::LogLikelihood(const std::vector<double> & p) {
    double logprob = 0.;
    double exp = 0.;

    int    nBins    = fExp->GetNbinsX();
    double qbb      = fExp->GetQbb();
    double sigma    = fExp->GetSigmaRes();

    for (int i = 0; i < nBins; ++i) {
        if (fHasSignal) exp = p[0]/nBins + p[1]*TMath::Gaus(fExp->GetBinCenter(i), qbb, sigma, true);
        else exp = p[0]/nBins;
        logprob += BCMath::LogPoisson(fExp->GetBinContent(i), exp);
    }
    return logprob;
}
