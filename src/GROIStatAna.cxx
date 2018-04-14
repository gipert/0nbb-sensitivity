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
    int B = fExp->GetBkgCounts();
    int S = fExp->GetSignalCounts();

    // define parameters
                   this->AddParameter("B", 0, 2*B);
    if (hasSignal) this->AddParameter("S", 0, 2*S);

    // set priors
    auto funcStr = "(x>0)*(exp(-0.5*((x-" + std::to_string(B) + ")/" + std::to_string(B*.1/2) + ")**2))";
    TF1 posgaus("posgaus", funcStr.c_str(), 0, 1);

                   this->SetPrior(0, posgaus);
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
