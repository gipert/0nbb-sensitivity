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
    /* [0] */  this->AddParameter("B", 0, 10);
    /* [1] */  this->AddParameter("S", 0, 10);

    auto B = fExp->GetBkgCounts();

    // set priors
    auto funcStr = "(x>0)*(exp(-0.5*((x-" + std::to_string(B) + ")/" + std::to_string(B*.1/2) + ")**2))";
    TF1 posgaus("posgaus", funcStr.c_str(), 0, 1);
    this->SetPrior(0, posgaus);
    this->SetPriorConstant(1);

    // eventually fix signal to zero
    if (!hasSignal) this->GetParameter(1).Fix(0);
}

double GROIStatAna::LogLikelihood(const std::vector<double> & p) {
    double logprob = 0.;

    int    nBins    = fExp->GetNbinsX();
    double qbb      = fExp->GetQbb();
    double sigma    = fExp->GetSigmaRes();

    for (int i = 0; i < nBins; ++i) {
        logprob += BCMath::LogPoisson(
                       fExp->GetBinContent(i),
                       p[0]/nBins + p[1]*TMath::Gaus(fExp->GetBinCenter(i), qbb, sigma, true)
                   );
    }
    return logprob;
}
