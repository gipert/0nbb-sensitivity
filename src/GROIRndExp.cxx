/* GROIRndExp.cxx
 *
 * Author: Luigi Pertoldi - pertoldi@pd.infn.it
 * Created: 11 Apr 2018
 *
 */

#include "GROIRndExp.h"

#include "TF1.h"

GROIRndExp::GROIRndExp(double exposure, double BI, double halfLife, double ROIWidth, int binning, double FWHM) :
    TH1D(),
    fExposure(exposure),
    fBI(BI),
    fHalfLife(halfLife),
    fROIWidth(ROIWidth),
    fBinning(binning),
    fFWHM(FWHM)
{
    // initialise base object
    this->SetName("ROIRndExp");
    this->SetTitle("Random experiment in ROI");
    this->SetBins(fROIWidth/fBinning, fQbb-fROIWidth/2, fQbb+fROIWidth/2);

    // calculate total number of bkg and signal counts
    int B = this->GetBkgCounts();
    int S = this->GetSignalCounts();

    // fill
    TF1 sModel("signal", "gaus", fQbb-fROIWidth/2, fQbb+fROIWidth/2);
    sModel.SetParameters(1, fQbb, fFWHM*0.4246);

    this->FillRandom("signal", S);
    this->FillRandom("pol0",   B);
}

int GROIRndExp::GetBkgCounts() const {
    return (int)std::round(fBI*fExposure*fROIWidth);
}

int GROIRndExp::GetSignalCounts() const {
    return fHalfLife != 0 ? (int)std::round(4.1615E24*fExposure/fHalfLife) : 0;
}
