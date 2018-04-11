/* GROIRndExp.cxx
 *
 * Author: Luigi Pertoldi - pertoldi@pd.infn.it
 * Created: 11 Apr 2018
 *
 */

#include "GROIRndExp.h"

#include "TF1.h"

GROIRndExp::GROIRndExp(double exposure, double BI, double halfLife, double ROIWidth, int binning, double FWHM) :
    fROIWidth(ROIWidth),
    fBinning(binning),
    fFWHM(FWHM),
    fSpectrum("ROIRndExp",
              "Random experiment in ROI",
              fROIWidth/fBinning,
              fQbb-fROIWidth/2,
              fQbb+fROIWidth/2)
{
    int B = BI * exposure * fROIWidth;
    int S = 4.1615E24 * exposure / halfLife;

    TF1 sModel("signal", "gaus", fQbb-fROIWidth/2, fQbb+fROIWidth/2);
    sModel.SetParameters(fQbb, fFWHM*0.4246);

    fSpectrum.FillRandom("signal", S);
    fSpectrum.FillRandom("pol0",   B);
}
