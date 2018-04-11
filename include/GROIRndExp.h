/* GROIRndExp.h
 *
 * Author: Luigi Pertoldi - pertoldi@pd.infn.it
 * Created: 11 Apr 2018
 *
 */

#ifndef ZERONBBSENSITIVITY_GROIRNDEXP_H
#define ZERONBBSENSITIVITY_GROIRNDEXP_H

#include <TH1D.h>

class GROIRndExp {

    public :
        // exposure in kg•yr
        // BI in cts/(keV•kg•yr)
        // half-live in yr
        GROIRndExp(double exposure,
                   double BI,
                   double halfLife,
                   double ROIWidth = 100,
                   int binning = 1,
                   double FWHM = 3);
        ~GROIRndExp();

    private :
        const double fQbb = 2039; // keV
        const int    fROIWidth;   // keV
        const int    fBinning;    // keV
        const double fFWHM;       // keV

        TH1D fSpectrum;
};

#endif
