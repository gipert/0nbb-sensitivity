/* GROIRndExp.h
 *
 * Author: Luigi Pertoldi - pertoldi@pd.infn.it
 * Created: 11 Apr 2018
 *
 */

#ifndef ZERONBBSENSITIVITY_GROIRNDEXP_H
#define ZERONBBSENSITIVITY_GROIRNDEXP_H

#include <TH1D.h>

class GROIRndExp : public TH1D {

    public :
        // delete default constructor
        GROIRndExp()                             = delete;
        ~GROIRndExp()                            = default;
        // delete copy constructor/assignement
        GROIRndExp           (GROIRndExp const&) = delete;
        GROIRndExp& operator=(GROIRndExp const&) = delete;
        // default move constructor/assignement, redundand
        GROIRndExp           (GROIRndExp&&)      = default;
        GROIRndExp& operator=(GROIRndExp&&)      = default;

        GROIRndExp(double exposure,
                   double BI,
                   double halfLife,
                   double ROIWidth = 100,
                   int    binning = 1,
                   double FWHM = 3);

        // getters
        double GetExposure() const {return fExposure;}
        double GetBI()       const {return fBI;}
        double GetHalfLife() const {return fHalfLife;}
        double GetQbb()      const {return fQbb;}
        double GetROIWidth() const {return fROIWidth;}
        double GetFWHM()     const {return fFWHM;}
        double GetSigmaRes() const {return fFWHM*0.4246;}

        int GetBkgCounts()    const {return (int)std::round(fBI*fExposure*fROIWidth);}
        int GetSignalCounts() const {return (int)std::round(4.1615E24*fExposure/fHalfLife);}

    private :
        const double fExposure;   // kg•yr
        const double fBI;         // cts/(keV•kg•yr)
        const double fHalfLife;   // yr
        const double fQbb = 2039; // keV
        const int    fROIWidth;   // keV
        const int    fBinning;    // keV
        const double fFWHM;       // keV
};

#endif
