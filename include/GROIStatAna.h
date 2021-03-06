/* GROIStatAna.h
 *
 * Author: Luigi Pertoldi - pertoldi@pd.infn.it
 * Created: 11 Apr 2018
 *
 */

#ifndef ZERONBBSENSITIVITY_GROISTATANA_H
#define ZERONBBSENSITIVITY_GROISTATANA_H

#include <BAT/BCModel.h>
#include "GROIRndExp.h"

class GROIStatAna : public BCModel {

    public :
        // delete default constructor
        GROIStatAna()                              = delete;
        ~GROIStatAna()                             = default;
        // delete copy constructor/assignement
        GROIStatAna           (GROIStatAna const&) = delete;
        GROIStatAna& operator=(GROIStatAna const&) = delete;
        // default move constructor/assignement, redundand
        GROIStatAna           (GROIStatAna&&)      = default;
        GROIStatAna& operator=(GROIStatAna&&)      = default;

        GROIStatAna(GROIRndExp* initExp, bool hasSignal = true, std::string name = "0nbbStatAna");
        // setters
//        void SetSpectrum(GROIRndExp* initExp) {fExp = initExp;}

        // getters
        const GROIRndExp* GetRndExp() const {return fExp;}
        bool HasSignal() const {return fHasSignal;}

        // methods from BCModel to be overloaded
        double LogLikelihood(const std::vector<double>& parameters);

    private :
        const bool fHasSignal;
        const double currentT12lowerLimit = 8.0E25;
        const GROIRndExp* fExp;
};

#endif
