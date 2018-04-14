/* discSensVsBI.cxx
 *
 * Author: Luigi Pertoldi - pertoldi@pd.infn.it
 * Create 12 Apr 2018
 *
 */
#include <iostream>
#include <fstream>
#include <map>

#include "GROIRndExp.h"
#include "GROIStatAna.h"

// BAT
#include <BAT/BCModelManager.h>

// other
//#include "tools/jsoncpp/json/json.h"
#include "json/json.h"

#include "TApplication.h"

int main(int argc, char** argv) {

    std::string filename;
    if (argc == 2) filename = argv[1];
    else {
        std::cerr << "I need the masterconf!\n";
        return 1;
    }

    std::ifstream fJSON(filename.c_str());
    if (!fJSON.is_open()) {
        std::cerr << filename << " does not exists!\n";
        return 1;
    }

    Json::Value J;
    fJSON >> J;

    // set silent BAT output
    BCLog::SetLogLevel(BCLog::nothing);

    // generate experiment
    GROIRndExp rndExp(J["exposure"].asDouble(),
                      J["BI"].asDouble(),
                      J["0nbb-halflife"].asDouble());

    // create B+S and B models
    GROIStatAna anaBS(&rndExp, true,  "B+S model");
    GROIStatAna anaB (&rndExp, false, "B model");
    // create model manager
    BCModelManager mgr;
    mgr.AddModel(&anaBS, J["BS-model"]["prior-probability"].asDouble());
    mgr.AddModel(&anaB ,  J["B-model"]["prior-probability"].asDouble());
    // comfy aliases
    auto mBS = mgr.GetModel(0);
    auto mB  = mgr.GetModel(1);

    // set integration settings, if available
    // set integration methods
    using IM = BCModel::BCIntegrationMethod;
    std::map<std::string,IM> iMeth = {
        {"kIntMonteCarlo", IM::kIntMonteCarlo},
        {"kIntCuba",       IM::kIntCuba}
    };
    mBS->SetIntegrationMethod(iMeth[J["BS-model"]["integration-method"].asString()]);
    mB ->SetIntegrationMethod(iMeth[ J["B-model"]["integration-method"].asString()]);

    auto ConfigureIntegrationModel = [&](BCModel* m, Json::Value JSONIntConf) {
        auto& IS = JSONIntConf;
        if (IS["kIntMonteCarlo"]) {
            auto& j = IS["kIntMonteCarlo"];
            if (j["niter-max"]) m->SetNIterationsMax(j["niter-max"].asInt64());
        }
        if (IS["kIntCuba"]) {
            auto& j = IS["kIntCuba"];
            if (j["CubaVegas"]) {
                auto& ji = j["CubaVegas"];
                auto o = mBS->GetCubaVegasOptions();
                if (ji["nstart"]) o.nstart = ji["nstart"].asInt64();
                m->SetCubaOptions(o);
            }
            if (j["CubaSuave"]) {
                auto& ji = j["CubaSuave"];
                auto o = mBS->GetCubaSuaveOptions();
                if (ji["neval"]) o.neval = ji["neval"].asInt64();
                m->SetCubaOptions(o);
            }
        }
    };

    // integration settings
    ConfigureIntegrationModel(mBS, J["BS-model"]["integrator-settings"]);
    ConfigureIntegrationModel(mB ,  J["B-model"]["integrator-settings"]);

    mgr.Integrate();

    std::cout << "B+S model Integral = "
              << mBS->GetIntegral() << " +- "
              << 100*mBS->GetError()/mBS->GetIntegral() << "%\n";
    std::cout << "B   model Integral = "
              << mB ->GetIntegral() << " +- "
              << 100*mB ->GetError()/mB ->GetIntegral() << "%\n";
    std::cout << "Bayes Factor = " + std::to_string(mgr.BayesFactor(0, 1)) << std::endl;

    return 0;
}
