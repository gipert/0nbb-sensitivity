/* discSensVsBI.cxx
 *
 * Author: Luigi Pertoldi - pertoldi@pd.infn.it
 * Create 12 Apr 2018
 *
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <chrono>

#ifdef GOPARALLEL
    #include "omp.h"
#endif

#include "GROIRndExp.h"
#include "GROIStatAna.h"

// BAT and ROOT
#include <BAT/BCModelManager.h>

// other
//#include "tools/jsoncpp/json/json.h"
#include "json/json.h"
#include "ProgressBar.h"

struct BF {
    double bayes_factor;
    double error;
};

BF GetBayesFactor(double BI, double hl, Json::Value J);
void ConfigureIntegrationModel(BCModel* m, Json::Value JSONIntConf);

int main(int argc, char** argv) {

#ifdef GOPARALLEL
    ROOT::EnableThreadSafety();
#endif

    if (argc != 3) {
        std::cerr << "USAGE: discSensVsBI <masterconf.json> <BI-value>\n";
        return 1;
    }

    std::string filename = argv[1];
    std::ifstream fJSON(filename.c_str());
    if (!fJSON.is_open()) {
        std::cerr << filename << " does not exists!\n";
        return 1;
    }

    double BI = std::atof(argv[argc-1]);

    Json::Value J;
    fJSON >> J;

    bool verbose = false;
    bool pbar = true;
    if (J["verbose"]) verbose = J["verbose"].asBool();
    if (J["progress-bar"]) pbar = J["progress-bar"].asBool();

    // set BAT output
    std::map<std::string,BCLog::LogLevel> mLog = {
        {"debug",   BCLog::debug},
        {"detail",  BCLog::detail},
        {"summary", BCLog::summary},
        {"warning", BCLog::warning},
        {"nothing", BCLog::nothing}
    };
    if (J["bat-verbose"]) BCLog::SetLogLevel(mLog[J["bat-verbose"].asString()]);

    int    nexp     = J["nexperiments"].asInt();
    double thBF     = J["threshold-bayesfactor"].asDouble();
    double eps      = J["root-search-precision"].asDouble();

    using clock  = std::chrono::steady_clock;
    using t_unit = std::chrono::seconds;

    if (verbose) std::cout << "BI = " << BI << " cts/(keV•kg•yr)\n" << std::flush;
    // search strategy: rude Bisection Method
    // initialise search boundaries for sensitivity
    double hl_low    = J["0nbb-halflife-range"][0].asDouble(); // x1
    double hl_mid    = 0;
    double hl_up     = J["0nbb-halflife-range"][1].asDouble(); // x2
    int    nsucc_low = 0;                                      // f(x1)
    int    nsucc_mid = 0;
    int    nsucc_up  = 0;                                      // f(x2)

    ProgressBar bar(nexp);
    if (!pbar) bar.ShowBar(false);
    // we need a first value for a boundary, e.g. hl_low
    if (verbose) std::cout << "Calculate starting point... ";
#pragma omp parallel for reduction(+:nsucc_low)
    for (int i = 0; i < nexp; ++i) {
#pragma omp critical
        if (verbose) bar.Update();
        if (GetBayesFactor(BI, hl_low, J).bayes_factor >= thBF) nsucc_low++;
    }
    while (true) {
        std::cout << "\nLooking into [" << hl_low << "," << hl_up << "]yr, ∆ = " << hl_up-hl_low << "yr, with " << nexp << " experiments...\n";
        // our next test point
        hl_mid = (hl_low+hl_up)/2;

        clock::time_point begin = clock::now();
        // generate the experiments
        bar.Reset();
#pragma omp parallel for reduction(+:nsucc_up)
        for (int i = 0; i < nexp; ++i) {
#pragma omp critical
            if (verbose) bar.Update();
            if (GetBayesFactor(BI, hl_mid, J).bayes_factor >= thBF) nsucc_mid++;
        }
        clock::time_point end = clock::now();
        if (verbose) std::cout << " time: " << std::chrono::duration_cast<t_unit>(end-begin).count() << "s\n";
        // determine direction of next search
        if ((nsucc_low-nexp/2)*(nsucc_mid-nexp/2) > 0
            and (fabs(nsucc_up-nsucc_low) > 2*eps
                 or fabs(nsucc_mid-nexp/2) > eps)) {

            if (verbose) {
                std::cout << "-----------------------------------------\n"
                          << "\thalf-life\t# succ. exp. / tot\n"
                          << "-----------------------------------------\n"
                          << "low\t" << std::setw(10) << hl_low << '\t' << nsucc_low << '\n'
                          << "mid\t" << std::setw(10) << hl_mid << '\t' << nsucc_mid << " <-- new low\n"
                          << "up\t"  << std::setw(10) << hl_up  << '\t' << nsucc_up  << " <-- new up\n"
                          << "-----------------------------------------\n";
            }

            hl_low = hl_mid;
            // hl_up = hl_up;
            nsucc_low = nsucc_mid;
        }
        else if ((nsucc_low-nexp/2)*(nsucc_mid-nexp/2) < 0
                  and (fabs(nsucc_up-nsucc_low) > 2*eps
                       or fabs(nsucc_mid-nexp/2) > eps)) {

            if (verbose) {
                std::cout << "-----------------------------------------\n"
                          << "\thalf-life\t# succ. exp. / tot\n"
                          << "-----------------------------------------\n"
                          << "low\t" << std::setw(10) << hl_low << '\t' << nsucc_low << " <-- new low\n"
                          << "mid\t" << std::setw(10) << hl_mid << '\t' << nsucc_mid << " <-- new up\n"
                          << "up\t"  << std::setw(10) << hl_up  << '\t' << nsucc_up  << "\n"
                          << "-----------------------------------------\n";
            }

            // hl_low = hl_low
            hl_up = hl_mid;
            // nsucc_low = nsucc_low;
            nsucc_up = nsucc_mid;
        }
        else {
            if (hl_low == J["0nbb-halflife-range"][0].asDouble()
                or hl_up == J["0nbb-halflife-range"][1].asDouble()) {
                std::cerr << "Warning: range boundaries reached!\n";
            }
            break;
        }
        nsucc_mid = 0;
    }
    if (verbose) std::cout << "Found sensitivity -> " << (hl_low+hl_up)/2 << " +- " << (hl_up-hl_low)/2 << "yr\n";
    std::ofstream out;
    out.open("results.txt", std::ofstream::out | std::ofstream::app);
    out << BI << '\t' << (hl_low+hl_up)/2 << '\t' << (hl_up-hl_low)/2 << '\n';
    out.close();
    return 0;
}

BF GetBayesFactor(double BI, double hl, Json::Value J) {
    // generate experiment
    GROIRndExp rndExp(J["exposure"].asDouble(), BI, hl);

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

    // set integration methods
    using IM  = BCIntegrate::BCIntegrationMethod;
    using CIM = BCIntegrate::BCCubaMethod;
    std::map<std::string,IM> iMeth = {
        {"kIntMonteCarlo", IM::kIntMonteCarlo},
        {"kIntGrid",       IM::kIntGrid},
        {"kIntLaplace",    IM::kIntLaplace},
        {"kIntCuba",       IM::kIntCuba}
    };
    std::map<std::string,CIM> iCMeth = {
        {"kCubaDivonne", CIM::kCubaDivonne},
        {"kCubaVegas",   CIM::kCubaVegas},
        {"kCubaSuave",   CIM::kCubaSuave},
        {"kCubaCuhre",   CIM::kCubaCuhre}
    };
    mBS->SetIntegrationMethod(iMeth[J["BS-model"]["integration-method"].asString()]);
    mB ->SetIntegrationMethod(iMeth[J["B-model"]["integration-method"].asString()]);
    if (J["BS-model"]["cuba-integration-method"]) {
        mBS->SetCubaIntegrationMethod(iCMeth[J["BS-model"]["cuba-integration-method"].asString()]);
    }
    if (J["B-model"]["cuba-integration-method"]) {
        mB ->SetCubaIntegrationMethod(iCMeth[J["B-model"]["cuba-integration-method"].asString()]);
    }

    // set integration settings, if available
    ConfigureIntegrationModel(mBS, J["BS-model"]);
    ConfigureIntegrationModel(mB ,  J["B-model"]);

    mgr.Integrate();
/*
    std::cout << "B+S model Integral = "
        << mBS->GetIntegral() << " +- "
        << 100*mBS->GetError()/mBS->GetIntegral() << "%\n";
    std::cout << "B   model Integral = "
        << mB ->GetIntegral() << " +- "
        << 100*mB ->GetError()/mB ->GetIntegral() << "%\n";
    std::cout << "Bayes Factor = " + std::to_string(mgr.BayesFactor(0, 1)) << std::endl;
*/
    double bf   = mgr.BayesFactor(0,1);
    double i0   = mBS->GetIntegral();
    double s_i0 = mBS->GetError();
    double i1   = mB ->GetIntegral();
    double s_i1 = mB ->GetError();

    if (i0 <= 0 or i1 <= 0) std::cerr << "GetBayesFactor: Warning: Integral = 0!\n";

    return {bf, bf*sqrt(s_i0*s_i0/(i0*i0) + s_i1*s_i1/(i1*i1))};
}

void ConfigureIntegrationModel(BCModel* m, Json::Value JSONModelSettings) {
    auto& MS = JSONModelSettings;
    auto& IS = MS["integrator-settings"];

    // use N iterations (min/max) from chosen integration method
    if (MS["integration-method"].asString() != "kIntCuba") {
        if (IS[MS["integration-method"].asString()]["niter-max"]) {
            m->SetNIterationsMax(IS[MS["integration-method"].asString()]["niter-max"].asInt64());
        }
        if (IS[MS["integration-method"].asString()]["niter-min"]) {
            m->SetNIterationsMin(IS[MS["integration-method"].asString()]["niter-min"].asInt64());
        }
    }
    else {
        if (IS["kIntCuba"][MS["cuba-integration-method"].asString()]["niter-max"]) {
            m->SetNIterationsMax(IS["kIntCuba"][MS["cuba-integration-method"].asString()]["niter-max"].asInt64());
        }
        if (IS["kIntCuba"][MS["cuba-integration-method"].asString()]["niter-min"]) {
            m->SetNIterationsMin(IS["kIntCuba"][MS["cuba-integration-method"].asString()]["niter-min"].asInt64());
        }
    }

    if (IS["kIntCuba"]) {
        auto& j = IS["kIntCuba"];
        if (j["kCubaVegas"]) {
            auto& ji = j["kCubaVegas"];
            auto o = m->GetCubaVegasOptions();

            if (ji["flags"])   o.flags = ji["flags"].asInt();
            if (ji["nstart"]) o.nstart = ji["nstart"].asInt64();

            m->SetCubaOptions(o);
        }
        if (j["kCubaSuave"]) {
            auto& ji = j["kCubaSuave"];
            auto o = m->GetCubaSuaveOptions();

            if (ji["flags"]) o.flags = ji["flags"].asInt();
            if (ji["neval"]) o.neval = ji["neval"].asInt64();

            m->SetCubaOptions(o);
        }
        if (j["kCubaDivonne"]) {
            auto& ji = j["kCubaDivonne"];
            auto o = m->GetCubaDivonneOptions();

            if (ji["flags"]) o.flags = ji["flags"].asInt();

            m->SetCubaOptions(o);
        }
    }
    return;
}
