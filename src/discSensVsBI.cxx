/* discSensVsBI.cxx
 *
 * Author: Luigi Pertoldi - pertoldi@pd.infn.it
 * Create 12 Apr 2018
 *
 */
#include <iostream>
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
#include "TThread.h"

// other
//#include "tools/jsoncpp/json/json.h"
#include "json/json.h"

double GetBayesFactor(double BI, double hl, Json::Value J);
void ConfigureIntegrationModel(BCModel* m, Json::Value JSONIntConf);

int main(int argc, char** argv) {

#ifdef GOPARALLEL
    TThread::Initialize();
#endif

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

    bool verbose = false;
    if (J["verbose"]) verbose = J["verbose"].asBool();

    // set BAT output
    std::map<std::string,BCLog::LogLevel> mLog = {
        {"debug",   BCLog::debug},
        {"detail",  BCLog::detail},
        {"summary", BCLog::summary},
        {"warning", BCLog::warning},
        {"nothing", BCLog::nothing}
    };
    if (J["bat-verbose"]) BCLog::SetLogLevel(mLog[J["bat-verbose"].asString()]);

    double BImin    = J["BI-range"][0].asDouble();
    double BImax    = J["BI-range"][1].asDouble();
    int    BIpoints = J["BI-npoints"].asInt();
    int    nexp     = J["nexperiments"].asInt();
    double thBF     = J["threshold-bayesfactor"].asDouble();
    double eps      = J["root-search-precision"].asDouble();

    std::ofstream outfile("results.txt");

    using clock  = std::chrono::steady_clock;
    using t_unit = std::chrono::microseconds;

    // main loop over x-axis values: BIs
#pragma omp parallel for schedule(auto)
    for (int j = 0; j < BIpoints; ++j) {
        double BI = BImin + (BImax-BImin)*j/BIpoints;
//    for (double BI = BImin; BI <= BImin; BI += (BImax-BImin)/BIpoints) {
#ifdef GOPARALLEL
#pragma omp critical
{
        if (verbose) std::cout << "Thread (" << omp_get_thread_num() << ") "
                               << "processing x = " << BI << " cts/(keV•kg•yr)\n" << std::flush;
}
#else
        if (verbose) std::cout << "processing x = " << BI << " cts/(keV•kg•yr)\n" << std::flush;
#endif
        // search strategy: rude Bisection Method
        // initialise search boundaries for sensitivity
        double hl_low    = J["0nbb-halflife-range"][0].asDouble(); // x1
        double hl_up     = J["0nbb-halflife-range"][1].asDouble(); // x2
        int    nsucc_low = 0;                                      // f(x1)
        int    nsucc_up  = 0;                                      // f(x2)
        // we need a first value for a boundary, e.g. hl_low
        for (int i = 0; i < nexp; ++i) {
            if (GetBayesFactor(BI, hl_low, J) >= thBF) nsucc_low++;
        }
        while (true) {
#ifdef GOPARALLEL
            if (verbose and omp_get_num_threads() == 1) std::cout << "Looking into [" << hl_low << "," << hl_up << "] ";
#else
            if (verbose) std::cout << "Looking into [" << hl_low << "," << hl_up << "] ";
#endif
            // our next test point
            auto hl_mid = (hl_low+hl_up)/2;

            clock::time_point begin = clock::now();
            // generate the experiments
            for (int i = 0; i < nexp; ++i) {
                if (GetBayesFactor(BI, hl_mid, J) >= thBF) nsucc_up++;
            }
            clock::time_point end = clock::now();
#ifdef GOPARALLEL
            if (verbose and omp_get_num_threads() == 1) std::cout << std::chrono::duration_cast<t_unit>(end-begin).count() << " s\n";
#else
            if (verbose) std::cout << std::chrono::duration_cast<t_unit>(end-begin).count() << " s\n";
#endif
            // determine direction of next search
            if ((nsucc_low-nexp/2)*(nsucc_up-nexp/2) > 0 and fabs(nsucc_low-nexp/2) > eps) {
                hl_low = hl_mid;
                // hl_up = hl_up;
                nsucc_low = nsucc_up;
            }
            else if ((nsucc_low-nexp/2)*(nsucc_up-nexp/2) < 0 and fabs(nsucc_low-nexp/2) > eps) {
                // hl_low = hl_low
                hl_up = hl_mid;
                // nsucc_low = nsucc_low;
            }
            else {
                if (hl_low == J["0nbb-halflife-range"][0].asDouble()
                    or hl_up == J["0nbb-halflife-range"][1].asDouble()) {
                    std::cout << "Warning: range boundaries reached!\n";
                }
                break;
            }
            nsucc_up = 0;
        }
#ifdef GOPARALLEL
#pragma omp critical
{
        if (verbose) std::cout << "Thread (" << omp_get_thread_num() << ") "
                               << "found sensitivity -> " << (hl_low+hl_up)/2 << std::endl;
        outfile << BI << '\t' << (hl_low+hl_up)/2 << '\n';
}
#else
        if (verbose) std::cout << "Found sensitivity -> " << (hl_low+hl_up)/2 << std::endl;
        outfile << BI << '\t' << (hl_low+hl_up)/2 << '\n';
#endif
    }
    return 0;
}

double GetBayesFactor(double BI, double hl, Json::Value J) {
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
/*    std::cout << "B+S model Integral = "
        << mBS->GetIntegral() << " +- "
        << 100*mBS->GetError()/mBS->GetIntegral() << "%\n";
    std::cout << "B   model Integral = "
        << mB ->GetIntegral() << " +- "
        << 100*mB ->GetError()/mB ->GetIntegral() << "%\n";
    std::cout << "Bayes Factor = " + std::to_string(mgr.BayesFactor(0, 1)) << std::endl;*/

    return mgr.BayesFactor(0,1);
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
