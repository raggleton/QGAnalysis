#pragma once

#include "UHH2/QGAnalysis/include/QGAddModules.h"
#include "UHH2/core/include/AnalysisModule.h"

#include "TUnfoldBinning.h"
#include "TH1.h"
#include "TH2.h"

/**
 * Helper class to handle filling all the hists of one lambda variable configuration,
 * both gen & reco
 *
 * TODO: store/setup the hists (all/split/jackknife/PDF etc)
 * however this requires the book<> method, to put it in the Context obj...not so easy
 */
class LambdaHistsFiller {
public:
    LambdaHistsFiller(const LambdaArgs & pfLambdaArgs,
                      TUnfoldBinning * recoBinning,
                      const LambdaArgs & genLambdaArgs,
                      TUnfoldBinning * genBinning,
                      const std::string & lambdaName);
    // setup for particular reco, gen jets - recalculate variable and TUnfold bin numbers
    void setPassReco(bool passReco);
    void setupReco(float recoJetPt,
                   const LambdaCalculator<PFParticle> & recoJetCalc,
                   bool passReco);
    void setPassGen(bool passGen);
    void setupGen(float genJetPt,
                  const LambdaCalculator<GenParticle> & genJetCalc,
                  bool passGen);
    // setup hists
    void setupRecoHists(uhh2::Context & ctx,
                        const std::string & dirname,
                        bool doSplitHists,
                        bool doFakeHists,
                        int nJackknifeVariations=0,
                        int nPDFVariations=0);

    void setupGenHists(uhh2::Context & ctx,
                       const std::string & dirname,
                       bool doSplitHists,
                       int nJackknifeVariations=0,
                       int nPDFVariations=0);

    void setupResponseHists(uhh2::Context & ctx,
                            const std::string & dirname,
                            bool doSplitHists,
                            int nJackknifeVariations=0,
                            int nPDFVariations=0);

    // Generic hist filling methods
    // These automatically account for e.g. check on # constituents, vai passReco()/passGen()
    void fillRecoTH1(TH1 * h, double weight);
    void fillRecoTH1GenBinning(TH1 * h, double weight);
    void fillGenTH1(TH1 * h, double weight);
    void fillResponseTH2(TH2 * h, double recoWeight, double genWeight);

    // getters
    float recoJetPt() { return recoJetPt_; }
    bool passReco() { return passReco_; }

    float genJetPt() { return genJetPt_; }
    bool passGen() { return passGen_; }

    TH1D * recoHist() { return recoHist_; }
    TH1D * recoSplitHist() { return recoSplitHist_; }
    TH1D * recoHistGenBinning() { return recoHistGenBinning_; }
    TH1D * recoSplitHistGenBinning() { return recoSplitHistGenBinning_; }
    TH1D * recoHistFakes() { return recoHistFakes_; }
    TH1D * recoSplitHistFakes() { return recoSplitHistFakes_; }
    TH1D * recoHistFakesGenBinning() { return recoHistFakesGenBinning_; }
    std::vector<TH1D*> & recoJackknifeVariations() { return recoJackknifeVariations_; }
    std::vector<TH1D*> & recoPDFVariations() { return recoPDFVariations_; }

    TH1D *genHist() { return genHist_; }
    TH1D *genSplitHist() { return genSplitHist_; }
    std::vector<TH1D*> & genJackknifeVariations() { return genJackknifeVariations_; }
    std::vector<TH1D*> & genPDFVariations() { return genPDFVariations_; }

    TH2D * responseHist() { return responseHist_; }
    TH2D * responseSplitHist() { return responseSplitHist_; }
    std::vector<TH2D*> & responseJackknifeVariations() { return responseJackknifeVariations_; }
    std::vector<TH2D*> & responsePDFVariations() { return responsePDFVariations_; }

private:
    /** \brief Create a histogram in the output file - copied from Hists base class
     *
     * HTYPE should be a root histogram type such as TH1D, TH2D, TH1F, etc. Then, use the method with the
     * same constructor arguments you would use for the constructor of the histogram, e.g.:
     * \code
     * // TH1D * hist1 = new TH1D("hist1_name", "hist1_title", 100, 0, 1)
     * // becomes:
     * TH1D * hist1 = book<TH1D>("hist1_name", "hist1_title", 100, 0, 1);
     *
     * // new TH2D("hist2d_name", "hist2d_title", 100, 0, 1, 200, 0, 2)
     * // becomes:
     * book<TH2D>("hist2d_name", "hist2d_title", 100, 0, 1, 200, 0, 2);
     * \endcode
     *
     * The histogram will be saved in the output file under the directory name given in the constructor. To keep
     * the convention that all histograms of one class are in the same directory in the output root file,
     * it is not allowed to specify further subdirectories at this point, and a runtime_error is thrown
     * if trying to add a histogram with a '/' in the name. (If you want to add histograms to subdirectories, use the
     * uhh2::Context::put method directly or create another Hist class with the subdirecory name as constructor argument.)
     *
     * Histograms created via this methods are handed over to the framework via uhh2::Context::put, so the
     * lifetime and memory is managed by the framework. So please do not call delete on the returned pointer.
     */
    template<typename HTYPE, typename... TARGS>
    HTYPE* book(uhh2::Context & ctx, const std::string & dirname, const char * name, TARGS... args){
        static_assert(std::is_base_of<TH1, HTYPE>::value, "for book<HTYPE>, HTYPE, must inherit from TH1!");
        std::string sname = name;
        if(sname.find('/')!=std::string::npos){
            throw std::runtime_error(" name '" + sname + "' illegal: '/' in histogram names not allowed in Hists (use uhh2::Context::put directly for putting histograms in subdirectories)");
        }
        HTYPE * h = new HTYPE(name, std::forward<TARGS>(args)...);
        h->SetName(sname.c_str());
        h->Sumw2();
        ctx.put(dirname, h);
        return h;
    }

    // utility function allowing to use strings as name and title for histogram construction.
    template<typename HTYPE, typename... TARGS>
    HTYPE* book(uhh2::Context & ctx,const std::string & dirname, const std::string & name, const std::string & title, TARGS... args){
        return book<HTYPE>(name.c_str(), title.c_str(), std::forward<TARGS>(args)...);
    }

    // methods to copy the binning from h, and create & book a new TH*D object
    // overload for different string types
    TH1D * copy_book_th1d(uhh2::Context & ctx, const std::string & dirname, TH1 * h, const std::string & newName);
    TH1D * copy_book_th1d(uhh2::Context & ctx, const std::string & dirname, TH1 * h, const TString & newName);
    TH2D * copy_book_th2d(uhh2::Context & ctx, const std::string & dirname, TH2 * h, const std::string & newName);
    TH2D * copy_book_th2d(uhh2::Context & ctx, const std::string & dirname, TH2 * h, const TString & newName);

    TString strWithLambdaName(const std::string & str);

private:
    LambdaArgs pfLambdaArgs_;
    TUnfoldBinning *recoBinning_;
    bool passReco_;
    LambdaArgs genLambdaArgs_;
    TUnfoldBinning *genBinning_;
    bool passGen_;
    std::string lambdaName_;
    int recoBin_, recoBinGenBinning_, genBin_;
    float recoJetPt_, genJetPt_;
    TH1D *recoHist_, *recoSplitHist_,
         *recoHistGenBinning_, *recoSplitHistGenBinning_,
         *recoHistFakes_, *recoSplitHistFakes_, *recoHistFakesGenBinning_;
    std::vector<TH1D*> recoJackknifeVariations_, recoPDFVariations_;
    TH1D *genHist_, *genSplitHist_;
    std::vector<TH1D*> genJackknifeVariations_, genPDFVariations_;
    TH2D *responseHist_, *responseSplitHist_;
    std::vector<TH2D*> responseJackknifeVariations_, responsePDFVariations_;
};