#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include <fstream>

#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace LHAPDF;
using namespace std;

int main(int argc, char* argv[]) {

    // const PDF* pdfCTEQ = mkPDF(10042); // CTEQ6L1 - yes this is the same as CTEQ6ll, I think
    const PDF* pdfCTEQ = mkPDF("cteq6l1"); // CTEQ6L1 - yes this is the same as CTEQ6ll, I think
    // const PDF* pdfNNPDF = mkPDF(263000); // NNPDF30_lo_as_0130
    const PDF* pdfNNPDF = mkPDF("NNPDF30_lo_as_0130"); // NNPDF30_lo_as_0130

    const double MINLOGX = -5;
    const double MAXLOGX = 0;
    const double DX = 0.1;
    vector<double> xValues = {};
    vector<double> exValues = {};
    for (double l=MINLOGX; l<MAXLOGX-DX; l+= DX){
        xValues.push_back(pow(10, l));
        exValues.push_back(0);
    }
    const int NX = xValues.size();

    const double MINLOGQ2 = 2;
    const double MAXLOGQ2 = 6;
    const double DQ2 = 2;
    vector<double> q2Values = {};
    for (double q=MINLOGQ2; q<=MAXLOGQ2; q+= DQ2){
        q2Values.push_back(pow(10, q));
    }

    vector<int> pids = pdfCTEQ->flavors();
    for (int pid : pids) {
        const string spid = lexical_cast<string>(pid);

        for (double q2 : q2Values) {
            // for each flavour, q2, do plot of xf(x) vs x for both PDFs
            vector<double> dataCTEQ;
            vector<double> edataCTEQ;
            vector<double> dataNNPDF;
            vector<double> edataNNPDF;
            for (double x : xValues) {
                dataCTEQ.push_back(pdfCTEQ->xfxQ2(pid, x, q2));
                edataCTEQ.push_back(0);
                dataNNPDF.push_back(pdfNNPDF->xfxQ2(pid, x, q2));
                edataNNPDF.push_back(0);
                // cout << pid << " : " << q2 << " : " << dataCTEQ.back() << " : " << dataNNPDF.back() << endl;
            }

           const string qstr = lexical_cast<string>(q2);
           TGraphErrors * grCTEQ = new TGraphErrors(NX, &xValues[0], &dataCTEQ[0], &exValues[0], &edataCTEQ[0]);
           TGraphErrors * grNNPDF = new TGraphErrors(NX, &xValues[0], &dataNNPDF[0], &exValues[0], &edataNNPDF[0]);
           grCTEQ->SetTitle(("PDGID: " + spid + " Q^{2}= " + qstr + " GeV^{2};x;xf(x)").c_str());
           TCanvas * c = new TCanvas("c", "", 800, 800);
           c->SetTicks(1, 1);
           c->SetLogx();
           gPad->SetLeftMargin(0.16);

           grCTEQ->SetMarkerColor(kRed);
           grCTEQ->SetMarkerStyle(22);
           grCTEQ->SetLineColor(kRed);
           grNNPDF->SetMarkerColor(kBlack);
           grNNPDF->SetMarkerStyle(23);
           grNNPDF->SetLineColor(kBlack);
           grCTEQ->Draw("ALP");

           grNNPDF->Draw("LP SAME");
           TLegend leg(0.6, 0.75, 0.88, 0.88);
           leg.AddEntry(grCTEQ, "CTEQ6L1", "LP");
           leg.AddEntry(grNNPDF, "NNPDF3.0 LO aS=0.130", "LP");
           leg.Draw();
           c->SaveAs(("pdfPlots/cteq_vs_nnpdf_"+spid+"_q2_"+qstr+".pdf").c_str());
           c->SetLogy();
           grCTEQ->SetMinimum(dataCTEQ.back() / 5.);
           grNNPDF->SetMinimum(dataCTEQ.back() / 5.);
           leg.SetX1NDC(0.2);
           leg.SetX2NDC(0.48);
           leg.SetY1NDC(0.15);
           leg.SetY2NDC(0.28);
           leg.Draw();
           c->SaveAs(("pdfPlots/cteq_vs_nnpdf_"+spid+"_q2_"+qstr+"_logy.pdf").c_str());

           delete grCTEQ;
           delete grNNPDF;
           delete c;
        } // end q2 loop
    } // end pid loop

    delete pdfCTEQ;
    delete pdfNNPDF;
    return 0;
}
