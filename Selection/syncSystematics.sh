#!/bin/bash -e
#
# Sync systematics files from indiv folders to one big one

set -u

VARIATIONS=("Up" "Down")
SYSTS=("chargedHadronShift" "neutralHadronShift" "photonShift" "jecsmear_direction" "jersmear_direction" "pileup_direction" "track_direction")
# SYSTS=("pileup_direction")
# SYSTS=("jecsmear_direction" "jersmear_direction")
# SYSTS=("track_direction")

OUTPUTDIR="workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet/"
OUTPUTDIR="workdir_ak4puppi_data_jetAsymCut_pt1RecoConstituents_V11JEC_JER_target0p6_wta_groomed_fwdcenDijet_ptAveBinning/"
OUTPUTDIR="workdir_ak4puppi_data_jetAsymCut_pt1RecoConstituents_V11JEC_JER_target0p5_wta_groomed_fwdcenDijet_ptAveBinning/"

OUTPUTDIR="workdir_ak4puppi_data_jetAsymCut_pt1RecoConstituents_V11JEC_JER_target0p5_wta_groomed_fwdcenDijet_ZReweight/"
OUTPUTDIR="workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet_Zreweight_noZjet2Cut_zPt30_zJetAyms0p4/"

NAFSTEM="/nfs/dust/cms/user/aggleton/QG/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/Selection/"
STEMSRCDIR="${NAFSTEM}/MGPythia/workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p6_noZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_ptAveBinning"
STEMSRCDIR="${NAFSTEM}/MGPythia/workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_noZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_ptAveBinning"
STEMSRCDIR="${NAFSTEM}/MGPythia/workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut_zPt30"
STEMSRCDIR="${NAFSTEM}/MGPythia/workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut_zPt30_zJetAsym0p4"

STEMSRCDIR="${NAFSTEM}/MGPythia/workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut_zPt30_newBinning"
OUTPUTDIR="workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet_Zreweight_noZjet2Cut_newBinning"

STEMSRCDIR="${NAFSTEM}/MGPythia/workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut_zPt30_newBinningFixed_trkSFMyFormula"
OUTPUTDIR="workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet_Zreweight_noZjet2Cut_newBinningFixed_trkSFMYFormula"

STEMSRCDIR="${NAFSTEM}/MGPythia/workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut_zPt30_trkSF_wtaAK_gt1Constit_passGenRsp"
OUTPUTDIR="workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_ZReweight_wta_groomed_fwdcenDijet_Zreweight_noZjet2Cut_zPt30_trkSF_wtaAK_gt1Constit_passGenRsp"

STEMSRCDIR="${NAFSTEM}/MGPythia/workdir_ak8puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut_zPt30_trkSF_wtaAK_gt1Constit_passGenRsp"
OUTPUTDIR="workdir_ak8puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_ZReweight_wta_groomed_fwdcenDijet_Zreweight_noZjet2Cut_zPt30_trkSF_wtaAK_gt1Constit_passGenRsp"

STEMSRCDIR="${NAFSTEM}/MGPythia/workdir_ak4puppi_mgpythia_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_trkSF_wtaAK_passGenRsp_sameGenCuts"
OUTPUTDIR="workdir_ak4puppi_data_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_trkSF_wtaAK_fixPrescales_sameGenCuts"

STEMSRCDIR="${NAFSTEM}/MGPythia/workdir_ak8puppi_mgpythia_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_trkSF_wtaAK_passGenRsp_sameGenCuts"
OUTPUTDIR="workdir_ak8puppi_data_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_trkSF_wtaAK_fixPrescales_sameGenCuts"

STEMSRCDIR="MGPythia/workdir_ak4puppi_mgpythia_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_trkSF_wtaAK_passGenRsp_sameGenCuts_fixPassGen"
OUTPUTDIR="workdir_ak4puppi_data_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_trkSF_wtaAK_fixPrescales_sameGenCuts_fixPassGen"

STEMSRCDIR="MGPythia/workdir_ak8puppi_mgpythia_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_trkSF_wtaAK_passGenRsp_sameGenCuts_fixPassGen"
OUTPUTDIR="workdir_ak8puppi_data_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_trkSF_wtaAK_fixPrescales_sameGenCuts_fixPassGen"

STEMSRCDIR="MGPythia/workdir_ak4puppi_mgpythia_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_trkSF_wtaAK_passGenRsp_sameGenCuts_fixPassGen_jackknife25"
OUTPUTDIR="workdir_ak4puppi_data_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_trkSF_wtaAK_fixPrescales_sameGenCuts_fixPassGen_jackknife25"

STEMSRCDIR="MGPythia/workdir_ak8puppi_mgpythia_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_trkSF_wtaAK_passGenRsp_sameGenCuts_fixPassGen_jackknife25"
OUTPUTDIR="workdir_ak8puppi_data_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_trkSF_wtaAK_fixPrescales_sameGenCuts_fixPassGen_jackknife25"

STEMSRCDIR="MGPythia/workdir_ak4puppi_mgpythia_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_wtaAK_passGenRsp_sameGenCuts_jackknife25_rapidity1p7_trkSFonlyLambda"
OUTPUTDIR="workdir_ak4puppi_data_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_wtaAK_sameGenCuts_jackknife25_rapidity1p7_trkSFonlyLambda"

STEMSRCDIR="MGPythia/workdir_ak8puppi_mgpythia_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_wtaAK_passGenRsp_sameGenCuts_jackknife25_rapidity1p7_trkSFonlyLambda"
OUTPUTDIR="workdir_ak8puppi_data_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_wtaAK_sameGenCuts_jackknife25_rapidity1p7_trkSFonlyLambda"

STEMSRCDIR="MGPythia/workdir_ak4puppi_mgpythia_target0p5_ZReweight_wta_groomed_fwdcenDijet_noPtHatCut_noZjet2Cut_zPt30_sameGenCuts_jackknife25_rapidity1p7_trkSFonlyLambda_constitPt0_fixGroomingConstitCheck"
OUTPUTDIR="workdir_ak4puppi_data_target0p5_ZReweight_wta_groomed_fwdcenDijet_noPtHatCut_noZjet2Cut_zPt30_sameGenCuts_jackknife25_rapidity1p7_trkSFonlyLambda_constitPt0_fixGroomingConstitCheck"

STEMSRCDIR="MGPythia/workdir_ak4puppi_mgpythia_target0p5_ZReweight_wta_groomed_fwdcenDijet_noPtHatCut_noZjet2Cut_zPt30_sameGenCuts_jackknife25_rapidity1p7_trkSFonlyLambda_constitPt1_fixGroomingConstitCheck"
OUTPUTDIR="workdir_ak4puppi_data_target0p5_ZReweight_wta_groomed_fwdcenDijet_noPtHatCut_noZjet2Cut_zPt30_sameGenCuts_jackknife25_rapidity1p7_trkSFonlyLambda_constitPt1_fixGroomingConstitCheck"

STEMSRCDIR="MGPythia/workdir_102X_v2_ak4puppi_mgpythia_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_fixLambda_newBinning2_zjAsym"
OUTPUTDIR="workdir_102X_v3data_v2mc_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_fixLambda_newBinning2_zjAsym"

# STEMSRCDIR="MGPythia/workdir_102X_v2_ak8puppi_mgpythia_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_fixLambda_newBinning2_zjAsym"
# OUTPUTDIR="workdir_102X_v3data_v2mc_ak8puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_fixLambda_newBinning2_zjAsym"

# STEMSRCDIR="MGPythia/workdir_102X_v2_ak4puppi_mgpythia_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_fixLambda_newBinning2_genjetGhostFlav_noBCpref"
# OUTPUTDIR="workdir_102X_v3data_v2mc_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_fixLambda_newBinning2_genjetGhostFlav_noBCpref"

# STEMSRCDIR="MGPythia/workdir_102X_v2_ak8puppi_mgpythia_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_fixLambda_newBinning2_genjetGhostFlav_noBCpref"
# OUTPUTDIR="workdir_102X_v3data_v2mc_ak8puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_fixLambda_newBinning2_genjetGhostFlav_noBCpref"

# STEMSRCDIR="MGPythia/workdir_102X_v2_ak4puppi_mgpythia_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_newBinning4"
# OUTPUTDIR="workdir_102X_v3data_v2mc_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_newBinning4"

# STEMSRCDIR="MGPythia/workdir_102X_v2_ak8puppi_mgpythia_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_newBinning4"
# OUTPUTDIR="workdir_102X_v3data_v2mc_ak8puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_newBinning4"

STEMSRCDIR="MGPythia/workdir_102X_v2_ak4puppi_mgpythia_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_jetIDSel_newBinning5"
OUTPUTDIR="workdir_102X_v3data_v2mc_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_jetIDSel_newBinning5"

# STEMSRCDIR="MGPythia/workdir_102X_v2_ak8puppi_mgpythia_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_jetIDSel_newBinning5"
# OUTPUTDIR="workdir_102X_v3data_v2mc_ak8puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_jetIDSel_newBinning5"

STEMSRCDIR="MGPythia/workdir_102X_v2_ak4puppi_mgpythia_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_jetIDSel_newBinning5_mergeZpJPt"
OUTPUTDIR="workdir_102X_v3data_v2mc_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_jetIDSel_newBinning5_mergeZpJPt"

STEMSRCDIR="MGPythia/workdir_102X_v2_ak8puppi_mgpythia_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_jetIDSel_newBinning5_mergeZpJPt"
OUTPUTDIR="workdir_102X_v3data_v2mc_ak8puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_jetIDSel_newBinning5_mergeZpJPt"


function copy {
    src="$1"
    dest="$2"
    rsync -avhPt "${src}" "${dest}" || echo -e "Missing ${src} \n"
}

echo "------------------------------------------------------------------"
echo ">>> MG+PYTHIA files"
echo "------------------------------------------------------------------"
ODIR="${OUTPUTDIR}/"
mkdir -p "${ODIR}"
copy "${STEMSRCDIR}/uhh2.AnalysisModuleRunner.MC.MC_QCD.root" "${ODIR}/"
copy "${STEMSRCDIR}/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root" "${ODIR}/"

# exit 0

echo "------------------------------------------------------------------"
echo ">>> HERWIG++"
echo "------------------------------------------------------------------"
base=$(basename $STEMSRCDIR)
base=${base/mgpythia/herwig}
echo $base
copy "Herwig/$base/uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD.root" "${ODIR}/"
copy "Herwig/$base/uhh2.AnalysisModuleRunner.MC.MC_HERWIG_DYJetsToLL_merged_PartonKtMin300.root" "${ODIR}/"

echo "------------------------------------------------------------------"
echo ">>> PYTHIA8"
echo "------------------------------------------------------------------"
base=$(basename $STEMSRCDIR)
base=${base/mgpythia/pythia}
echo $base
copy "Pythia/$base/uhh2.AnalysisModuleRunner.MC.MC_PYTHIA-QCD.root" "${ODIR}/"

for syst in ${SYSTS[@]}; do
    for var in ${VARIATIONS[@]}; do
        echo "------------------------------------------------------------------"
        echo ">>> ${syst}${var}"
        echo "------------------------------------------------------------------"
        ODIR="${OUTPUTDIR}/systematics_files/${syst}${var}"
        mkdir -p "${ODIR}"
        copy "${STEMSRCDIR}_${syst}${var}/uhh2.AnalysisModuleRunner.MC.MC_QCD.root" "${ODIR}/"
        copy "${STEMSRCDIR}_${syst}${var}/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root" "${ODIR}/"
    done
done

# Just list these ones explicitly, cba to write a loop
SYSTS=(
    "ScaleVariationMuRNominal_ScaleVariationMuFUp"
    "ScaleVariationMuRNominal_ScaleVariationMuFDown"
    "ScaleVariationMuRUp_ScaleVariationMuFNominal"
    "ScaleVariationMuRDown_ScaleVariationMuFNominal"
    "ScaleVariationMuRUp_ScaleVariationMuFUp"
    "ScaleVariationMuRDown_ScaleVariationMuFDown"
    "PDFvariationsTrue"
)

for syst in ${SYSTS[@]}; do
    echo "------------------------------------------------------------------"
    echo ">>> ${syst}"
    echo "------------------------------------------------------------------"
    ODIR="${OUTPUTDIR}/systematics_files/${syst}"
    mkdir -p "${ODIR}"
    copy "${STEMSRCDIR}_${syst}/uhh2.AnalysisModuleRunner.MC.MC_QCD.root" "${ODIR}/"
    copy "${STEMSRCDIR}_${syst}/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root" "${ODIR}/"
done

