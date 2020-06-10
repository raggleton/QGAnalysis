#!/bin/bash -e
#
# Sync systematics files from indiv folders to one big one

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

# CMD="rsync -ahvzP"
CMD="cp"
CMD="rsync -avhPt"

echo "------------------------------------------------------------------"
echo ">>> Main files"
echo "------------------------------------------------------------------"
ODIR="${OUTPUTDIR}/"
mkdir -p "${ODIR}"
$CMD "${STEMSRCDIR}/uhh2.AnalysisModuleRunner.MC.MC_QCD.root" "${ODIR}/"
$CMD "${STEMSRCDIR}/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root" "${ODIR}/"

echo "------------------------------------------------------------------"
echo ">>> HERWIG++"
echo "------------------------------------------------------------------"
base=$(basename $STEMSRCDIR)
base=${base/mgpythia/herwig}
echo $base
$CMD "Herwig/$base/uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD.root" "${ODIR}/"
$CMD "Herwig/$base/uhh2.AnalysisModuleRunner.MC.MC_HERWIG_DYJetsToLL_merged_PartonKtMin300.root" "${ODIR}/"

echo "------------------------------------------------------------------"
echo ">>> PYTHIA8"
echo "------------------------------------------------------------------"
base=$(basename $STEMSRCDIR)
base=${base/mgpythia/pythia}
echo $base
$CMD "PythiaOnlyFlat/$base/uhh2.AnalysisModuleRunner.MC.MC_PYTHIA-QCD.root" "${ODIR}/"

for syst in ${SYSTS[@]}; do
    for var in ${VARIATIONS[@]}; do
        echo "------------------------------------------------------------------"
        echo ">>> ${syst}${var}"
        echo "------------------------------------------------------------------"
        ODIR="${OUTPUTDIR}/systematics_files/${syst}${var}"
        mkdir -p "${ODIR}"
        $CMD "${STEMSRCDIR}_${syst}${var}/uhh2.AnalysisModuleRunner.MC.MC_QCD.root" "${ODIR}/"
        $CMD "${STEMSRCDIR}_${syst}${var}/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root" "${ODIR}/"
    done
done

# Just list these ones explicitly, cba to write a loop
SYSTS=(
    "ScaleVariationMuRNom_ScaleVariationMuFUp"
    "ScaleVariationMuRNom_ScaleVariationMuFDown"
    "ScaleVariationMuRUp_ScaleVariationMuFNom"
    "ScaleVariationMuRDown_ScaleVariationMuFNom"
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
    $CMD "${STEMSRCDIR}_${syst}/uhh2.AnalysisModuleRunner.MC.MC_QCD.root" "${ODIR}/"
    $CMD "${STEMSRCDIR}_${syst}/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root" "${ODIR}/"
done

