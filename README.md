# Quark-Gluon analysis

Does analysis of quark/gluon jets in data & MC. Uses UHH2 & SFrame frameworks.
Based on 80X, using Moriond samples

## Installation

**Need my fork of UHH2**: https://github.com/raggleton/UHH2/tree/RunII_80X_v3_QGJetStudies

- First install SFrame & UHH2 via the usual instructions on their wiki
- Clone this repo inside `CMSSW_X_Y_Z/src/UHH2`
- Compile with `make -j10`

## Make NTuples from MiniAOD

This is required to add in extra info to the Jets and GenJets.

TODO

## Running analysis code

Either locally:

```
sframe_main config/QGAnalaysis.xml
```

Or on batch system using SFrameBatch:

```
cd config
sframe_batch -s QGAnalysis.xml
```

(use my fork as fixes project name, job name)


