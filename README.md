# GEMCSCBendingAnalyzer

## how to check out cmssw and this package
cmsrel CMSSW_10_1_5  

cd CMSSW_10_1_5/src/

cmsenv

git cms-init

git clone https://github.com/tahuang1991/GEMCSCBendingAnalyzer.git

scram b -j 9


## SliceTestAnalysis.cc: package for GEM residual 
### run on data
cmsRun runSliceTestAnalysis.py

### run on MC
cmsRun runSliceTestAnalysis_MC.py

The packaged is inherited from Jason Lee's MuonPerformance and it is used for GEM related analysis at TAMU group. The major target is to analyse
GEM-CSC bending angle in real data from CMS. 
