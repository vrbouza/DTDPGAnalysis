# DTDPGAnalysis
Code for CMS DT offline analysis
This package contains the code needed for CMS DT Prompt Offline Analysis and for DT root-ple production.

To install it and run DTNtuple production:

```
cmsrel CMSSW_9_4_0
cd CMSSW_9_4_0/src/

cmsenv

git clone https://github.com/cmsdtoffline/DTDPGAnalysis UserCode/DTDPGAnalysis

scramv1 b -j 5

cd UserCode/DTDPGAnalysis/test

cmsRun RunTree_collisions_cfg.py 

```

