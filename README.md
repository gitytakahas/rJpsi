Analysis flow: https://docs.google.com/presentation/d/1oV9VAjN-gjE4RI6PWifJsd0fmfeRklfViCrAe6K3imk/edit#slide=id.gca82e4db1e_0_0

To use this package, I assume Ntuples already exist (Yuta takes care this part)

https://github.com/UZHCMS/EXOVVNtuplizerRunII/tree/BcMu_10210_NoLepProducers

You can check ntuplizer outputs at 
```
uberftp -ls gsiftp://storage01.lcg.cscs.ch//pnfs/lcg.cscs.ch/cms/trivcat/store/user/ytakahas/
```

BEAWARE:
In some places, you have to change by hand (e.g. output directory), which is currently set to my environment!!


### 0. Setup package
```
cmsrel CMSSW_10_2_10
cd CMSSW_10_2_10/src
cmsenv
git clone https://github.com/gitytakahas/rJpsi.git
```


### 1. Make a flat ntuples 

```
export RJPSI=$PWD
cd $RJPSI/flatter
```

Inspect runTauDisplay_BcJpsiTauNu and see what is doing here. 
When you want to run on a single file, do, 


python runTauDisplay_BcJpsiTauNu.py <-o outputfilename> <-p priority><-t type><-y year><-f inputfile>

-p option: you can choose either "pt" or "multiple". "pt" option will pick up the highest in pT that satisfies vertex prob > 10% and flight significance > 3 sigma (this can be changed). "multiple" option will pick up all the triplets per event

-t option: can be "data", "bg", "signal"

-y option: year. will be used for the pileup reweighting

-f option: input files. you can specify multiple files by comma separated way file1.root,file2.root,... (no space between "," and the next file)


If you are confident, you can submit jobs through T3 batch. How to do it? You can see, for instance, 

```
sh do.sh 
```

In this script, the output of the job submission will be directed to "/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/" but you can choose wherever you want. 


### 2. Update ntuples 

```
cd $RJPSI/updateTuple;
```

Here, we do, 

(1) Hadd flat-ntuples produced in the previous step ---> This take some time, especially for data !! be patient !!
(2) split the flat-ntuples into 20% (for training analysis-BDT) and 80% (for the analysis).

```
sh do.sh
```


### 3. Building MVA based on the xgboost 
```
cd $RJPSI/mva;
sh do.sh 
```

Once finish, your MVA is ready. Inspect Plots_*/* to see the performance.


### 4. add analysis MVA, build in the previous step, and produce the final ROOT files

```
cd $RJPSI/updateTuple; 
sh application.sh 
```

This takes some time, especially for the data. Be Patient!!


In the application step,, the timestamp of the input files will be printed out. Keep an eye on it if it is reasonable.
You should also make sure the "features" as defined in the application.py should match to that of the $RJPSI/mva/TrainModel_XGB.py

this will produce the final tree at the directory, $RJPSI/updateTuple/final_root_* for instance. 


### 5. plotting 

```
cd $RJPSI/anal/dev
```

If you don't have it yet, copy, 
```
cp /t3home/ytakahas/.rootrc ~/.rootrc
cp /t3home/ytakahas/rootlogon.C ~/
mkdir ~/tool
cp -r /t3home/ytakahas/tool ~/tool 
```

Then inspect draw.py what is going to be plotted, and then do, 

```
python draw.py 
```

the list of variables to be plotted is defined in varConfig.py.

To create the datacard and systematics comparison plots, do, 
```
python createFinalDatacard.py
```
