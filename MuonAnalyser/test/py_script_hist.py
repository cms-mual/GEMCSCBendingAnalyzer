import ROOT
import os
F=ROOT.TFile.Open('histo.root');

def plot_hist(branch_name):
	c1=ROOT.TCanvas("c1","New Graph",0,0,1600,1000);
	h=F.Get("SliceTestAnalysis/MuonData");
	h.Draw(branch_name);
	k=str(branch_name)+".png"
	os.chdir("/afs/cern.ch/user/m/mkhurana/CMSSW_10_1_5/src/GEMCSCBendingAnalyzer/MuonAnalyser/test/Graphs")

	c1.Print(k);
	
	os.chdir("/afs/cern.ch/user/m/mkhurana/CMSSW_10_1_5/src/GEMCSCBendingAnalyzer/MuonAnalyser/test")


branch_list=["lumi","run","event","muonpt","muoneta","muonphi","muoncharge","muonendcap","has_TightID","isGood_GE11","has_GE11","has_ME11","rechit_phi_GE11",
"prop_phi_GE11","rechit_phi_ME11","prop_phi_ME11"]

for name in branch_list:
	plot_hist(name);
	


#print(h.Print());
#for bname in branch_list:
#	print bname
#print(len(branch_list))
#h.MakeClass("my class")
#k=h.GetListOfBranches()
#print(k.Print());
#c1.Divide(4,4);
#c1.cd(2);
#h.Draw("has_ME11");
#c1.cd(3);
#h.Draw("muonphi");
#c1.cd(4);
#h.Draw("muonpt");
#c1.Print("output_hist.png")

