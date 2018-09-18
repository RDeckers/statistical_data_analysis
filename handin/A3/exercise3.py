import sys
import ROOT
rng_gen = ROOT.TRandom3(1)
def exp_rng():
    return rng_gen.Exp(1)
def uni_rng():
    return rng_gen.Rndm()
#cauchy is not cauchy, it's
def cauchy_rng():
    return (1-2*uni_rng())/(1-2*uni_rng())

histograms = {}
n_entries = 10**5
rng = [uni_rng, exp_rng, cauchy_rng]
types = ["uniform", "exponential", "cauchy"]
N = [1,2,10,50,100,200,300,400,500]
file_name = "3762971.root"

def fill_histograms():
    for j,type in enumerate(types):
        for entry in range(n_entries):
            if(0 == (entry % 10000)):
                print(type, entry)
            sum = 0
            old_n = 0
            for i in range(len(N)):
                for n in range(N[i]-old_n):
                    sum += rng[j]()
                old_n = N[i]
                histograms[type][i].Fill(sum/N[i])

def normalize_histograms():
    for type in types:
        for hist in histograms[type]:
            hist.Scale(1/(hist.GetBinWidth(1)*n_entries))

def restore_from_file():
    file = ROOT.TFile(file_name)
    for type in types:
        histograms[type] = [file.Get(type + "/" + type+"_"+str(n)) for n in N]
    return file

def create_file():
    file = ROOT.TFile(file_name,"recreate");
    for type in types:
        file.mkdir(type)
        file.cd(type)
        histograms[type] = [ROOT.TH1D(type+"_"+str(n),type+"_"+str(n)+"_title", 300, -6, 6 ) for n in N]
        file.cd("..")
    fill_histograms()
    normalize_histograms()
    file.Write()
    return file



#we keep the file that is returned so that ROOT does not destroy it
if("--cached" in sys.argv):
    file = restore_from_file()
else:
    file = create_file()
print("after create")
canv = ROOT.TCanvas("canv","plots for SDA course - ROOT 2", 2000, 2000)
canv.Divide(2,2)
for i,type in enumerate(types):
    canv.cd(i+1)
    histograms[type][0].Draw("HIST")
canv.Update()
canv.SaveAs("distributions.png")
