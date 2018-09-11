import ROOT
rng = ROOT.TRandom3(1)
def exp_rng():
    return rng.Exp(1)
def uni_rng():
    return rng.Rndm()
def cauchy_rng():
    return (1-2*uni_rng())/(1-2*uni_rng())


h1 = ROOT.TH1D("h1","uni", # name and title
                             200, -0.5, 1.5 )
h2 = ROOT.TH1D("h2","exp", # name and title
                             200, 0, 6 )
h3 = ROOT.TH1D("h3","cauchy", # name and title
                             200, -3, 3 )
h_unisum = ROOT.TH1D("h3","(x1+x2)", # name and title
                             200, -3, 3 )
n_entries = 10**5
for i in range(n_entries):
    h2.Fill(exp_rng())
    h3.Fill(cauchy_rng())
    x1 = uni_rng()
    h1.Fill(x1)
    h_unisum.Fill((x1+uni_rng())/2)


canv = ROOT.TCanvas("canv","plots for SDA course - ROOT 2", 2000, 2000)
canv.Divide(2,2)
canv.cd(1)
h1.Scale(1.0/n_entries);
h1.Draw("HIST")
canv.cd(2)
h2.Scale(1.0/n_entries);
h2.Draw("HIST")
canv.cd(3)
h3.Scale(1.0/n_entries);
h3.Draw("HIST")
canv.Update()
canv.SaveAs("distributions.png")

canv.Clear()
canv.Divide(1,2)
canv.cd(1)
h1.Draw("HIST")
canv.cd(2)
h_unisum.Scale(1.0/n_entries)
h_unisum.Draw("HIST")
canv.Update()
canv.SaveAs("a.png")
