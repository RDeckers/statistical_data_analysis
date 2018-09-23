#Roel Deckers, Utrecht University, 3762971
import sys
import array
import ROOT
import math

#define our basic generator seeded for constiency
rng_gen = ROOT.TRandom3(314159)
#define our three generator functions
def exp_rng():
    return rng_gen.Exp(1)
def uni_rng():
    return rng_gen.Rndm()
def cauchy_rng(): #not a real cauchy
    return (1-2*uni_rng())/(1-2*uni_rng())

#predefine a dictionary for the data
data = {}
#predefine a dictionary for the histograms
histograms = {}
#how many entries we generate for each generator
n_entries = 10**5
#an array of all our rng functions so we can loop over it.
rngs = [uni_rng, exp_rng, cauchy_rng]
types = ["uniform", "exponential", "cauchy"]
#We use a lot more N than needed so that we can make a nice line plot later
# does not increase running time significantly because we only compute max(N) samples per
# type of histogram as opposed to sum(N) that the naive solution would have.
N = [1,2,5,10,25,50,75,100]
#we cache our results in a root file to save time
#use my student number for the filename
file_name = "3762971.root"

# file stores raw data, so that it doesn't need to be recreated but histograms etc can be freely made.
# each type of generator is a top level folder.
# for each type it stores one TArray for each element of N (named after the element)
# top level also stores one array named N, which unsurprisingly contains N.

#computes the number of bins according to the Freedman-Diaconis rule
def compute_n_bins(x_min, x_max, data):
    #first we need to fetch the 0.25 and 0.75 quantile for the IRQ
    quantiles = array.array('d', [0,0]);
    ROOT.TMath.Quantiles(data.GetSize(), 2, data.GetArray(), quantiles, array.array('d',[0.25, 0.75]), False)
    return int((x_max-x_min)/2*pow(float(data.GetSize()), 1.0/3.0)/(quantiles[1]-quantiles[0]))

def create_data(rng, entries, N):
    #allocate the arrays
    data = [ROOT.TArrayD(entries) for _ in N]
    for entry in range(entries):
        #here we do a little trickery so that we don't
        # take sum(N) samples per entry but max(N)
        # by continouing our sum for N[i] from N[i-1]
        if (entry % (5*10**3)) == 0:
            print entry
        sum = 0
        old_n = 0
        for new_n,array in zip(N,data):
            for n in range(new_n-old_n):
                sum += rng()
            old_n = new_n
            array[entry] = sum/new_n
    return data


# Create and fill a file `filename` with all our data
def create_and_fill_datafile():
    file = ROOT.TFile(file_name,"recreate");
    rootN = ROOT.TArrayI(len(N))
    for i,n in enumerate(N):
        rootN[i] = n
    file.WriteObject(rootN, "N")
    for rng,type in zip(rngs,types):
        file.mkdir(type)
        file.cd(type)
        print("creating data for %s") % type
        data[type] = create_data(rng, n_entries, N)
        print("created data for %s") % type
        for n,array in zip(N,data[type]):
            ROOT.gDirectory.WriteObject(array, str(n))
        file.cd("..")
    file.Write()
    return file


# Load our root file and set histograms using that
# n.b. only works if the file being loaded is a contains at least
# all the N and types that we specified.
def restore_from_file():
    file = ROOT.TFile(file_name)
    N = file.Get("N")
    for type in types:
        dir = file.Get(type)
        data[type] = [dir.Get(str(n)) for n in N]
    n_entries = data["uniform"][0].GetSize()
    print("Restored file with %d entries and N = %s") % (n_entries, [n for n in N])
    return file

#creates histograms and fills them with the data
def create_histograms():
    for type in types:
        #what's the x-range on our histograms per type
        x_range = {"uniform": [-0.25, 1.25], "exponential" : [0.0,6.0], "cauchy" : [-6.0, 6.0]}
        histograms[type] = [
        ROOT.TH1D(
            type+str(n),
            type+str(n),
            #autocompute the number of bins
            compute_n_bins(x_range[type][0], x_range[type][1], array),
            x_range[type][0],
            x_range[type][1]
        ) for n,array in zip(N,data[type])]
        #now fill them with data
        for hist,array in zip(histograms[type],data[type]):
            for x in array:
                hist.Fill(x)
            #normalize the histogram
            hist.Scale(1/(hist.GetBinWidth(1)*n_entries))


def mean_and_stddev(array):
    sum = 0.0
    sum_sq = 0.0
    n = len(array)
    for x in array:
        sum += x
        sum_sq += x*x
    sum /= len(array)
    sum_sq /= len(array)
    return (sum, math.sqrt(sum_sq-sum*sum))

#restore from file if --cached is parsed.
#we keep the file that is returned so that ROOT does not destroy it
if("--cached" in sys.argv):
    file = restore_from_file()
else:
    file = create_and_fill_datafile()

#print the mean and variance for each type
for type in types:
    print(type, mean_and_stddev(data[type][0]))

create_histograms()
canv = ROOT.TCanvas("canv","plots for SDA course - ROOT 3", 2000, 2000)

#plot the stdevs
canv.Clear()
canv.Divide(2,2)
#workaround for roots ownership model.
#if we don't keep this outside of the loop different iteration will overwrite each others result.
graphs = [None for _ in types]
for i,type in enumerate(types):
    canv.cd(i+1)
    D = mean_and_stddev(data[type][0])[1]
    y = array.array('d', [(h.GetStdDev()*math.sqrt(float(n)))/D for h,n in zip(histograms[type],N)])
    x = array.array('d', N)
    graphs[i] = ROOT.TGraphErrors(len(N),x,y)
    graphs[i].SetMarkerColor(ROOT.kRed)
    graphs[i].SetMarkerStyle(21)
    graphs[i].SetMinimum(0)
    graphs[i].SetMaximum(1.2)
    graphs[i].SetTitle(type + ";N;#sqrt{(N)}#sigma_{D,N}/#sigma_{D}")
    graphs[i].Draw("ALP")
    canv.Update()
canv.SaveAs("stddev.png")

#plot the
canv.Clear()
canv.Divide(2,2)
graphs = [None for _ in types]
for i,type in enumerate(types):
    canv.cd(i+1)
    D = mean_and_stddev(data[type][0])[0]
    y = array.array('d', [h.GetMean() for h,n in zip(histograms[type],N)])
    x = array.array('d', N)
    graphs[i] = ROOT.TGraphErrors(len(N),x,y)
    graphs[i].SetMarkerColor(ROOT.kRed)
    graphs[i].SetMarkerStyle(21)
    graphs[i].SetMinimum(-0.1)
    graphs[i].SetMaximum(1.5)
    line = ROOT.TLine(0,D,N[len(N)-1],D);
    graphs[i].SetTitle(type + ";N;mean(D_{N})")
    graphs[i].Draw("ALP")
    canv.Update()
canv.SaveAs("mean.png")


canv.Clear()
canv.Divide(2,2)
#make a plot of our underlying distributions, i.e. with n=1.
for i,type in enumerate(types):
    canv.cd(i+1)
    histograms[type][0].SetLineWidth(3)
    histograms[type][0].SetFillColor(6)
    histograms[type][0].SetAxisRange([-0.5,0,-5][i], [1.5,5,5][i],"X");
    histograms[type][0].SetTitle(type + ";x;P(x)")
    histograms[type][0].Draw("HIST")
canv.Update()
canv.SaveAs("distributions.png")

#then we make a plot for each distribution showing the convergence to a normal
#distribution for increasing N.
ROOT.gStyle.SetOptStat(0)# Do not show statistics box
canv.Clear()
canv.Divide(2,2)
legends = [None for _ in types]
for i,type in enumerate(types):
    canv.cd(i+1)
    legends[i] = ROOT.TLegend(0.8,0.7,0.9,0.9);
    legends[i].SetHeader("N = ...","C");# // option "C" allows to center the header
    print([[n for n in N].index(k) for k in [1,2,10,100]])
    for j in [[n for n in N].index(k) for k in [1,2,10,100]]:
        histograms[type][j].SetLineWidth(3)
        histograms[type][j].SetLineColor(j+1)
        histograms[type][j].SetFillColor(0)
        histograms[type][j].SetAxisRange([-0.25,0,-5][i], [1.25,5,5][i],"X");
        histograms[type][j].SetAxisRange(0, [15,4.5,0.5][i],"y");
        histograms[type][j].SetTitle(type + ";x;P(x)")
        histograms[type][j].Draw("HIST" if j == 0 else "HIST;SAME")
        legends[i].AddEntry(type+str(N[j]), str(N[j]), "l")
    legends[i].Draw()
canv.Update()
canv.SaveAs("dist2.png")
