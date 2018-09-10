# author : Aart Heijboer
# student: Roel Deckers

import ROOT             # This will make root available in your script
import os               # To create output folders when needed
from math import *      # Use math functions (cos, sin, etc)

# -----------------------------
#  make root look a bit better
# -----------------------------
ROOT.gStyle.SetOptStat(0)           # Do not show statistics box
ROOT.gStyle.SetPadLeftMargin(0.14)  # make room for y-title, adjust with pad.SetLeftMargin()
ROOT.gStyle.SetPadRightMargin(0.2)  # make room for y-title, adjust with pad.SetLeftMargin()
ROOT.gStyle.SetTitleOffset(1.8,"y") # adjust with histogram.GetYaxis().SetTitleOffset)
ROOT.gStyle.SetTitleOffset(1.5,"z") # adjust with histogram.GetYaxis().SetTitleOffset)

#used for computing the sigma in the ratios
def binomial_sigma(p, n):
    return sqrt(n*p*(1-p))/n

def to_four_vector(pt, theta, phi, m):
	p = pt/sin(theta)
	px = pt*sin(phi)
	py = pt*cos(phi)
	pz = p*cos(theta)
	E = sqrt(m*m+p*p)
	return (E,px,py,pz)

def invariant_mass_from_4vec(four_vec_1, four_vec_2):
	return sqrt((four_vec_1[0]+four_vec_2[0])**2-
 	 (four_vec_1[1]+four_vec_2[1])**2+
	 (four_vec_1[2]+four_vec_2[2])**2+
	 (four_vec_1[3]+four_vec_2[3])**2)

def invariant_mass(pt1, theta1, phi1, pt2, theta2, phi2):
	return invariant_mass_from_4vec(
      #rest mass of muon hard coded
	  to_four_vector(pt1, theta1, phi1, 0.105),
	  to_four_vector(pt2, theta2, phi2, 0.105)
	)

# We start by making a ROOT canvas, which is a placeholder for plots
# and other graphics. This is on object of type TCanvas
# (https://root.cern.ch/root/html/TCanvas.html).

plot_width = 900 #slightly bigger than the height because of larger offset
plot_height = 800
canv = ROOT.TCanvas("canv","plots for SDA course - ROOT 1", 2*plot_width,2*plot_height)
canv.Divide(2,2)  # See link above for documentation of the Divide method.

#no titles because we'll be embeding in a document
histogram_of_pt1_plus_pt2 = ROOT.TH1D("histogram_of_pt1_plus_pt2","", # name and title
                             500, 0, 6000 )
histogram_of_invariant_mass = ROOT.TH1D("histogram_of_invariant_mass","", # name and title
                             500, 0, 6000 )
histogram_of_pt1_v_pt2 = ROOT.TH2D("histogram_of_pt1_v_pt2","", # name and title
                             500, 0, 4000, 500, 0, 4000 )

# We are ready to open the datafile and loop over the lines/event
# (remember each line !that is not a comment! is an event).
input_file = open("/home/roel/repos/statistical_data_analysis/datasets/dimuon-dataset.txt")       # open the file

line_counter = 0
pt1_gt_1tev = 0
pt2_gt_1tev = 0
pt1_or_pt2_gt_1tev = 0
pt1_and_pt2_gt_1tev = 0
for line in input_file :   # loop over every line in the input file
    # skip any comment lines which start with '#'
    if line.startswith("#") : continue
    line_counter += 1

    pt1, theta1, phi1, pt2, theta2, phi2 = [ float(x) for x in line.split() ]
    histogram_of_pt1_plus_pt2.Fill( pt1 + pt2 )
    histogram_of_invariant_mass.Fill(invariant_mass(pt1, theta1, phi1, pt2, theta2, phi2))
    histogram_of_pt1_v_pt2.Fill(pt1, pt2)
    if pt1 > 1000:
	pt1_gt_1tev += 1
    if pt2 > 1000:
	pt2_gt_1tev += 1
    #pt2 OR pt1 > 1000
    if max(pt2,pt1) > 1000:
	pt1_or_pt2_gt_1tev += 1
    #pt1 AND pt2 > 1000
    if min(pt1, pt2) > 1000:
	pt1_and_pt2_gt_1tev += 1

    if line_counter % 10000 ==0 : print "processed", line_counter," events."

histogram_of_pt1_plus_pt2.SetLineWidth(1)   # nice thin line so we can clearly see any bumps
histogram_of_pt1_plus_pt2.SetFillColor(5)   # yellow fill color
histogram_of_pt1_plus_pt2.SetXTitle("transverse momentum of #mu_{1} + #mu_{2} (GeV)")
histogram_of_pt1_plus_pt2.SetYTitle("events per bin")


histogram_of_pt1_v_pt2.SetXTitle("transverse momentum of #mu_{1} (GeV)")
histogram_of_pt1_v_pt2.SetYTitle("transverse momentum of #mu_{2} (GeV)")
histogram_of_pt1_v_pt2.SetZTitle("events per bin")

histogram_of_invariant_mass.SetLineWidth(1)
histogram_of_invariant_mass.SetFillColor(5)
histogram_of_invariant_mass.SetXTitle("Invariant mass of #mu_{1} and #mu_{2} (GeV)")
histogram_of_invariant_mass.SetYTitle("events per bin")

canv.cd(1)
histogram_of_pt1_plus_pt2.Draw()
canv.cd(2)
#use colors for drawing. Locally we have ROOT6 (7 actually) installed so this actually
#produces a non-abominable color plot.
histogram_of_pt1_v_pt2.Draw("colz")
canv.cd(3)
histogram_of_invariant_mass.Draw()
# Actually write everything to the canvas
canv.Update()

#save our images in seperate pngs
if not os.path.exists("graphs"):
    os.makedirs("graphs")
img= ROOT.TImage.Create()
#trim a little bit of the graph because it has excessive padding
img.FromPad(canv, 0, 0, int(plot_width*0.9), plot_height);
img.WriteImage("graphs/pt_sum.png")
img.FromPad(canv, plot_width,0,plot_width,plot_height);
img.WriteImage("graphs/2d_hist.png")
img.FromPad(canv, 0,plot_height,int(plot_width*0.9),plot_height);
img.WriteImage("graphs/invariant_mass.png")


def prob(x):
	return float(x)/line_counter
def sigma(x):
    return binomial_sigma(prob(x), line_counter)

print "pt1 > 1TeV:\t\t", prob(pt1_gt_1tev), sigma(pt1_gt_1tev)
print "pt2 > 1TeV:\t\t", prob(pt2_gt_1tev), sigma(pt2_gt_1tev)
print "pt1|pt2 > 1 TeV:\t", prob(pt1_or_pt2_gt_1tev), sigma(pt1_or_pt2_gt_1tev)
print "pt1&pt2 > 1 TeV:\t", prob(pt1_and_pt2_gt_1tev), sigma(pt1_and_pt2_gt_1tev)
print "pt1&pt2 > 1 TeV:\t", prob(pt1_gt_1tev)+prob(pt2_gt_1tev)-prob(pt1_or_pt2_gt_1tev), binomial_sigma(prob(pt1_gt_1tev)+prob(pt2_gt_1tev)-prob(pt1_or_pt2_gt_1tev), line_counter)
