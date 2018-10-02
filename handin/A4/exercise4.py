import ROOT
import array
from math import *


#define our basic generator seeded for constiency
rng_gen = ROOT.TRandom3()


def random_u():
    (x1,x2) = (rng_gen.Rndm(), rng_gen.Rndm())
    if (1.0-4.0/3.0*x1*(1.0-x1)) >= x2:
        return x1
    else:
        return random_u()

rho0 = 1.225
a = 8420.0

def rho(h):
    return rho0*exp(-h/a)

def rho_int(h):
    return -a*rho0*exp(-h/a)

def rho_int_inv(y):
    return a*log(-(a*rho0)/y)

xpp = 380
xbrem = 263

def gen_interaction(x, h1):
    Xreal = rng_gen.Exp(x)
    return rho_int_inv(rho_int(h1)-Xreal)


# provide better print functionality for ROOT TVector3
ROOT.TVector3.__str__ = lambda v : "({:g},{:g},{:g})".format( v.X(), v.Y(), v.Z() )


class Particle(object):

    "A class to store the information of a particle in our air shower."

    def __init__( self ) :

        self.direction  = ROOT.TVector3(0,0,0)
        self.energy     = 0
        self.kind       = 0 # 1 = gamma, 2 = electron
        self.start_pos  = ROOT.TVector3(0,0,0)
        self.end_pos    = ROOT.TVector3(0,0,0)
        self.col = 1

    def __str__(self):

        s = " Particle "
        s+= " start: " +str( self.start_pos )
        s+= " end: " +str( self.end_pos )
        s+= " Energy: " +str(self.energy)
        return s

    def __repr__(self):
        return self.__str__()

    def color(self):
        return [9, 46, 8][self.kind]
    def is_electron(self):
        return self.kind == 2
    def is_gamma(self):
        return self.kind == 1
    def momentum(self):
        if self.is_electron():
            return sqrt(self.energy**2-0.510**2)*self.direction
        else:
            return self.direction*self.energy
    def create_end_pos(self):
        #we compute the height difference and then use that as the _vector_ length between
        #start an end. The alternative would be to use this as the difference in Z an compute x,y from there
        # but this is an approximation (particle moves straight down) which should be offset
        # slightly by using dZ as a vector length instead. (+it's easier)
        delta = self.start_pos[2]-gen_interaction([xpp,xpp,xbrem][self.kind], self.start_pos[2])
        self.end_pos = self.start_pos + delta*self.direction
        return self.end_pos
    def split(self):
        x1 = random_u()
        p1 = Particle()
        p2 = Particle()

        p1.energy = self.energy*x1
        p2.energy = self.energy*(1-x1)

        p1.start_pos = p2.start_pos = self.end_pos
        if self.is_gamma():
            p1.kind = p2.kind = 2
        else:
            p1.kind = 1 if rng_gen.Rndm() > 0.5 else 2
            p2.kind = p1.kind ^ 3
        phi = rng_gen.Rndm() * 2 * pi
        p1.direction = direction_at_angle(self.direction, 0.510/p1.energy, phi)
        p2.direction = direction_at_angle(self.direction, 0.510/p2.energy, pi-phi)
        p1.create_end_pos()
        p2.create_end_pos()
        # delta_momentum =  self.momentum() - p1.momentum() - p2.momentum()
        # print(self.is_gamma(), delta_momentum.Mag()/self.momentum().Mag())
        return [p1, p2]
def direction_at_angle( initial_direction, theta, phi ):

    """
    Return a vector that makes an angle theta with
    respect to the TVector3 initial_direction.
    The remaining degree of freedom is the angel phi.
    """
    v = ROOT.TVector3( sin(theta)* cos(phi), sin(theta)*sin(phi), cos(theta ))
    v.RotateY( initial_direction.Theta() )
    v.RotateZ( initial_direction.Phi()  )
    return v


def plot_shower( particles ,
                 title = "Worst title ever.",
                 xysize = 100,
                 zsize = 50000,
                 min_E = 85):

    """
    Plot a list of particles.
    You may use and modify this function to make plots of your showers
    """
    ROOT.gStyle.SetOptStat(0)

    # We use a dummy 3d histogram to define the
    # 3d-axes and coordinate ranges.

    h = ROOT.TH3D(title,title,
                  1, -xysize,xysize,
                  1, -xysize,xysize,
                  1, 0, zsize)

    h.GetXaxis().SetTitleOffset(1.7)
    h.SetXTitle("x (m)")
    h.SetYTitle("y (m)")
    h.SetZTitle("Height from surface (m)")

    h.DrawCopy();

    for p in particles :
        # create a line object to draw
        if p.energy > min_E:
            line = ROOT.TPolyLine3D(2)
            line.SetPoint(0, p.start_pos.X(), p.start_pos.Y(), p.start_pos.Z() )
            line.SetPoint(1, p.end_pos.X(),   p.end_pos.Y(),   p.end_pos.Z()   )
            line.SetLineColor( p.color() )
            line.Draw()

            # The follwing tells python that line should not be cleaned up.
            ROOT.SetOwnership(line, False )
    ROOT.SetOwnership( h, False )




def create_shower(energy):
    particles = [] # start with empty list
    p = Particle()
    p.energy = energy
    p.start_pos = ROOT.TVector3( 0,0, 50000 ) # start at 50 km height
    p.direction = direction_at_angle( ROOT.TVector3(0,0,-1), 0, 0 ) #straight down
    p.create_end_pos()
    p.kind = 1
    particles.append(p) # put our particle in the list
    #function to make children of a generation
    def add_generation(parents):
        particles = []
        for p in parents:
            #if enough energy and above ground:
            if (p.energy >= 85) and (p.end_pos[2] >= 0):
                particles += p.split()
        return particles

    #first split
    children = add_generation(particles)
    #while the split produces any particles
    while(children):
        #append this generation to the list of particles
        particles += children
        #and make a new generation
        children = add_generation(children)
    return particles

#computes the number of bins according to the Freedman-Diaconis rule
def compute_n_bins(x_min, x_max, data):
    #first we need to fetch the 0.25 and 0.75 quantile for the IRQ
    quantiles = array.array('d', [0,0]);
    ROOT.TMath.Quantiles(data.GetSize(), 2, data.GetArray(), quantiles, array.array('d',[0.25, 0.75]), False)
    bins = int((x_max-x_min)/2*pow(float(data.GetSize()), 1.0/3.0)/(quantiles[1]-quantiles[0]))
    return bins


canv = ROOT.TCanvas("canv","plots for SDA course - ROOT 4", 2000, 2000)

#height of first interaction
n_entries = 500000
interactions = ROOT.TArrayD(n_entries)
for i in range(n_entries):
    interactions[i] = gen_interaction(xpp, float('inf'))

x_range = [0, 150000]
h = ROOT.TH1D(
    "height of first split",
    ";height of first split (m); probablity",
    #autocompute the number of bins
    compute_n_bins(x_range[0], x_range[1], interactions),
    x_range[0],
    x_range[1]
)
for x in interactions:
    h.Fill(x)
h.Scale(1.0/(h.GetBinWidth(1)*n_entries))
h.Draw("HIST")
line_x = rho_int_inv(rho_int(float('inf'))-xpp)
print(line_x)
line = ROOT.TLine(line_x,0, line_x, h.GetMaximum())
line.SetLineColor(ROOT.kRed);
line.Draw("SAME")
canv.Update()
canv.SaveAs("a.png")

#energies to simulate
E = [10**5, 10**6, 10**7]
showers = [create_shower(e) for e in E]
z_range = [0,20*10**3]

#plot n_particles vs height
canv.Clear()
canv.Update()
canv.Divide(2,2)
x_range = [0,20]
H = [ROOT.TArrayD(len(shower)) for shower in showers]
for h,shower in zip(H,showers):
    print("computing H for shower")
    for i,p in enumerate(shower):
        h[i] = 0.5*(p.start_pos[2]+p.end_pos[2])/1000

#we use the same number of bins for each shower so that there is a fair comparison to be made.
ROOT.gStyle.SetTitleOffset(1.8,"y") # adjust with histogram.GetYaxis().SetTitleOffset)
n_bins = compute_n_bins(x_range[0], x_range[1], H[1])
print(n_bins)
h2 = [ROOT.TH1D(
    "shower at E = " + ["100 GeV", "1 TeV", "10 TeV"][i],
    ";Height from surface (km);Number of charged particles",
    n_bins,
    x_range[0],
    x_range[1]
) for (i,e),h in zip(enumerate(E),H)]
for hist,h in zip(h2,H):
    print("filling...")
    for x in h:
        hist.Fill(x)
for i,h in enumerate(h2):
    canv.cd(i+1)
    h.Draw()
canv.Update()
canv.SaveAs("c.png")

ROOT.gStyle.SetTitleOffset(1.8,"z") # adjust with histogram.GetYaxis().SetTitleOffset)
#plot showers
canv.Clear()
canv.Divide(2,2)
for i,shower in enumerate(showers):
    canv.cd(i+1)
    cutoff = 85
    plot_shower(shower , "shower at E = " + ["100 GeV", "1 TeV", "10 TeV"][i],
                 xysize = [100,100,100][i],
                 zsize = z_range[1],
                 min_E = cutoff)
canv.Update()
canv.SaveAs("b.png")
#hangs here for some reasons...
print("done with B")


#compute Xmav vs energy
# L = []
# H = []
# means = []
# for e in range(10,14):
#     E = 10**(0.5*e)
#     L = [0.5*(p.start_pos[2]+p.end_pos[2]) for p in create_shower(E)]
#     h2[0].Clear()
#     for l in L:
#         h2[0].Fill(l)
#     means.append(h2[0].GetMean())
# print(means)
