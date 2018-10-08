#roel deckers

import inspect
import ROOT
import array
import math
from copy import copy
from operator import add,mul,sub

#contour drawing from example 2
def draw_contour( mean, cov ) :
    def matrix_mult( m1, m2 ) :

        if type(m2) == float or  type(m2) == int :
            r = copy(m1)
            r*=m2
            return r

        elif type(m2) == ROOT.TVectorD :
            # turn it into a matrix and back again...
            m2_ = ROOT.TMatrixD( m2.GetNrows(),1 , m2.GetMatrixArray())
            prod = matrix_mult(m1,m2_)
            return ROOT.TVectorD( m2.GetNrows(), prod.GetMatrixArray() )

        if  m2.GetNrows() != m1.GetNcols() :
            print "error in matrix multiplication: rows and columns don't match."
            return None

        r=  ROOT.TMatrixD(m1.GetNrows(), m2.GetNcols() )
        r.Mult( m1, m2 )
        return r


    def matrix_rmult( m1, m2 ) :

        if type(m2) == float or  type(m2) == int :
            r = copy(m1)
            r*=m2
            return r

    def matrix_radd( m1, m2 ) :
        return m1+m2

    ROOT.TMatrixD.__mul__  = matrix_mult
    ROOT.TMatrixD.__rmul__ = matrix_rmult # for float * matrix
    #ROOT.TMatrixD.__radd__ = matrix_radd # for float + matrix
    ROOT.TVectorD.__rmul__ = matrix_rmult
    """ Draw an error-ellipse.
        X is a 2-d TVectorD describing the central values.
        C is the 2x2 covariance matrix """
    C = ROOT.TMatrixD(2,2, array.array('d', cov))
    X = ROOT.TVectorD(2, array.array('d', mean))
    M_ = ROOT.TMatrixDEigen( C )
    ME = M_.GetEigenVectors()
    VE = M_.GetEigenValues()

    x,y,V = [],[], ROOT.TVectorD(2)

    # The idea is that there exist coordinates a,b where the
    # elipse is 'horizontal'. That is:
    # (a,b) form an elipse if we let phi run and do
    # a = ea * sin(phi), b = eb * cos(phi), where ea and eb are
    # the errors in these coordinates. They are the eigenvalues
    # of C. The eigenvectors tell us how a and b are related to
    # x and y.

    for phi in range(361):

        V[0] = math.sqrt(VE[0][0]) * math.sin(phi/180.0*math.pi)
        V[1] = math.sqrt(VE[1][1]) * math.cos(phi/180.0*math.pi)

        U = ME * V
        x.append(U[0]+X[0])
        y.append(U[1]+X[1])

    g = ROOT.TGraph( len(x), array.array('d',x), array.array('d',y) )
    g.SetLineWidth(5)
    g.SetLineColor(3)
    g.Draw("LP")
    ROOT.SetOwnership(g, False )
    return g


#numerical optimizers

# MUCH slower than newton_rhapson even if each step is heavier for NR
def momentum_solver(objective, init, delta, abs_limit = 1e-6):
    step = 2*abs_limit
    dimensions = 2
    ddx = [0 for _ in range(dimensions)]
    ddx2 = [0 for _ in range(dimensions)]
    objectives = [[0 for _1 in range(3)] for _2 in range(3)]
    steps = [ROOT.TVectorD(2, array.array('d',init))]
    v = ROOT.TVectorD(2,array.array('d', [0, 0]))
    grad = ROOT.TVectorD(2,array.array('d', [0, 0]))
    while step > abs_limit:
        for i in range(2):
            argp = [steps[-1][0], steps[-1][1]]
            argm = [steps[-1][0], steps[-1][1]]
            argp[i] = argp[i] + delta
            argm[i] = argp[i] - delta
            grad[i]= (objective(*argp)-objective(*argm))/(2*delta)
        v = 0.5*v + grad
        next = steps[-1]-v
        step = grad.Norm2Sqr()
        steps.append(next)
    return steps

#solve a 1D objective with NR
def newton_rhapson_1D(objective, init, delta, abs_limit = 1e-6):
    step = 2*abs_limit
    dimensions = 1
    objectives = [0 for _1 in range(3)]
    steps = [init]
    while step > abs_limit:
        arg = copy(steps[-1])
        for i in range(3):
            darg = arg + (i-1)*delta
            objectives[i] = objective(darg)
        ddx = (objectives[2]-objectives[0])/(2*delta)
        invddx2 = (delta*delta)/(objectives[2]-2*objectives[1]+objectives[0])
        step = ddx*invddx2
        v = 0.25*step
        next = steps[-1]-v
        steps.append(next)
    return (steps, invddx2)

#minimize a 2d objective with NR
def newton_rhapson_2D(objective, init, delta, abs_limit = 1e-6):
    gamma_new = 6.0/5.0
    step = 2*abs_limit
    dimensions = 2
    ddx = [0 for _ in range(dimensions)]
    ddx2 = [0 for _ in range(dimensions)]
    objectives = [[0 for _1 in range(3)] for _2 in range(3)]
    steps = [init]
    v = [0,0]
    while step > abs_limit:
        gamma = gamma_new
        arg = [steps[-1][0], steps[-1][1]]
        for i in range(3):
            arg[1] = steps[-1][1]  + (i-1) * delta
            for j in range(3):
                arg[0] = steps[-1][0] + (j-1) * delta
                objectives[i][j] = objective(*arg)

        ddx = (objectives[1][2]-objectives[1][0])/(2*delta)
        ddy = (objectives[2][1]-objectives[0][1])/(2*delta)

        ddx2 = (objectives[1][2]-2*objectives[1][1]+objectives[1][0])/(delta*delta)
        ddy2 = (objectives[2][1]-2*objectives[1][1]+objectives[0][1])/(delta*delta)
        dxdy = (objectives[2][2]-objectives[2][1]-objectives[1][2]+2*objectives[1][1]-objectives[0][1]-objectives[1][0]+objectives[0][0])/(2*delta*delta)

        invH = [ddy2, -dxdy, -dxdy, ddx2]
        c = (ddx2*ddy2-dxdy*dxdy)*gamma
        if c < 1e-12:
            print "can't invert Hessian matrix safely!"
            return steps
        invH = [h/c for h in invH]
        grad =  [ddx, ddy]
        v = [invH[0]*grad[0]+invH[2]*grad[1], invH[1]*grad[0]+invH[3]*grad[1]]
        if any([math.isnan(x) for x in v]):
            v = [0,0]
            gamma_new = gamma*2
            print "got NaN values, increasing gamma to ", gamma
            if gamma > 10000:
                print "Error: unstable? gamma grew large and still no convergence"
                return (steps, [gamma*h for h in invH])
            if len(steps) > 1:
                steps.pop()
        else:
            next = map(sub, steps[-1], v)
            step = math.sqrt(grad[0]*grad[0]+grad[1]*grad[1])
            steps.append(next)
            gamma_new = max(gamma*0.95, 1)
    return (steps, [gamma*h for h in invH])

#utility wrapper around 1d and 2d.
def newton_rhapson(objective, init, delta, abs_limit = 1e-6):
    dim = len(init)
    if dim == 1:
        return newton_rhapson_1D(objective, init[0], delta, abs_limit)
    elif dim == 2:
        return newton_rhapson_2D(objective, init, delta, abs_limit)
    else:
        print("ERROR: NR not defined for dimensions > 2")
        exit(-1)


#convert back from y1 y3 to a and b
def y1y3_to_ab(y1,y3):
    return  ((y3-y1)/2, (3*y1-y3)/2)
#original model, only valid b > -1/3a (dependent variables)
def model_ab(m,a,b,w = 1):
    return (a*m+b)*w
#coordinate changed model, valid on y1 > 0, y3 > 0 (independent variables)
def model_y1y3(m,y1,y3,w = 1):
    return ((3-m)*y1+(m-1)*y3)*0.5*w

#coordinate changed model, valid on y1 > 0, y3 > 0 (independent variables)
def model_abcd(m,a,b,c,d,w = 1):
    return (a+m*(b+m*(c+m*d)))*w

ROOT.gStyle.SetPadRightMargin(0.2)  # make room for y-title, adjust with pad.SetLeftMargin()
ROOT.gStyle.SetPadLeftMargin(0.125)  # make room for y-title, adjust with pad.SetLeftMargin()
ROOT.gStyle.SetTitleOffset(1.4,"z") # adjust with histogram.GetYaxis().SetTitleOffset)
canv = ROOT.TCanvas("canv","plots for SDA course - ROOT 5", 2000, 2000)
infile = ROOT.TFile("dataset/data.root")

original = infile.Get('hdata')
ROOT.gStyle.SetOptStat(0);
original.SetTitle(";Mass (GeV);Event Count")
original.SetLineWidth(6);
original.Draw("HIST")
canv.Update()
canv.SaveAs("original.png")


#draw the histogram with a guess for a and b drawn on top
def draw_guess(h, title, a, b):
    y = array.array('d', [model_ab(h.GetBinCenter(i), a, b, h.GetBinWidth(1)) for i in range(1,h.GetSize()-1)])
    x = array.array('d', [h.GetBinCenter(i) for i in range(1,h.GetSize()-1)])
    markers = ROOT.TGraphErrors(h.GetSize()-2,x,y)
    markers.SetLineWidth(6)
    markers.SetLineColor(ROOT.kRed)
    h.Draw("HIST")
    markers.Draw("SAME")
    canv.Update()
    canv.SaveAs(title)
    canv.Clear()

#get the likelihood that bin i of hist h is x
def get_log_likelihood(k, x):
    #p(k) = e^(-x)*x^k/k!
    #log(p(k)) = -x + log(x^k) - log(k!)
    #  = -x + k*log(x) - LnGamma(k+1)
    return -x + k*ROOT.TMath.Log(x) - ROOT.TMath.LnGamma(k+1)

def get_total_log_likelihood_ab(h,a,b):
    return sum([get_log_likelihood(h.GetBinContent(i), model_ab(h.GetBinCenter(i), a, b, h.GetBinWidth(1)))  for i in range(1,h.GetSize()-1)])

#gets the log likelihood that a and b match the data in h.
def get_total_log_likelihood_y1y3(h,y1,y3):
    return sum([get_log_likelihood(h.GetBinContent(i), model_y1y3(h.GetBinCenter(i), y1, y3, h.GetBinWidth(1)))  for i in range(1,h.GetSize()-1)])

#gets the log likelihood that a and b match the data in h.
def get_total_log_likelihood_abcd(h,a,b,c,d):
    return sum([get_log_likelihood(h.GetBinContent(i), model_abcd(h.GetBinCenter(i), a,b,c,d, h.GetBinWidth(1)))  for i in range(1,h.GetSize()-1)])

#objective function to minimize for finding optimal y1, y3
def objective_func(y1,y3):
    return -get_total_log_likelihood_y1y3(original, y1, y3)

def objective_func_ab(a,b):
    return -get_total_log_likelihood_ab(original, a, b)

def objective_func_abcd(a,b,c,d):
    return -get_total_log_likelihood_abcd(original, a, b, c, d)

def objective_func_b(b):
    return -get_total_log_likelihood_ab(original, 0, b)

(minpath,cov) = newton_rhapson(objective_func_b, init = [200], delta = 1e-4, abs_limit = 1e-12)
print "b only guess:", minpath[-1], "+-", math.sqrt(cov)
draw_guess(original, "b.png", 0,minpath[-1])


def make_2d_hist(h,title, y1, y3, init = [0,111]):
    max = -float('inf')
    min = -max
    h2d = ROOT.TH2D(title, title, y1[0], y1[1], y1[2], y3[0], y3[1], y3[2])
    print h2d.GetXaxis().GetNbins(), h2d.GetYaxis().GetNbins()
    for i in range(1,h2d.GetXaxis().GetNbins()+1):
        for j in range(1,h2d.GetYaxis().GetNbins()+1):
            a = h2d.GetXaxis().GetBinCenter(i)
            b = h2d.GetYaxis().GetBinCenter(j)
            z = -get_total_log_likelihood_ab(h, a, b)
            if not math.isnan(z):
                if z > max:
                    max = z
                if z < min:
                    min = z
                h2d.SetBinContent(i,j,z)
    h2d.SetMinimum(min)
    h2d.SetMaximum(max)
    h2d.SetTitle(";a (dimensionless);b (GeV); -log L_{M_{a,b}}")
    print min, max
    canv.Clear()
    h2d.Draw("COLZ1")
    canv.Update()
    (minpath,cov) = newton_rhapson(objective_func_ab, init = init, delta = 1e-4, abs_limit = 1e-14)
    a_path = array.array('d', [m[0] for m in minpath])
    b_path = array.array('d', [m[1] for m in minpath])
    markers = ROOT.TGraphErrors(len(minpath),a_path,b_path)
    markers.SetLineWidth(3)
    markers.SetLineColor(ROOT.kBlack)
    markers.SetMarkerStyle(21);
    markers.SetMarkerColor(ROOT.kRed);
    markers.Draw("LP")
    draw_contour(minpath[-1], cov)
    canv.SaveAs(title)
    return (minpath,cov)

(minpath, cov) = make_2d_hist(original, "2d.png", [50,-100,10], [50, 100,350])
print len(minpath)
(minpath, cov) = make_2d_hist(original, "2d2.png", [100,minpath[-1][0]-15,minpath[-1][0]+15], [100, minpath[-1][1]-30,minpath[-1][1]+30], [minpath[-1][0]+6, minpath[-1][1]+20])
(a,b) = (minpath[-1][0], minpath[-1][1])
print cov
print a, " +- ", math.sqrt(cov[0])
print b, " +- ", math.sqrt(cov[3])
draw_guess(original, "best_guess.png", a,b)

def minimize(function, start_values = None, ranges = None, maxcalls = 10000, tollerance = 1e-6):
    args = inspect.getargspec(function)[0]
    npar = len(args)
    if not start_values : start_values = [0.0] * npar
    if not ranges : ranges = [[0,0]]*npar
    step_size = [0.1]*npar
    minuit = ROOT.TMinuit(npar)
    for i in range(npar):
        ierflg = ROOT.Long()
        minuit.mnparm(i, args[i], start_values[i], step_size[i], ranges[i][0], ranges[i][1], ierflg)

    def fcn( _npar, gin, f, par, iflag):
        parlist = [par[i] for i in range(npar)]
        f[0]= function(*parlist)

    minuit.SetFCN(fcn)
    arglist = array.array('d', [maxcalls, tollerance])
    minuit.mnexcm("MIGRAD", arglist, len(arglist), ierflg)
    r, errs = [],[]
    for i in range(npar):
        r.append( ROOT.Double())
        errs.append(ROOT.Double())
        minuit.GetParameter(i, r[-1],errs[-1])
    return r,errs

r,errs = minimize(objective_func_ab, start_values = [0,111])
#gives the same fit values as our minimizer however the errors are different
print r

r,errs = minimize(objective_func_abcd, start_values = [1,1,1,1])
print r

#draw the histogram with a guess for abcd drawn on top
def draw_guess(h, title, a, b, c, d):
    y = array.array('d', [model_abcd(h.GetBinCenter(i), a, b, c, d, h.GetBinWidth(1)) for i in range(1,h.GetSize()-1)])
    x = array.array('d', [h.GetBinCenter(i) for i in range(1,h.GetSize()-1)])
    markers = ROOT.TGraphErrors(h.GetSize()-2,x,y)
    markers.SetLineWidth(6)
    markers.SetLineColor(ROOT.kRed)
    h.Draw("HIST")
    markers.Draw("SAME")
    canv.Update()
    canv.SaveAs(title)
    canv.Clear()

draw_guess(original, "3rd_degree.png", r[0], r[1], r[2], r[3])
