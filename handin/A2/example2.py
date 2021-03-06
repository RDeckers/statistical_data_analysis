#Roel Deckers
import ROOT
import array
from math import *
import matrices

from matrices import make_TMatrixD
from matrices import make_TVectorD

#seems to be missing from matrices.py on the server even though it was pulled in before?
#add it just in case
def matrix_transposed( M ) :
    return ROOT.TMatrixD(ROOT.TMatrixD.kTransposed, M)

x    = [ 812., 1800., 400., 464., 818., 356., 289., 301., 164., 393. ]
errx = [ 41.,  50.,   33.,  35.,  33.,  22.,  27.,  29.,  19.,  24.  ]
for i in range( len(x) ):
    print "x%d  = %f +/- %f " % (i,x[i],errx[i] )

#uncorrelated, so we can use this formula directly for the expectation values
def R_gen(x):
    return make_TVectorD(4, [x[2]/x[0], x[3]/x[1], (x[6]-x[8]*x[2]/x[0])/x[4], (x[7] - x[9] * x[3]/x[1])/x[5]])
R = R_gen(x)
#we can use basic formulas for propagating the errors onto R, but later we need the Covariance for
# computing the errors on u/d
# this is the R used in the RC(R^T) exression for propagated errors
R_matrix = matrices.numerical_derivative(R_gen, x)
#create a (zero-initialized) matrix
Cx = ROOT.TMatrixD(10,10)
for i in range( len(x) ):
    Cx[i][i] = errx[i]*errx[i] #set it with the variance on the diagonal
CR = R_matrix*Cx*matrix_transposed(R_matrix)
print CR
for i in range( len(R) ):
    # The variance of the elements of R are on the diagonal of CR
    print "R%d  = %f +/- %f " % (i,R[i],sqrt(CR[i][i]) )

H =  make_TMatrixD(4,4, [  0.675 , -0.607 , -0.119 ,  0.010,
       -0.282 ,  1.331 ,  0.027 , -0.049,
       -0.133 ,  0.060 ,  0.477 , -0.078,
        0.024 , -0.299 , -0.186 ,  0.185 ])
#propagate the (correlated) errors of R
Cud = H*CR*matrix_transposed(H)
print Cud
ud = H*R
for i in range( len(ud) ):
    print "ud%d  = %f +/- %f " % (i,ud[i],sqrt(Cud[i][i]) )

#make a new canvas, and draw the contours
ROOT.gStyle.SetPadLeftMargin(0.12)  # make room for y-title
canv = ROOT.TCanvas("canv","plots for SDA course - ROOT 2", 800, 800)
g = matrices.draw_contour( make_TVectorD(2,   [ud[0],ud[1]]),
                       make_TMatrixD(2,2, [Cud[0][0],Cud[0][1],
                                           Cud[1][0],Cud[1][1]] ))
g.SetTitle(";u_{L}^{2};d_{L}^{2};");
#add a point with the errors in ud[0] and ud[1] to the graph
n1 = 1;
x1  = array.array('d', [ud[0]])
y1  = array.array('d', [ud[1]])
ex1 = array.array('d', [sqrt(Cud[0][0])])
ey1 = array.array('d', [sqrt(Cud[1][1])])
gr1 = ROOT.TGraphErrors(n1,x1,y1,ex1,ey1)
gr1.SetMarkerColor(ROOT.kRed)
gr1.SetMarkerStyle(21)
gr1.Draw("P")
canv.Update()
canv.SaveAs("result.png")
