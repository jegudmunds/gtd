import numpy as np
# import scipy as sci
# from scipy import fftpack
# from scipy import stats
# from scipy.interpolate import splprep, splev
# from scipy import optimize as optimize
# import scipy.special as ss


'''
Collection of generatl GTD functions

'''

### Global variables
dtr = np.pi/180.
cee = 3.e8 #m/sec


def gtd_airy(x):
    arg = -x*np.power(3.,-1./3.)
    gtd_afa, gtd_dafa, gtd_afb, gtd_dafb = ss.airy(arg)
    gtd_afa = np.pi*np.power(3.,-1./3.)*gtd_afa
    gtd_dafa = -np.pi*np.power(3.,-2./3.)*gtd_dafa
    return gtd_afa, gtd_dafa


def airy_a(x):
    afa, dafa = gtd_airy(x)
    xmin = afa*afa
    return xmin

def airy_da(x):
    afa, dafa  = gtd_airy(x)
    xmin = dafa*dafa
    return xmin

def norm_vec(xc, yc, zc):

    leng = leng = np.sqrt((xc[1]-xc[0])**2+(yc[1]-yc[0])**2+(zc[1]-zc[0])**2)
    nv = np.array([xc[1]-xc[0], yc[1]-yc[0], zc[1]-zc[0]])/leng

    return nv

### The following three functions are for walking around the top of the rim.
def ang2_3pts(theta, gs_diam, gs_dist, ap_diam):
    '''
    Given an angle from the center of the aperture, return three points on the
    ground screen rim that correspond to a projection.

    See diagram in JEG's notebook to better understand geometry
    '''

    th_c = np.pi - theta
    rrr = gs_diam/2.
    apr = ap_diam/2.
    bbb = rrr - gs_dist # Distance from gs center to aperture
    th2 = np.arcsin(bbb*np.sin(th_c)/rrr) #Angle between radii at gs rim
    th_gc = np.pi - th_c - th2 #Angle from gs center to rim
    r_prime = rrr * np.sin(th_gc)/np.sin(th_c)
    if (theta == 0):
        r_prime = gs_dist

    xgs_c = r_prime*np.cos(theta)+bbb
    ygs_c = r_prime*np.sin(theta)
    slo = np.tan(theta)
    # print 'theta and slope',theta/dtr,slo
    off_t = -bbb*np.tan(theta) + apr
    off_b = -bbb*np.tan(theta) - apr
    pts = np.array([xgs_c, ygs_c, r_prime, slo, off_b, off_t, th2, th_gc] )

    return pts


def circline2xy(slo, off, rgs):
    '''
    Given the slope, offset and radius of the GS, find the intercept on the top of the GS.
    '''

    # print slo,off,rgs
    aa = 1 + slo * slo
    bb =  2. * slo * off
    cc = off*off - rgs*rgs
    aa2 = aa
    bb2 = -2. * off
    cc2 = off * off - rgs*rgs*slo*slo
    xc = (-bb + np.sqrt(bb*bb-4.*aa*cc))/(2.*aa)
    yc = (-bb2 + np.sqrt(bb2*bb2-4.*aa2*cc2))/(2.*aa2)

    return xc, yc

def toppt2ang(xt, yt, gs_diam, gs_dist, ap_diam):
    '''
    Go from the top point to an angle for the center of the next point up.
    '''

    xdist = xt - (gs_diam/2.-gs_dist)
    theta = np.arctan((yt+ap_diam/2.)/xdist)

    return theta

#
# Now make a plot from the top. Now the y-axis in the above becomes the z-axis and is out of
# the page. This loop sets up the geometry for the integral so it is mandatory.