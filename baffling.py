'''
An attempt at framing Lyman's GTD calculations using classes.
This helps us identify what properties need to be shared between objects
and what attributes can remain private.
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import time
import argparse as ap
import numpy as np
import gtd

opj = os.path.join

cee = 2.998e8 #m/sec
freq_90 = 90. #GHz
freq_150 = 150. #GHz
wl_90 = cee/(freq_90*1.e9)
wl_150 = cee/(freq_150*1.e9)
kwave_90 = 2.*np.pi/wl_90
kwave_150 = 2.*np.pi/wl_150

def gain_model(nu):
    if nu == 90:
        return 8e-4
    elif nu == 150:
        return 5.14e-4

def kwave(freq):

    wl = cee/(freq * 1.e9)

    return 2.*np.pi/wl

class GroundScreen(object):

    def __init__(self, diameter=7.5, height=4.2, slope=10, units='m'):
        '''
        Ground screen object has height, diameter, and slope
        '''

        self.diameter = diameter
        self.r = diameter/2.0
        self.height = height
        self.slope = slope
        self.units = units

    @property
    def diameter(self):
        return self._diameter

    @diameter.setter
    def diameter(self, value):
        self._diameter = value


    @property
    def height(self):
        return self._height

    @height.setter
    def height(self, value):
        self._height = value


    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, value):
        self._units = value


    def __str__(self):

        return 'Ground screen object with diameter {} {} and height {} {}'.\
            format(self.diameter, self.units, self.height, self.units)


    def circline2xy():
        '''
        Given the slope, offset, and radius of the GS, find the intercept on
        the top of the GS.
        '''

        aa = 1 + self.slope * self.slope
        bb =  2. * self.slope * off
        cc = off*off - (self.diameter/2.)*(self.diameter/2.)
        aa2 = aa
        bb2 = -2. * off
        cc2 = off * off - (self.diameter/2.)*(self.diameter/2.)*self.slope*self.slope
        xc = (-bb + np.sqrt(bb*bb-4.*aa*cc))/(2.*aa)
        yc = (-bb2 + np.sqrt(bb2*bb2-4.*aa2*cc2))/(2.*aa2)

        return xc, yc

class Telescope(object):

    def __init__(self, x=0.9, y=2.8, aperture=0.42, elevation=50., fov=35.,
        units='m', fwhm_rt=0.5):
        '''
        Telescope object, has a location relative to the center of the ground
        screen, aperture diameter, elevation, and a field of view (FOV)
        '''

        self.x = x
        self.y = y

        self.aperture = aperture
        self.r = aperture/2.0
        self.elevation = np.radians(elevation)
        self.fov = np.radians(fov)
        self.fwhm_rt = np.radians(fwhm_rt)
        self.units = units

    def __str__(self):
        pass

class Forebaffle(object):

    def __init__(self, length=1, angle=20, units='m', diameter=0.5,
        Temperature=300, reflective=True, Nl=11, Nphi):
        '''
        Mostly self explainatory attributes
        diameter : diameter of baffle at the base (closest to aperture)

        Nl   : number of meshing elements along the length of the baffle
               (int)
        Nphi : number of meshing elements along the azimuthal direction
               (int)

        '''



        self.length = length
        self.angle = angle
        self.diameter = diameter
        self.units = units

    @property
    def length(self):
        return self._length

    @length.setter
    def length(self, value):
        self._length = value

    @property
    def angle(self):
        return self._angle

    @angle.setter
    def angle(self, value):
        self._angle = value

    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, value):
        self._units = value

    @property
    def diameter(self):
        return self._diameter

    @diameter.setter
    def diameter(self, value):
        self._diameter = value


    def angle2baffletip():
        '''

        Calculates the anglue from the center of the aperture to the tip of the baffle

        '''

        pass




    def __str__(self):
        pass


class BaffleGeometry(object):
    '''
    This class combines information from a GroundScreen and Telescope object
    and calculates various geometrical aspects that help us draw figures of the
    setup and calculate various properties that depend on the joint
    geometrical relations.

    ToDo: Add forebaffle
    '''

    def __init__(self, gs=None, t=None, oa_len=4, alpha=45, nth=361,
        ngst=26, delta_h=0.5, freqs=[90, 150], verbose=False):
        '''
        gs: Ground screen object
        t: Telescope object
        oa_len: length for lines for edge "rays." Just for plotting purposes.

        ngst: Number of integration steps along the surface of the ground screen

        '''

        self.freqs = freqs
        self.nfreqs = len(freqs)
        self.gs = GroundScreen() if gs is None else gs
        self.t = Telescope() if t is None else t

        # Defining geometries that help frame the system
        self.oa_len = oa_len

        # print(self.gs.r)
        # print(self.gs.height)

        self.gs_dist = self.gs.diameter/2.0 - self.t.x

        self.xang = np.pi/2 - self.t.elevation
        self.xlenap = np.cos(self.xang) * self.t.aperture
        self.ylenap = np.sin(self.xang) * self.t.aperture

        self.x_ap = np.array([self.t.x - self.xlenap/2., self.t.x + self.xlenap/2.])
        self.y_ap = np.array([self.t.y + self.ylenap/2., self.t.y - self.ylenap/2.])

        self.x_oa = np.array([self.t.x, self.t.x \
            + np.cos(self.t.elevation)*oa_len])

        self.y_oa = np.array([self.t.y, self.t.y \
            + np.sin(self.t.elevation)*oa_len])

        self.rt_ang = self.t.elevation + self.t.fov/2. #top ray angle
        self.rb_ang = self.t.elevation - self.t.fov/2. #bottom ray angle

        # Top
        self.x_oa_t = np.array([self.t.x - self.xlenap/2., self.t.x - \
            self.xlenap/2. + self.oa_len * np.cos(self.rt_ang)])
        self.y_oa_t = np.array([self.t.y + self.ylenap/2., self.t.y + \
            self.ylenap/2. + self.oa_len * np.sin(self.rt_ang)])

        # Bottom
        self.x_oa_b = np.array([self.t.x + self.xlenap/2., self.t.x + \
            self.xlenap/2. + self.oa_len * np.cos(self.rb_ang)])
        self.y_oa_b = np.array([self.t.y - self.ylenap/2., self.t.y - \
            self.ylenap/2. + self.oa_len * np.sin(self.rb_ang)])

        self.x_gs = np.array([self.gs_dist + self.t.x, self.gs_dist + self.t.x])
        self.y_gs = np.array([0, self.gs.height])

        # Make the outlines of three beams, cyan, green, magenta
        rc_ang = self.t.elevation + self.t.fov/2. #top ray angle
        x_oa_c = np.array([self.t.x, self.t.x + np.cos(rc_ang)*oa_len])
        y_oa_c = np.array([self.t.y, self.t.y + np.sin(rc_ang)*oa_len])
        rt_ang = self.t.elevation + self.t.fov/2. + self.t.fwhm_rt/2. #top ray angle
        x_oa_t_c = np.array([self.t.x-self.xlenap/2., self.t.x-self.xlenap/2.+oa_len*np.cos(rt_ang)])
        y_oa_t_c = np.array([self.t.y+self.ylenap/2., self.t.y+self.ylenap/2.+oa_len*np.sin(rt_ang)])
        rb_ang = self.t.elevation + self.t.fov/2.-self.t.fwhm_rt/2. #bottom ray angle
        x_oa_b_c = np.array([self.t.x+self.xlenap/2., self.t.x+self.xlenap/2.+oa_len*np.cos(rb_ang)])
        y_oa_b_c = np.array([self.t.y-self.ylenap/2., self.t.y-self.ylenap/2.+oa_len*np.sin(rb_ang)])

        rc_ang = self.t.elevation #top ray angle
        x_oa_g = np.array([self.t.x, self.t.x + np.cos(rc_ang)*oa_len])
        y_oa_g = np.array([self.t.y, self.t.y + np.sin(rc_ang)*oa_len])
        rt_ang = self.t.elevation + self.t.fwhm_rt/2. #top ray angle
        x_oa_t_g = np.array([self.t.x-self.xlenap/2., self.t.x-self.xlenap/2.+oa_len*np.cos(rt_ang)])
        y_oa_t_g = np.array([self.t.y+self.ylenap/2., self.t.y+self.ylenap/2.+oa_len*np.sin(rt_ang)])
        rb_ang = self.t.elevation - self.t.fwhm_rt/2. #bottom ray angle
        x_oa_b_g = np.array([self.t.x+self.xlenap/2., self.t.x+self.xlenap/2.+oa_len*np.cos(rb_ang)])
        y_oa_b_g = np.array([self.t.y-self.ylenap/2., self.t.y-self.ylenap/2.+oa_len*np.sin(rb_ang)])

        rc_ang = self.t.elevation - self.t.fov/2. #top ray angle
        x_oa_m = np.array([self.t.x, self.t.x + np.cos(rc_ang)*oa_len])
        y_oa_m = np.array([self.t.y, self.t.y + np.sin(rc_ang)*oa_len])
        rt_ang = self.t.elevation-self.t.fov/2. + self.t.fwhm_rt/2. #top ray angle
        x_oa_t_m = np.array([self.t.x-self.xlenap/2., self.t.x-self.xlenap/2.+oa_len*np.cos(rt_ang)])
        y_oa_t_m = np.array([self.t.y+self.ylenap/2., self.t.y+self.ylenap/2.+oa_len*np.sin(rt_ang)])
        rb_ang = self.t.elevation-self.t.fov/2.-self.t.fwhm_rt/2. #bottom ray angle
        x_oa_b_m = np.array([self.t.x+self.xlenap/2., self.t.x+self.xlenap/2.+oa_len*np.cos(rb_ang)])
        y_oa_b_m = np.array([self.t.y-self.ylenap/2., self.t.y-self.ylenap/2.+oa_len*np.sin(rb_ang)])

        # Alpha
        self.alpha = np.radians(-alpha)
        self.lena = self.y_gs[1]/np.sin(-self.alpha)
        self.x_inc = np.array([self.x_gs[1] + self.lena * np.cos(self.alpha), self.x_gs[1]])
        self.y_inc = np.array([0, self.y_gs[1]]) #start all incident rays on the ground

        # Find the distance and angle of the center of the aperture to the top of the ground screen
        self.x_gs_ap = np.array([self.x_gs[1], self.t.x])
        self.y_gs_ap = np.array([self.y_gs[1], self.t.y])
        self.rr = np.sqrt(np.square(self.x_gs_ap[0] - self.x_gs_ap[1]) + \
            np.square(self.y_gs_ap[0] - self.y_gs_ap[1]) )

        self.phi = np.arctan((self.y_gs_ap[1] - self.y_gs_ap[0])\
            /(self.x_gs_ap[1]-self.x_gs_ap[0] )) + np.pi

        if verbose:
            print('Nominal phi (deg): {}, Nominal s (m): {:.3f}'.format(np.degrees(self.phi), self.rr))

        # Aperture.
        self.th = np.radians(np.arange(nth))
        self.xl = self.t.x + (self.t.r) * np.sin(self.th) * np.sin(self.t.elevation)
        self.yl = (self.t.r) * np.cos(self.th)
        self.zl = self.t.y

        # Ground screen
        self.x_gs_rim = self.gs.r * np.cos(self.th)
        self.y_gs_rim = self.gs.r * np.sin(self.th)
        self.z_gs_rim = self.gs.height

        self.ngst = ngst
        self.delta_h = delta_h
        self.npatch = np.int(self.gs.height/self.delta_h) - 2 # don't go all the way down to the ground
        self.gs_hts = self.gs.height - self.delta_h/2. - self.delta_h * np.arange(self.npatch)
        self.dhts = np.repeat(self.delta_h, self.npatch)

        #dist from ap center to ground screen
        self.ppp = self.gs.r - self.t.x

        # Set up useful unit vectors along three optical axes for the following.
        # For Cyan, Green, and Magenta rays
        xl_oa_c = x_oa_c
        yl_oa_c = np.array([0,0])
        zl_oa_c = y_oa_c
        self.oa_c_vec = gtd.norm_vec(xl_oa_c, yl_oa_c, zl_oa_c)

        xl_oa_g = x_oa_g
        yl_oa_g = np.array([0,0])
        zl_oa_g = y_oa_g
        self.oa_g_vec = gtd.norm_vec(xl_oa_g, yl_oa_g, zl_oa_g)

        xl_oa_m = x_oa_m
        yl_oa_m = np.array([0,0])
        zl_oa_m = y_oa_m
        self.oa_m_vec = gtd.norm_vec(xl_oa_m, yl_oa_m, zl_oa_m)

        oa_vecs = np.zeros((3, 3))
        oa_vecs[0,:] = self.oa_c_vec
        oa_vecs[1,:] = self.oa_g_vec
        oa_vecs[2,:] = self.oa_m_vec
        self.oa_vecs = oa_vecs

    def irimwalk(self, plot_points=False, fignum=2, verbose=False):
        '''
        Walk along the top of the ground screen going counter clockwise to come up
        with boundaries where the aperture projects onto the rim. For each point
        we get a projected width of the cylinder of waves that hit the aperture.
        By construction, each cylinder fully intersects with the aperture. The
        integral is performed for each of the following, but the integral is not heavy.

        In this option we also support computing the emission from the ground screen

        Exposed parameters:

        gs_sa : Ground screen solid angle

        psi_c, psi_g, psi_m : The angle to the top of the ground screen for three
            different detector beams corresponding to min, mid, and max elevations

        rim_ang : The angle of incidence of a ray from the aperture to the rim
            (see th2 in JEG's notebook)

        rim_fac : The relative projection of the cylinder beam on the rim.
            What's the difference between this and the rim_fac defined in iuconst?

        gs_psi :

        sss :

        ang_to_gs : I think this corresponds to theta in JEG's notebook

        '''

        ptsgs = np.zeros([self.ngst, 8])
        trios = np.zeros([self.ngst, 6])
        cylwidths = np.zeros(self.ngst) #The projected widths of the cylinders of rays
        self.sss = np.zeros(self.ngst)
        phi = np.zeros(self.ngst)
        self.psi_c = np.zeros(self.ngst)
        self.psi_g = np.zeros(self.ngst)
        self.psi_m = np.zeros(self.ngst)
        self.rim_fac = np.zeros(self.ngst)
        self.rim_ang = np.zeros(self.ngst)
        thgsxywalk = np.zeros(self.ngst + 1)
        ang_to_gs = np.zeros(self.ngst)
        gs_areas = np.zeros([self.ngst, self.npatch])
        self.gs_sa = np.zeros([self.ngst, self.npatch])
        gs_surf_thinc = np.zeros([self.ngst, self.npatch])
        gs_norm = np.zeros([self.ngst, 3, self.npatch])
        gs_pos = np.zeros([self.ngst, 3, self.npatch])
        gs_sss = np.zeros([self.ngst, self.npatch]) #distances to patch centers
        self.gs_psi = np.zeros([self.ngst, 3, self.npatch]) # three psi's for cyan, green, magenta

        for i in range(self.ngst):

            theta = thgsxywalk[i]
            ptsgs[i, :] = gtd.ang2_3pts(theta, self.gs.diameter, self.gs_dist,
                self.t.aperture)

            xgs_b, ygs_b = gtd.circline2xy(ptsgs[i, 3], ptsgs[i, 4], self.gs.r)
            xgs_t, ygs_t = gtd.circline2xy(ptsgs[i, 3], ptsgs[i, 5], self.gs.r)

            trios[i, 0] = xgs_b
            trios[i, 1] = ygs_b
            trios[i, 2] = ptsgs[i, 0]
            trios[i, 3] = ptsgs[i, 1]
            trios[i, 4] = xgs_t
            trios[i, 5] = ygs_t
            self.rim_ang[i] = ptsgs[i, 6]

            if verbose:
                print('{:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.1f}'.\
                    format(trios[i, 0], trios[i, 1], trios[i, 2], trios[i, 3],
                        trios[i, 4], trios[i, 5], np.degrees(thgsxywalk[i])))

            if plot_points:
                plt.figure(2)
                plt.plot(trios[i, 2], trios[i, 3], 'go')

            # The following are cylinder widths along the GS. To get to perp need to
            # multiply by np.cos(rim_ang) which is very near unity.
            cylwidths[i] = np.sqrt((trios[i, 5] - trios[i, 1])**2 + \
                (trios[i, 4] - trios[i, 0])**2)

            self.rim_fac[i]= cylwidths[i] / self.t.aperture
            rad_dist = np.sqrt((trios[i, 2] - self.t.x)**2 + trios[i, 3]**2)
            self.sss[i] = np.sqrt(rad_dist**2 + (self.gs.height - self.t.y)**2)
            phi[i] = np.pi + np.arctan((self.gs.height - self.t.y)/rad_dist) # Same as phi in Fig 1

            unit_v_2_rimx = (trios[i,2] - self.t.x)/self.sss[i]
            unit_v_2_rimy = trios[i, 3]/self.sss[i]
            unit_v_2_rimz = (self.gs.height - self.t.y)/self.sss[i]

            cosang_c = self.oa_c_vec[0]*unit_v_2_rimx + self.oa_c_vec[1]*unit_v_2_rimy + \
                       self.oa_c_vec[2]*unit_v_2_rimz
            self.psi_c[i] = np.arccos(cosang_c)

            cosang_g = self.oa_g_vec[0]*unit_v_2_rimx + self.oa_g_vec[1]*unit_v_2_rimy + \
                       self.oa_g_vec[2]*unit_v_2_rimz
            self.psi_g[i] = np.arccos(cosang_g)

            cosang_m = self.oa_m_vec[0]*unit_v_2_rimx + self.oa_m_vec[1]*unit_v_2_rimy + \
                       self.oa_m_vec[2]*unit_v_2_rimz
            self.psi_m[i] = np.arccos(cosang_m)

            # Now get the information for the emission from the ground screen.
            # The last array element is for the different heights of the ground screen.
            # In this case we want the integral of the solid angle.
            areas = cylwidths[i]*np.cos(self.rim_ang[i])*self.dhts
            gs_areas[i, :] = areas
            gs_norm[i, 0, :] = np.cos(ptsgs[i,7]) #outward normal.
            gs_norm[i, 1, :] = np.sin(ptsgs[i,7])
            gs_norm[i, 2, :] = 0. #The screen is vertical
            gs_sss[i, :] = np.sqrt(rad_dist**2 + (self.gs_hts[:]-self.t.y)**2)
            gs_pos[i, 0, :] = (trios[i,2] - self.t.x)/gs_sss[i]  #unit vec towards gs
            gs_pos[i, 1, :] = trios[i,3]/gs_sss[i]
            gs_pos[i, 2, :] = (self.gs_hts[:] - self.t.y)/gs_sss[i]

            for j in range(3):
               cosang = self.oa_vecs[j,0]*gs_pos[i,0,:] + self.oa_vecs[j,1]*gs_pos[i,1,:] + \
                       self.oa_vecs[j,2]*gs_pos[i,2,:]
               self.gs_psi[i, j, :] = np.arccos(cosang) # three psi's for cyan, green, magenta

            self.gs_sa[i,:] = gs_areas[i,:]/np.square(gs_sss[i,:]) #solid angles
            cosang = gs_norm[i,0,:]*gs_pos[i,0,:] + gs_norm[i,1,:]*gs_pos[i,1,:] + \
                       gs_norm[i,2,:]*gs_pos[i,2,:]
            gs_surf_thinc[i,:] = np.arccos(cosang) # three psi's for cyan, green, magenta

            # This is the place to put in the projections for
            # e perpendicular and e parallel to tne plane of incidence. We can
            # make a guess that the net effect will have E_perp in the circumferential
            # direction and e parallel will be vertical.
            #
            # Get the next angle
            thgsxywalk[i+1] = gtd.toppt2ang(trios[i,4], trios[i,5], self.gs.diameter,
                self.gs_dist, self.t.aperture)

        self.ang_to_gs = thgsxywalk[:self.ngst]




    def iuconst(self, fignum=2, add2plot=False):
        '''
        Now make a construction that finds sections of cylindrical wave fronts.
        Maked by the aperture width

        Integrate using the construction of a single line. The alternative is
        the contribution from walking around the rim

        phi : The angle from the top of the ground screen to the aperture center
        (see Figure 1/2 in paper manuscript.)

        rim_fac :
        psi_m :
        psi_g :
        psi_c :

        '''

        plt.figure(fignum)

        x_const = np.array([self.gs.r, self.gs.r])
        y_const = np.array([-10.*self.gs.diameter, 10.*self.gs.diameter])

        if add2plot:
            plt.plot(x_const, y_const, color='grey', lw =0.5)

        gst_const = np.arange(self.ngst) #ground screen top construction
        yc_const = gst_const * self.t.aperture
        yc_top_const =  self.t.r + gst_const * self.t.aperture
        yc_bot_const = -self.t.r + gst_const * self.t.aperture
        xc_const = np.zeros(self.ngst)
        xc_const[0:] = self.gs.r

        if add2plot:
            plt.plot(xc_const, yc_const, 'bo')
            # plt.plot(xc_const, yc_bot_const, 'go')

        ang1_c = np.arctan(yc_const / self.ppp)
        ang_to_gs = ang1_c
        opang1_c = np.pi - ang1_c
        ang2_c = np.arcsin(np.sin(opang1_c) * self.t.x / self.gs.r)
        ang3_c = np.pi - opang1_c - ang2_c
        x_rim_c = self.gs.r * np.cos(ang3_c)
        y_rim_c = self.gs.r * np.sin(ang3_c)

        # Also get the positions on the rim for the bounding rays (top/bot)
        ang1_b   = np.arctan(yc_bot_const / self.ppp)
        opang1_b = np.pi - ang1_b
        ang2_b   = np.arcsin(np.sin(opang1_b) * self.t.x / self.gs.r)
        ang3_b   = np.pi - opang1_b - ang2_b
        x_rim_b  = self.gs.r * np.cos(ang3_b)
        y_rim_b  = self.gs.r * np.sin(ang3_b)
        ang1_t   = np.arctan(yc_top_const/self.ppp)
        opang1_t = np.pi - ang1_t
        ang2_t   = np.arcsin(np.sin(opang1_t) * self.t.x / self.gs.r)
        ang3_t   = np.pi - opang1_t - ang2_t
        x_rim_t  = self.gs.r * np.cos(ang3_t)
        y_rim_t  = self.gs.r * np.sin(ang3_t)
        rim_len  = np.sqrt((x_rim_t - x_rim_b)**2 + (y_rim_t - y_rim_b)**2)
        self.rim_fac  = rim_len / self.t.aperture
        rim_ang  = ang2_c
        gamma    = ang1_c #Can use the above + geometry or this gamma from the construction.
        rad_dist = np.sqrt((x_rim_c - self.t.x)**2 + y_rim_c**2)
        sss = np.sqrt(rad_dist**2 + (self.gs.height - self.t.y)**2)
        self.phi = np.pi + np.arctan((self.gs.height - self.t.y) / rad_dist)

        # Now find the angle from the three optical axes.
        # We called them cyan (c), green(g), and magenta (m)
        unit_v_2_rimx = (x_rim_c - self.t.x)/sss
        unit_v_2_rimy = y_rim_c / sss
        unit_v_2_rimz = (self.gs.height - self.t.y) / sss

        cosang_c = self.oa_c_vec[0] * unit_v_2_rimx + self.oa_c_vec[1] * unit_v_2_rimy + \
            self.oa_c_vec[2] * unit_v_2_rimz
        self.psi_c = np.arccos(cosang_c)
        cosang_g = self.oa_g_vec[0] * unit_v_2_rimx + self.oa_g_vec[1] * unit_v_2_rimy + \
            self.oa_g_vec[2] * unit_v_2_rimz
        self.psi_g = np.arccos(cosang_g)
        cosang_m = self.oa_m_vec[0] * unit_v_2_rimx + self.oa_m_vec[1] * unit_v_2_rimy + \
            self.oa_m_vec[2] * unit_v_2_rimz
        self.psi_m = np.arccos(cosang_m)

        # plt.plot(x_rim_c, y_rim_c, 'ro')
        if add2plot:
            for k in range(self.ngst):
                xg = np.array([self.t.x, xc_const[k]])
                yg = np.array([0., yc_const[k]])
                plt.plot(xg, yg, color='grey', lw =.5 )


    def pickup(self, nang=85, phip=206, alpha=30, verbose=False):
        '''
        Calculating power pickup

        fss_u0  : Geometrical scattering factor (u=0)
        fss_du0 : Geometrical scattering factor (du=0)
        alp_s : The alpha in LP's Eq. 2 from gtd_3d_screen_4_sac.pdf, see also Fig. 1

        Depends on beam solid angle (forward gain)

        #### Exposed parameters:

        u0_int_c150
        du0_int_c150
        u0_int_g150
        du0_int_g150
        u0_int_m150
        du0_int_m150
        # Similarly for 90 GHz

        '''

        fss_u0 =  ( 1./np.cos((np.radians(phip)+np.radians(alpha))/2.) +
            1./np.sin((np.radians(phip)-np.radians(alpha))/2.))/(2.*np.power(2.*np.pi, 0.5))

        fss_du0 = (1./np.cos((np.radians(phip)+np.radians(alpha))/2.) -
            1./np.sin((np.radians(phip)-np.radians(alpha))/2.))/(2.*np.power(2.*np.pi, 0.5))

        if verbose:
            print('Typical values of f(theta, phi):')
            print('    for alpha = {} deg, phi = {} deg (u=0, du=0 amp & power): {:.3f}, {:.3f}, {:.3f}, {:.3f}'.\
                format(alpha, phip, fss_u0, fss_du0, fss_u0**2, fss_du0**2))

        alp_s = np.radians(1. + np.arange(nang)) #in degrees
        theta_sv = 1. + np.arange(nang)
        fsints_u0 = np.zeros(nang) #integrand gets a sin^2
        fsints_du0 = np.zeros(nang) #integrand gets a sin^2

        zcv  = np.zeros(self.ngst)

        self.u0_int_c = np.zeros((self.nfreqs, self.ngst))
        self.du0_int_c = np.zeros((self.nfreqs, self.ngst))
        self.u0_int_g = np.zeros((self.nfreqs, self.ngst))
        self.du0_int_g = np.zeros((self.nfreqs, self.ngst))
        self.u0_int_m = np.zeros((self.nfreqs, self.ngst))
        self.du0_int_m = np.zeros((self.nfreqs, self.ngst))

        for i, freq in enumerate(self.freqs):

            tm_u0_c = np.zeros(self.ngst)
            tm_du0_c = np.zeros(self.ngst)
            tm_u0_g = np.zeros(self.ngst)
            tm_du0_g = np.zeros(self.ngst)
            tm_u0_m = np.zeros(self.ngst)
            tm_du0_m = np.zeros(self.ngst)

            for j in range(self.ngst):

                ### Now do the integral over alpha. In this equation \alpha is positive.
                fss =  np.sin(alp_s) * (1./np.cos((self.phi[j] + alp_s)/2.) + \
                    1./np.sin((self.phi[j] - alp_s)/2.)) / (2. * np.sqrt(2*np.pi))
                # fss =  ( 1./np.cos((self.phi*dtr+alp_s)/2.) +
                #   1./np.sin((self.phi*dtr-alp_s)/2.))/(2.*np.power(2.*np.pi,0.5))

                fsints_u0 = np.trapz(fss * fss, alp_s)
                # if (i == 0): print 'i=0 fsints_u0: ',fsints_u0

                fss = np.sin(alp_s) * (1./np.cos((self.phi[j] + alp_s)/2.) -\
                    1./np.sin((self.phi[j] - alp_s)/2.)) / (2.*np.sqrt(2.*np.pi))
                # fss = (1./np.cos((self.phi*dtr+alp_s)/2.) -
                #   1./np.sin((self.phi*dtr-alp_s)/2.))/(2.*np.power(2.*np.pi,0.5))
                fsints_du0 = np.trapz(fss * fss, alp_s)
                # if (i == 0): print 'i=0 fsints_du0: ',fsints_du0

                # print 'loc: ',gst, '| zc (cm): ','%.1f'%(zc), '| s (cm):','%.1f'%(sss),'| psi (deg):','%.1f'%(psi/dtr),\
                #   '| self.phi (deg):','%.1f'%(self.phi/dtr),'| beta (deg):','%.1f'%(beta/dtr)

                ### Now get the other factors

                gainfac_c = gain_model(freq)/np.sin(self.psi_c[j])**3
                gainfac_g = gain_model(freq)/np.sin(self.psi_g[j])**3
                gainfac_m = gain_model(freq)/np.sin(self.psi_m[j])**3
                gammafac = np.cos(self.rim_ang) * self.rim_fac
                # print 'gamma', gammafac

                tgeofac = 270*np.pi/(4.*kwave(freq))

                tm_u0_c[j] = fsints_u0*gainfac_c*gammafac[j]*tgeofac/self.sss[j]
                tm_du0_c[j] = fsints_du0*gainfac_c*gammafac[j]*tgeofac/self.sss[j]
                tm_u0_g[j] = fsints_u0*gainfac_g*gammafac[j]*tgeofac/self.sss[j]
                tm_du0_g[j] = fsints_du0*gainfac_g*gammafac[j]*tgeofac/self.sss[j]
                tm_u0_m[j] = fsints_u0*gainfac_m*gammafac[j]*tgeofac/self.sss[j]
                tm_du0_m[j] = fsints_du0*gainfac_m*gammafac[j]*tgeofac/self.sss[j]

            tm_u0_c[1:] = 2. * tm_u0_c[1:]
            tm_du0_c[1:] = 2. * tm_du0_c[1:]
            tm_u0_g[1:] = 2. * tm_u0_g[1:]
            tm_du0_g[1:] = 2. * tm_du0_g[1:]
            tm_u0_m[1:] = 2. * tm_u0_m[1:]
            tm_du0_m[1:] = 2. * tm_du0_m[1:]

            self.u0_int_c[i] = np.cumsum(tm_u0_c)
            self.du0_int_c[i] = np.cumsum(tm_du0_c)
            self.u0_int_g[i] = np.cumsum(tm_u0_g)
            self.du0_int_g[i] = np.cumsum(tm_du0_g)
            self.u0_int_m[i] = np.cumsum(tm_u0_m)
            self.du0_int_m[i] = np.cumsum(tm_du0_m)

        if verbose:
            for i, freq in enumerate(self.freqs):
                print('\nTotal ground pickup at {} GHz'.format(freq))
                print('Cyan E_hor E_vert (uK)', '%.1f'%(self.u0_int_c[i, -1]*1.e6),\
                    '%.1f'%(self.du0_int_c[i, -1]*1.e6))
                print('Green E_hor E_vert (uK)', '%.1f'%(self.u0_int_g[i, -1]*1.e6),\
                    '%.1f'%(self.du0_int_g[i, -1]*1.e6))
                print('Magenta E_hor E_vert (uK)', '%.1f'%(self.u0_int_m[i, -1]*1.e6),\
                    '%.1f'%(self.du0_int_m[i, -1]*1.e6))

    def draw_geometry_top(self, ymin=-4, ymax=4, xmin=0, xmax=8, fign=2):

        fig = plt.figure(num=fign)
        ax = plt.axes()

        plt.plot(self.xl, self.yl, 'b-', lw=1)
        plt.plot(self.x_gs_rim, self.y_gs_rim, 'k-', lw=2)

        rayxl = 6. #m

        itry = 0
        rays_1_x_t = np.array([self.t.x,self.t.x+rayxl])
        rays_1_y_t = np.array([self.t.r, self.t.r + rayxl*np.tan(self.ang_to_gs[itry]) ])
        plt.plot(rays_1_x_t, rays_1_y_t, color='grey', linewidth =.5 )
        rays_1_x_c = np.array([self.t.x,self.t.x+rayxl])
        rays_1_y_c = np.array([0, rayxl*np.tan(self.ang_to_gs[itry]) ])
        plt.plot(rays_1_x_c, rays_1_y_c, color='g', linewidth =.5 )
        rays_1_x_b = np.array([self.t.x,self.t.x+rayxl])
        rays_1_y_b = np.array([-self.t.r, -self.t.r + rayxl*np.tan(self.ang_to_gs[itry]) ])
        plt.plot(rays_1_x_b, rays_1_y_b, color='grey', linewidth =.5 )
        itry = 1
        rays_1_x_t = np.array([self.t.x,self.t.x+rayxl])
        rays_1_y_t = np.array([self.t.r, self.t.r + rayxl*np.tan(self.ang_to_gs[itry]) ])
        plt.plot(rays_1_x_t, rays_1_y_t, color='grey', linewidth =.5 )
        rays_1_x_c = np.array([self.t.x,self.t.x+rayxl])
        rays_1_y_c = np.array([0, rayxl*np.tan(self.ang_to_gs[itry]) ])
        plt.plot(rays_1_x_c, rays_1_y_c, color='g', linewidth =.5 )
        rays_1_x_b = np.array([self.t.x,self.t.x+rayxl])
        rays_1_y_b = np.array([-self.t.r, -self.t.r + rayxl*np.tan(self.ang_to_gs[itry]) ])
        plt.plot(rays_1_x_b, rays_1_y_b, color='grey', linewidth =.5 )
        itry = 8
        rays_1_x_t = np.array([self.t.x,self.t.x+rayxl])
        rays_1_y_t = np.array([self.t.r, self.t.r + rayxl*np.tan(self.ang_to_gs[itry]) ])
        plt.plot(rays_1_x_t, rays_1_y_t, color='grey', linewidth =.5 )
        rays_1_x_c = np.array([self.t.x,self.t.x+rayxl])
        rays_1_y_c = np.array([0, rayxl*np.tan(self.ang_to_gs[itry]) ])
        plt.plot(rays_1_x_c, rays_1_y_c, color='g', linewidth =.5 )
        rays_1_x_b = np.array([self.t.x,self.t.x+rayxl])
        rays_1_y_b = np.array([-self.t.r, -self.t.r + rayxl*np.tan(self.ang_to_gs[itry]) ])
        plt.plot(rays_1_x_b, rays_1_y_b, color='grey', linewidth =.5 )

        plt.xlabel('x (m)',fontsize=14)
        plt.ylabel('y (m)', fontsize=14)
        plt.ylim(ymin, ymax)
        plt.xlim(xmin, xmax)
        plt.savefig('img/gtd_3d_top_geom_sac.png')

    def draw_geometry_side(self, ymin=2, ymax=4.75, xmin=0, xmax=5.,
        units='m', fign=1):

        fig = plt.figure(fign)
        ax = plt.gca()

        ax.arrow(self.x_gs[1], self.y_gs[1], self.x_gs_ap[1] - self.x_gs_ap[0],
            self.y_gs_ap[1] - self.y_gs_ap[0], length_includes_head=True,
            head_width=0.08, head_length=0.16, fc='k', ec='k')

        ax.arrow(self.x_inc[0], self.y_inc[0], self.x_inc[1] - self.x_inc[0],
            self.y_inc[1] - self.y_inc[0], length_includes_head=True,
            head_width=0.08, head_length=0.16, fc='k', ec='k')

        # Aperture
        plt.plot(self.x_ap, self.y_ap, 'b-', lw=1)
        plt.plot(self.x_oa, self.y_oa, 'b--', lw=1)
        plt.plot(self.x_oa_t, self.y_oa_t, 'b--', lw=1)
        plt.plot(self.x_oa_b, self.y_oa_b, 'b--', lw=1)

        # Ground screen and angle arrows
        plt.plot(self.x_gs, self.y_gs, 'k-', lw=2)

        plt.xlabel('x (m)', fontsize=14)
        plt.ylabel('z (m)', fontsize=14)

        plt.ylim(ymin, ymax)
        plt.xlim(xmin, xmax)

        plt.savefig('img/gtd_3d_side_geom_sac.png', bbox_inches='tight')

    def plot_cum(self, fignum=3, fs=12):

        #
        # Aperture
        #

        for i, freq in enumerate(self.freqs):

            fig = plt.figure(num=fignum, figsize=(6,6))
            ax = plt.axes()
            ymin = 0
            ymax = 400 #cm
            xmin = 0
            xmax = 90 #cm

            plt.plot(np.degrees(self.ang_to_gs), 1.e6*self.u0_int_c[i], 'c--', lw=2)
            plt.plot(np.degrees(self.ang_to_gs), 1.e6*self.du0_int_c[i], 'c-', lw=4)
            plt.plot(np.degrees(self.ang_to_gs), 1.e6*self.u0_int_g[i], 'g--', lw=2)
            plt.plot(np.degrees(self.ang_to_gs), 1.e6*self.du0_int_g[i], 'g-', lw=4)
            plt.plot(np.degrees(self.ang_to_gs), 1.e6*self.u0_int_m[i], 'm--', lw=2)
            plt.plot(np.degrees(self.ang_to_gs), 1.e6*self.du0_int_m[i]/10, 'm-', lw=4)
            # plt.plot(x_oa_t, y_oa_t, 'b--', lw = 1 )
            # plt.plot(x_oa_b, y_oa_b, 'b--', lw = 1 )
            plt.text(40,248,'signal/10 ',color='m',fontsize=14)
            plt.text(10,375,'Vertical polarization, du/dn=0 (thick)',color='k', fontsize=fs)
            plt.text(10,355,'Horizontal polarization, u=0 (thin dash)',color='k', fontsize=fs)
            plt.text(10,335,'Magenta (bottom) ',color='m', fontsize=fs)
            plt.text(10,315,'Green (middle) ',color='g', fontsize=fs)
            plt.text(10,295,'Cyan (top) ',color='c', fontsize=fs)
            plt.title('Cum contribution at {} GHz, 42 cm ap diam'.format(freq),
                fontsize=14)
            plt.xlabel('Angle along top of ground screen from ap (deg)', fontsize=14)
            plt.ylabel(r'Measured temp for $T_g=270$K ($\mu$K)',  fontsize=14)
            plt.ylim(ymin, ymax)
            plt.xlim(xmin, xmax)
            plt.savefig('img/cum_pickup_sac_{}.png'.format(freq))
            plt.close()


        print(' ')
        gainfacs = np.zeros([self.ngst, 3, self.npatch]) # num of positions arounf rim, 3 beams, n patches
        gainfacs = 5.14e-4/np.power(np.sin(self.gs_psi), 3)
        temit = 0.35 #K
        for i in range(3):
            ssum = 0
            sum_sa = 0
            for j in range(self.ngst):
                for k in range(self.npatch):
                    sum_sa = sum_sa + self.gs_sa[j,k]
                    ssum = ssum + self.gs_sa[j,k]*gainfacs[j,i,k]*temit/(4.*np.pi)
            print('For i the sum is (K) ', i, ssum, sum_sa,sum_sa/(4.*np.pi))

class Beam(object):

    def __init__():

        pass

def compare_with_lp_results():
    '''
    A test to see if the code is still producing results that are consistent with
    LP's original calculations. This should be the case until we identify
    some error in the original calculation
    '''
    pass

def basic():

    import seaborn as sns

    font = {'family' : 'normal',
        'size'   : 16}

    matplotlib.rc('font', **font)

    # sns.set_style("whitegrid")
    # sns.set_style("ticks")


    freqs = [90, 150]
    elevations = np.array([30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90,
        100, 110, 120, 130, 140, 150, 160, 170, 180])


    u0_c = np.nan * np.ones((len(freqs), len(elevations)))
    du0_c = np.nan * np.ones((len(freqs), len(elevations)))

    u0_g = np.nan * np.ones((len(freqs), len(elevations)))
    du0_g = np.nan * np.ones((len(freqs), len(elevations)))

    u0_m = np.nan * np.ones((len(freqs), len(elevations)))
    du0_m = np.nan * np.ones((len(freqs), len(elevations)))

    for e, elevation in enumerate(elevations):


        t = Telescope(elevation=elevation)
        bg = BaffleGeometry(t=t, freqs=freqs)

        t1 = time.time()

        # Geometrical calculations
        bg.irimwalk()
        bg.iuconst()

        # Loading calculation
        bg.pickup()

        for i, freq in enumerate(bg.freqs):
            u0_c[i, e] = bg.u0_int_c[i, -1]
            u0_g[i, e] = bg.u0_int_g[i, -1]
            u0_m[i, e] = bg.u0_int_m[i, -1]


        t2 = time.time()
        print('Time to calculate geometry: {} ms'.format(1e3*(t2-t1)))


    for i, freq in enumerate(bg.freqs):
        plt.semilogy(elevations, 1e6*u0_c[i, :], label='{} GHz (top)'.format(freq))
        plt.semilogy(elevations, 1e6*u0_g[i, :], label='{} GHz (middle)'.format(freq))
        plt.semilogy(elevations, 1e6*u0_m[i, :], label='{} GHz (bottom)'.format(freq))


    plt.legend(ncol=2)
    plt.ylabel('Signal amplitude [uK]')
    plt.xlabel('Elevation [deg]')
    plt.savefig(opj('img/', 'elevations.png'), dpi=300)
    plt.close()

def main():

    def_outdir = opj(os.getcwd())

    parser = ap.ArgumentParser(formatter_class=\
        ap.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--outdir', default=def_outdir, required=False,
                        help="Output directory", type=str)


    parser.add_argument('--dry', action='store_true', dest='dry',
                        default=False)
    parser.add_argument('--test', action='store_true', dest='test',
                        default=True)

    if test:

        compare_with_lp_results()
        return


    basic_calculation()

if __name__ == '__main__':

    main()

    # bg.draw_geometry_side()
    # bg.draw_geometry_top()
    # bg.plot_cum()
