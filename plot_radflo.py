from pylab import *
import numpy as np
import imageio
import matplotlib.animation as animation
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['figure.max_open_warning'] = 0

def set_globals():
   global sielabl, tevlabl, gamlabl, uclabl
   global hnulabl, tklabl, uklabl
   global nrho, nkt, nhv, nhvp
   
   gamlabl = 'Gamma -1 '
   sielabl = 'SIE (erg/gm)'
   tklabl  = 'T (K)'
   tevlabl = 'T (eV)'
   uklabl  = 'UK (cm**2/gm)'
   hnulabl = 'HNU (eV)'
   uclabl  = 'UC (m/s)'
   
   nrho = 7
   nkt  = 100
   nhv  = 51
   nhvp = nhv + 1

# HYCHEM2013AF read ptdata
def rdptdata(fname):
   mydmp = open(fname) 
   mydmp.readline()
   pl_title = mydmp.readline()
   mydmp.readline()

   line   = mydmp.readline()
   values = line.split()
   npts   = int(values[0]) 

   for i in range(0, 2, 1): mydmp.readline()

   tpwr     = np.zeros([npts])
   pow_band = np.zeros([npts,16])
   totpow   = np.zeros([npts])
   pblue    = np.zeros([npts])

   for iline in range(0,npts,1):
      line   = mydmp.readline()
      values = [float(s) for s in line.split()]
      tpwr[iline]           = values[0]
      pow_band[iline,0:16]  = values[1:17]
      pblue[iline]          = values[12]
      totpow[iline]   = 0. 
      for k in range(1,17,1):
         totpow[iline] = totpow[iline] + values[k]
   mydmp.close()
   return pl_title, npts, tpwr, pow_band, totpow, pblue
   
def rdpowpts(fname):
   mydmp = open(fname)
   for i in range(0,3,1): mydmp.readline()
   line = mydmp.readline()
   values = line.split()
   npts   = int(values[0])
   mydmp.close()
   return npts 
   
def rdpower(npt, fname): 
   tpwr     = np.zeros([npt])
   pow_band = np.zeros([npt,16])
   totpow   = np.zeros([npt])
   sipow    = np.zeros([npt])
   wgt = [0., 0., 0., 0., 0., 0., 0.164, 0.943, 0.846, 0.629, \
          0.384, 0.117, 0., 0., 0., 0.]
   mydmp = open(fname)
   mydmp.readline()
   pl_title = mydmp.readline()
   for i in range(0,4,1): mydmp.readline()

#  wgt index starts at 0

   for iline in range(0,npt,1):
      line   = mydmp.readline()
      values = [float(s) for s in line.split()]
      tpwr[iline]           = values[0]
      pow_band[iline,0:15]  = values[1:16]
      totpow[iline]   = 0.
      sipow[iline]    = 0. 
      for k in range(0,16,1):
         totpow[iline] = totpow[iline] + pow_band[iline,k]
         sipow[iline]  = sipow[iline]  + pow_band[iline,k] * wgt[k]
   mydmp.close()
   return pl_title, tpwr, pow_band, totpow
   
def plotpwr(tpwr, pow_band, totpow, npts, pl_title):

   xmax  = 1.1 * max(tpwr)
   xmin  = 5.e-7

   ymin = 1.e30
   ymax = 0.
   for i in range(0,npts):
     if(tpwr[i] > xmin):
       if(totpow[i] > ymax): ymax = totpow[i]
       if(totpow[i] < ymin): ymin = totpow[i]

   yminchk = 1.e-5 * ymax
   if(ymin < yminchk): ymin = yminchk
   ymax = 1.2 * ymax
   ymin = 0.8 * ymin
   
   thyld = 0.
   for i in range(0,npts-1):
     thyld = thyld + (tpwr[i+1]-tpwr[i])*(totpow[i+1]+totpow[i])/2.
   thyld = thyld / 4.185e12
   tyld_labl = ("Thermal Yield (kt) = %.3f" % (thyld))   
   hvtitl    = pl_title
   
   fig = plt.figure(figsize=(9.0, 6.5))
   plt.title(hvtitl)
   plt.xscale('log')
   plt.yscale('log')
   plt.xlim(xmin, xmax)
   plt.ylim(ymin, ymax)
   plt.xlabel('T (s)')
   plt.ylabel('Thermal Power (W)')
   plt.plot(tpwr, totpow, color='blue')
   plt.grid(True)
   plt.gcf().text(.35, .89, tyld_labl, fontsize=12, color='blue')
   pdf_pages.savefig(fig)
 
def plotpwr2(tpwr, sipow, npts, tpwr2, sipow2, npts2, lably, hvtitl):
   xmax1  = max(tpwr)
   xmax2  = max(tpwr2)
   xmax   = 1.1 * max(xmax1, xmax2)
   xmin   = 1.e-6

   ymax = 0.
   ymin = 1.e20
   for i in range(0,npts,1):
     if (tpwr[i] > xmin):
       if(sipow[i] > ymax): ymax = sipow[i]
       if(sipow[i] < ymin): ymin = sipow[i]

   for i in range(0,npts2,1):
     if (tpwr2[i] > xmin):
       if(sipow2[i] > ymax): ymax = sipow2[i]
       if(sipow2[i] < ymin): ymin = sipow2[i]
   ymax = 1.5 * ymax
   ymin = 0.9 * ymin
   
   fig = plt.figure(figsize=(9.0, 6.5))
   plt.title(hvtitl)
   plt.xscale('log')
   plt.yscale('log')
   plt.xlim(xmin, xmax)
   plt.ylim(ymin, ymax)
   plt.xlabel('T (s)')
   plt.ylabel(lably)
   plt.plot(tpwr, sipow, color='red')
   plt.plot(tpwr2, sipow2,  color='blue')
   plt.grid(True)
   pdf_pages.savefig(fig)  
    
def rdrtvar(tpwr, ptherm, ke, ie, erad, etot, rshk, ritmx, titmxm, 
        titmx, titmxp, rhorat, flupow, npts, fname):
   mydmp = open(fname)
   for i in range(0,5,1): mydmp.readline()
   for iline in range(0,npts,1):
      line   = mydmp.readline()
      values = [float(s) for s in line.split()]
      tpwr[iline]   = values[0]
      ptherm[iline] = values[1]
      ke[iline]     = values[2]
      ie[iline]     = values[3]
      erad[iline]   = values[4] 
      etot[iline]   = values[5]
      rshk[iline]   = values[6]
      ritmx[iline]  = values[7]
      titmxm[iline] = values[8]
      titmx[iline]  = values[9]
      titmxp[iline] = values[10]
      rhorat[iline] = values[11] 
      flupow[iline] = values[12] 
   mydmp.close()

def plotrthyd(tpwr, ptherm, ke, ie, erad, etot, rshk, ritmx, titmxm, 
        titmx, titmxp, rhorat, npts):

   ismooth = int(25) 
   nsmooth = int((npts - ismooth)) 
    
   vshk    = np.zeros((nsmooth))
   tshk    = np.zeros((nsmooth))
   vfb     = np.zeros((nsmooth))
   for i in range(0,nsmooth,1):
     tshk[i] = 0.5 * (tpwr[i] + tpwr[i+ismooth])
     vshk[i] = 1.e-5 * ( rshk[i+ismooth]  - rshk[i])/(tpwr[i+ismooth] - tpwr[i])
     vfb[i]  = 1.e-5 * (ritmx[i+ismooth] - ritmx[i])/(tpwr[i+ismooth] - tpwr[i])

   hvtitl = 'Total (blue), Internal (green), Kinetic (red)'

   ymax  = 1.05 * max(etot)
   ymin  = 0.9 * min(ke)
   ymin2 = max(ymin, 1.e-6*ymax)
   xmax  = 1.1 * max(tpwr)
   xmin  = tpwr[0]   
   fig   = plt.figure(figsize=(9.0, 6.5))
   plt.title(hvtitl)
   plt.xscale('log')
   plt.yscale('linear')
   plt.xlim( xmin, xmax)
   plt.ylim(0., ymax)
   plt.xlabel('T (s)')
   plt.ylabel('Energy (kt)')
   plt.plot(tpwr, etot, color ='blue')
   plt.plot(tpwr, ie,   color = 'green')
   plt.plot(tpwr, ke,   color = 'red')
   plt.grid(True)
   pdf_pages.savefig(fig)

   hvtitl = 'Thermal Power With No Averaging'

   ymax  = 1.5*ptherm.max()
   ymin  = 0.8*ptherm.min()
   if(ymin < 1.e-5*ymax): ymin = 1.e-5*ymax
   fig   = plt.figure(figsize=(9.0, 6.5))
   plt.title(hvtitl)
   plt.xscale('log')
   plt.yscale('log')
   plt.ylim(ymin, ymax)
   plt.xlabel('T (s)')
   plt.ylabel('Power (W)')
   plt.plot(tpwr, ptherm,  color='blue')
   plt.grid(True)
   pdf_pages.savefig(fig)
   
   hvtitl = 'R Shock in blue,  R(itmx) in red'

   ritmx = 0.01 * ritmx		# convert to m
   rshk  = 0.01 * rshk		# convert to m   
   ymax  = 1.05 * max(rshk)
   ymin  = 0.0
   fig   = plt.figure(figsize=(9.0, 6.5))
   plt.title(hvtitl)
   plt.xscale('log')
   plt.yscale('log')
   plt.ylim(0.10, ymax)
   plt.xlabel('T (s)')
   plt.ylabel('Radii (m)')
   plt.plot(tpwr, rshk,  color='blue')
   plt.plot(tpwr, ritmx, color='red')
   plt.grid(True)
   pdf_pages.savefig(fig)
   
   hvtitl = "V Shock in blue,  V(itmx) in red,   smooth points = %3.0f" %ismooth
   
   xmin    = 1.e-6
   xmax    = 1.1 * max(tshk) 
   vshkmax = 0.0
   vfbmax  = 0.0
   for i in range(0,nsmooth,1):
     if(tshk[i] > xmin):
       if (vshk[i] > vshkmax): vshkmax = vshk[i]
       if (vfb[i]  >  vfbmax): vfbmax  = vfb[i]
   ymax  = 1.05 * vshkmax
   ymax2 = 1.05 * vfbmax
   if(ymax2 > ymax): ymax = ymax2
  
   fig   = plt.figure(figsize=(9.0, 6.5))
   plt.title(hvtitl)
   plt.xscale('log')
   plt.yscale('log')
   plt.xlim(xmin, xmax)
   plt.ylim(0.2, ymax)
   plt.xlabel('T (s)')
   plt.ylabel('Velocity (km/s)')
   plt.plot(tshk, vshk,  color='blue')
   plt.plot(tshk, vfb,   color='red')
   plt.grid(True)
   pdf_pages.savefig(fig)

   fig   = plt.figure(figsize=(9.0, 6.5))
   plt.title(hvtitl)
   plt.xscale('log')
   plt.yscale('linear')
   plt.xlim(xmin, xmax)
   plt.ylim(0.01, ymax)
   plt.xlabel('T (s)')
   plt.ylabel('Velocity (km/s)')
   plt.plot(tshk, vshk,  color='blue')
   plt.plot(tshk, vfb,   color='red')
   plt.grid(True)
   pdf_pages.savefig(fig)

   hvtitl = 'T(itmx) in blue, T(itmx+1) in green, T(itmx+2) in red '

   ymax  = 1.5 * max(titmxm)
   ymin  = 0.9 * min(titmxp)
   ymin2 = max(ymin, 1.e-6*ymax)   
   fig   = plt.figure(figsize=(9.0, 6.5))
   plt.title(hvtitl)
   plt.xscale('log')
   plt.yscale('log')
   plt.xlim( xmin, xmax)
   plt.ylim(ymin2, ymax)
   plt.xlabel('t (s)')
   plt.ylabel('T (eV)')
   plt.plot(tpwr, titmxm, color='blue')
   plt.plot(tpwr, titmx,  color='green')
   plt.plot(tpwr, titmxp, color='red')
   plt.grid(True)
   pdf_pages.savefig(fig) 

   hvtitl = 'Max RHO ratio in Air '

   ymax  = 1.1 * max(rhorat)
   ymin  = 0.9    
   fig   = plt.figure(figsize=(9.0, 6.5))
   plt.title(hvtitl)
   plt.xscale('log')
   plt.yscale('linear')
   plt.xlim(xmin, xmax)
   plt.ylim(ymin, ymax)
   plt.xlabel('T (s)')
   plt.ylabel('RHO Ratio')
   plt.plot(tpwr, rhorat, color='blue')
   plt.grid(True)
   pdf_pages.savefig(fig)     
   
def rd_header():
   label = mydmp.readline()
   mydmp.readline()
   time   = float(mydmp.readline())
   ncycle = int(mydmp.readline())
   ndc    = int(mydmp.readline())   
   itmax  = int(mydmp.readline())
   irmax  = int(mydmp.readline())
   imax   = int(mydmp.readline())
   mydmp.readline()   
   mydmp.readline()
   
   return imax, label, time, ncycle 

def rdxdep(fname):
   mydmp = open(fname)

#  xdep.owt is written twice - coarse grid and full grid

#  coarse grid

   mydmp.readline()
   pl_title = mydmp.readline()
   for i in range(3,13): mydmp.readline()

   time   = float(mydmp.readline())
   xyld   = float(mydmp.readline())
   radyld = float(mydmp.readline())

   for i in range(0,7,1): mydmp.readline()

   npts   = int(mydmp.readline())
   rc     = np.zeros((npts))
   rho    = np.zeros((npts))
   sie    = np.zeros((npts))
   uc     = np.zeros((npts))
   temp   = np.zeros((npts))
   de     = np.zeros((npts))
   mass   = np.zeros((npts))

   mydmp.readline()
   mydmp.readline()
   for iline in range(0,npts,1):
      line   = mydmp.readline()
      values = [float(s) for s in line.split()]
      rc[iline]   = values[1]
      rho[iline]  = values[2]
      sie[iline]  = values[3]
      uc[iline]   = values[4]
      temp[iline] = values[5] 
      mass[iline] = values[6]

#  full grid

   mydmp.readline()
   pl_title = mydmp.readline()
   for i in range(3,13): mydmp.readline()

   time2   = float(mydmp.readline())
   xyld2   = float(mydmp.readline())
   radyld2 = float(mydmp.readline())

   for i in range(0,7,1): mydmp.readline()

   npts2   = int(mydmp.readline())
   rc2     = np.zeros((npts2))
   rho2    = np.zeros((npts2))
   sie2    = np.zeros((npts2))
   uc2     = np.zeros((npts2))
   temp2   = np.zeros((npts2))
   de2     = np.zeros((npts2))
   mass2   = np.zeros((npts2))

   mydmp.readline()
   mydmp.readline()
   for iline in range(0,npts2,1):
      line   = mydmp.readline()
      values = [float(s) for s in line.split()]
      rc2[iline]   = values[1]
      rho2[iline]  = values[2]
      sie2[iline]  = values[3]
      uc2[iline]   = values[4]
      temp2[iline] = values[5] 
      mass2[iline] = values[6]


   mydmp.close()
   return pl_title, rc, temp, npts, rc2, temp2, npts2

def plotxdep(rc, temp, npts, rc2, temp2, npts2, pl_title):
 
   hvtitl    = pl_title 
   subtitl = 'Coarse Grid (blue) and Full Grid (red) at end of X-ray Deposition'

   for i in range(npts2-1,0,-1):
      if(temp2[i] > 9000.):
         r9000 = rc2[i]
         break
   
   fig = plt.figure(figsize=(9.0, 6.5))
   plt.title(hvtitl)
   plt.xscale('log')
   plt.yscale('log')
   plt.xlabel('RC (cm)')
   plt.ylabel('Temp (K)')
   plt.plot(rc,  temp,  color='blue')
   plt.plot(rc2, temp2, color='red')
   xpl = [r9000, r9000]
   ypl = [300., 10000.]
   plt.plot(xpl, ypl, color='black', linestyle='dashed')
   plt.grid(True)
   plt.gcf().text(.25, .89, subtitl, fontsize=11, color='black')
   pdf_pages.savefig(fig)
   
def rd_rflo_dmp(fname):

   mydmp  = open(fname,"r")  
   runlbl = mydmp.readline()
   runlbl = runlbl.strip()
   mydmp.readline()
   time   = float(mydmp.readline())
   ncycle = int(mydmp.readline())
   ndc    = int(mydmp.readline())
   itmx   = int(mydmp.readline())
   irmx   = int(mydmp.readline())
   icmax  = int(mydmp.readline())

   rc    = np.zeros((icmax))  
   dr    = np.zeros((icmax))
   rho   = np.zeros((icmax))
   uc    = np.zeros((icmax))
   pr    = np.zeros((icmax))
   tkel  = np.zeros((icmax))
   mass  = np.zeros((icmax))
   sie   = np.zeros((icmax))
   cs    = np.zeros((icmax))   
   gamma = np.zeros((icmax))

   mydmp.readline()
   
   for iline in range(0,icmax,1):
      line = mydmp.readline()
      values = [float(s) for s in line.split()]
      rc[iline]    = values[0]
      dr[iline]    = values[1]
      rho[iline]   = values[2]
      uc[iline]    = values[3]
      pr[iline]    = values[4]
      tkel[iline]  = values[5]
      mass[iline]  = values[6]
      sie[iline]   = values[7] 
      cs[iline]    = values[8]     
      gamma[iline] = values[9]

   mydmp.close()
   return icmax, rc, dr, rho, uc, pr, tkel, mass, sie, cs, gamma, time, ncycle, runlbl
              
def rdplanck(fname):
   mydmp = open(fname)
   mydmp.readline()
   
# HNUR(nhvp) table
#  The data file is written with 7 values per line for etable

   mydmp.readline()
   ibeg = 0
   for iline in range(0,8,1):
     iend = min(ibeg+7, nhvp)
     line = mydmp.readline()
     values = [float(s) for s in line.split()]
     hnur[ibeg:iend] = values[0:7]
     ibeg += 7
     
#   print(hnur) 
   
# KT(nkt) table
#  The data file is written with 7 values per line for etable

   mydmp.readline()
   ibeg = 0
   for iline in range(0,15,1):
     iend = min(ibeg+7, nkt+1)
     line = mydmp.readline()
     values = [float(s) for s in line.split()]
     kt[ibeg:iend] = values[0:7]
     ibeg += 7
     
#   print(kt) 

# B(nkt, nhv) table 
# there are eight values per line

   for mval in range(0,nhv):
      mydmp.readline()
      ibeg = 0
      for iline in range(0,13,1): 
         iend = min(ibeg+8, nkt+1)
         line = mydmp.readline()
         values = [float(s) for s in line.split()]
         b[ibeg:iend,mval] = values[0:8]
         ibeg += 8 
         
   mydmp.close()
   return b, hnur, kt  
   
def plot_b(b, hnur, kt, hmid):

   colorvar = cm.rainbow(np.linspace(0,1,100))

   ymax  = b.max()
   ymin  = b.min()
   if(ymin < 10.): ymin = 10.
   
   hvtitl   = "pi * Integral[B(kt,hv)]dhv / delta(hv) for various Temperatures"
   
   fig = plt.figure(figsize=(9.0, 6.5))
   plt.title(hvtitl)
   plt.xscale('log')
   plt.yscale('log')
   plt.ylim(ymin, ymax)
   plt.xlabel('Photon Energy (eV)')
   plt.ylabel('B Function / delta(hv)')

   for j in range(3,100,16):
      yval = b[j,:]
      if (j == 3): plt.plot(hmid, yval, color=colorvar[j], label=str(kt[j]) + "  eV ")
      if (j > 3):  plt.plot(hmid, yval, color=colorvar[j], label=kt[j])

   plt.legend(loc='center left',bbox_to_anchor=(0.05,0.82), fontsize="small")

#   plt.grid(True)                 
   pdf_pages.savefig(fig)
         
def rdeos(fname):
   mydmp = open(fname,"r")
   mydmp.readline()
   nr  = int(mydmp.readline())
   nk  = int(mydmp.readline())
   
#  RHOTABL

   mydmp.readline()
   line = mydmp.readline()
   rhotabl = [float(s) for s in line.split()]
   
#   print rhotabl
   
#  ETABLE
#  The data file is written with 7 values per line for etable
#  nkt = etable size

   mydmp.readline()
   ibeg = 0
   for iline in range(0,15,1):
     iend = min(ibeg+7, nkt+1)
     line = mydmp.readline()
     values = [float(s) for s in line.split()]
     etable[ibeg:iend] = values[0:7]
     ibeg += 7
     
#   print(etable)
  
# GT array
#  The data file is written with 7 values per line for etable
#  gt(nkt,nrho)

   mydmp.readline()
   for iline in range(0,nkt,1):
      line = mydmp.readline()
      values = [float(s) for s in line.split()]
      gt[iline][0:7] = values[0:7]
      
#   print(gt[0][0], gt[99][6])
       
# FP array
#  The data file is written with 7 values per line for etable
#  fp(nkt,nrho)

   mydmp.readline()
   for iline in range(0,nkt,1):
      line = mydmp.readline()
      values = [float(s) for s in line.split()]
      fp[iline][0:7] = values[0:7]
      
#   print(fp[0][0], fp[99][6]) 
   
# T(ev) array
#  The data file is written with 7 values per line for etable
#  tev(nkt,nrho)

   mydmp.readline()
   for iline in range(0,nkt,1):
      line = mydmp.readline()
      values = [float(s) for s in line.split()]
      tev[iline][0:7] = values[0:7]
      
#   print(tev[0][0], tev[99][6])      
   
   mydmp.close()
   return rhotabl, etable, fp, gt, tev
   
def ploteos(etable, rhotabl, fp, tev):

   colorvar = cm.rainbow(np.linspace(0,1,7))
   
# plot gamma-1 vs sie for various densities

   ymax  = 0.7
   ymin  = 0.0
   
   fig = plt.figure(figsize=(9.0, 6.5))
   plt.xscale('log')
   plt.yscale('linear')
   plt.ylim(ymin, ymax)
   plt.xlabel('SIE (erg/gm)')
   plt.ylabel('Gamma - 1')

   for j in range(0,7,1):
      yval = fp[:,j]
      if (j == 0): plt.plot(etable, yval, color=colorvar[j], label=str(rhotabl[j]) + " gm/cc")
      if (j > 0):  plt.plot(etable, yval, color=colorvar[j], label=rhotabl[j])

   plt.legend(loc='center left',bbox_to_anchor=(0.80,0.40), fontsize="small")

   plt.grid(True)
   pdf_pages.savefig(fig)

# plot (gamma+1)/(gamma-1) vs sie for various densities

   ymax  = 25.0
   ymin  = 1.0
   
   fig = plt.figure(figsize=(9.0, 6.5))
   plt.xscale('log')
   plt.yscale('linear')
   plt.ylim(ymin, ymax)
   plt.xlabel('SIE (erg/gm)')
   plt.ylabel('(Gamma+1) / (Gamma-1)')

   for j in range(0,7,1):
      yval = (fp[:,j] + 2.)/fp[:,j]
      if (j == 0): plt.plot(etable, yval, color=colorvar[j], label=str(rhotabl[j]) + " gm/cc")
      if (j > 0):  plt.plot(etable, yval, color=colorvar[j], label=rhotabl[j])

   plt.legend(loc='center left',bbox_to_anchor=(0.80,0.40), fontsize="small")

   plt.grid(True)
   pdf_pages.savefig(fig)
   
# Plot Tev versus Sie for various densities

   ymax  = 1.1 * tev.max()
   ymin  = 0.9 * tev.min()
   
   fig = plt.figure(figsize=(9.0, 6.5))
   plt.xscale('log')
   plt.yscale('log')
   plt.ylim(ymin, ymax)
   plt.xlabel('SIE (erg/gm)')
   plt.ylabel('T (eV)')

   for j in range(0,7,1):
      yval = tev[:,j]
      if (j == 0): plt.plot(etable, yval, color=colorvar[j], label=str(rhotabl[j]) + " gm/cc")
      if (j > 0):  plt.plot(etable, yval, color=colorvar[j], label=rhotabl[j])

   plt.legend(loc='center left',bbox_to_anchor=(0.80,0.40), fontsize="small")

   plt.grid(True)
   pdf_pages.savefig(fig) 
 
# Plot Sound Speed versus Sie for various densities

   ymax  = 10000000.
   ymin  = 100.
   
   fig = plt.figure(figsize=(9.0, 6.5))
   plt.xscale('log')
   plt.yscale('log')
   plt.ylim(ymin, ymax)
   plt.xlabel('SIE (erg/gm)')
   plt.ylabel('Sound Speed (m/s)')

   for j in range(0,7,1):
      yval = 0.01 * np.sqrt(fp[:,j] * (1. + fp[:,j]) * etable)
      if (j == 0): plt.plot(etable, yval, color=colorvar[j], label=str(rhotabl[j]) + " gm/cc")
      if (j > 0):  plt.plot(etable, yval, color=colorvar[j], label=rhotabl[j])

   plt.legend(loc='center left',bbox_to_anchor=(0.80,0.40), fontsize="small")

   plt.grid(True)
   pdf_pages.savefig(fig) 
     
def rdopac(fname, array_size, array_size2):

   my_uk   = np.zeros(array_size)
   my_hnur = np.zeros(array_size2)
   mydmp = open(fname,"r")
   for i in range(0,4,1): mydmp.readline()
   
# RHOTABL(nrho)

   mydmp.readline()
   line = mydmp.readline()
   rhotabl = [float(s) for s in line.split()]
   
#   print(rhotabl)
   
# KT(nkt) table
#  The data file is written with 7 values per line for etable

   mydmp.readline()
   ibeg = 0
   for iline in range(0,15,1):
     iend = min(ibeg+7, nkt+1)
     line = mydmp.readline()
     values = [float(s) for s in line.split()]
     kt[ibeg:iend] = values[0:7]
     ibeg += 7
     
#   print(kt)      
   
# HNUR(nhvp) table
#  The data file is written with 7 values per line for etable

   mydmp.readline()
   ibeg = 0
   for iline in range(0,8,1):
     iend = min(ibeg+7, nhvp)
     line = mydmp.readline()
     values = [float(s) for s in line.split()]
     my_hnur[ibeg:iend] = values[0:7]
     ibeg += 7
     
#   print(hnur) 
   
# UK(nkt,nrho,nhv)

   mydmp.readline()
   
   for m in range(0 ,51, 1):
     mydmp.readline()   
     
     for iline in range(0,100,1):
        line = mydmp.readline()
        values = [float(s) for s in line.split()]
        my_uk[iline][0][m] = values[0]
        my_uk[iline][1][m] = values[1] 
        my_uk[iline][2][m] = values[2]
        my_uk[iline][3][m] = values[3]  
        my_uk[iline][4][m] = values[4] 
        my_uk[iline][5][m] = values[5]
        my_uk[iline][6][m] = values[6]  
        
#   print(uk[0][0][0], uk[nkt-1][nrho-1][nhv-1])
   
   mydmp.close()
   return rhotabl, kt, my_hnur, my_uk

def plot_uk(uk, kt, hmid, rhotabl, tevlabl, uklabl, nhv, nrho, nkt, opactitl):

# Plot uk(nkt,nrho,nhv) versus hv  

   colorvar = cm.rainbow(np.linspace(0,1,7))
   dum = np.zeros((nrho,nhv))

   for i in range(1,nkt,5):
      for j in range(0,7,1): dum[j,:] = uk[i,j,:]                  
      
      ymax  = 2.0 * dum.max()
      ymin2 = 0.8 * dum.min()
      ymin  = max(ymin2, 1.e-12*ymax)
   
      hvtitl   = opactitl + "  T(eV)= %9.3f" %(kt[i])
   
      fig = plt.figure(figsize=(9.0, 6.5))
      plt.title(hvtitl)
      plt.xscale('log')
      plt.yscale('log')
      plt.ylim(ymin, ymax)
      plt.xlabel('Photon Energy (eV)')
      plt.ylabel('Rosseland Mean Opacity (cm^2/gm)')

      for j in range(0,7,1):
         yval = dum[j,:]
         if (j == 0): plt.plot(hmid, yval, color=colorvar[j], label=str(rhotabl[j]) + " gm/cc")
         if (j > 0):  plt.plot(hmid, yval, color=colorvar[j], label=rhotabl[j])

      plt.legend(loc='center left',bbox_to_anchor=(0.80,0.80), fontsize="small")

      plt.grid(True)
      pdf_pages.savefig(fig)

def plot_uk_compare(uk, hmid, uk2, hmid2, uk3, hmid3, kt, rhotabl, nhv, nrho, nkt):

# Plot uk(nkt,nrho,nhv) versus hv  

   colorvar = cm.rainbow(np.linspace(0,1,7))
   dum = np.zeros((nrho,nhv))

   for i in range(89,nkt,1):
      for j in range(2,7,3):                 

         yval  = uk[i,j,:] 
         yval2 = uk2[i,j,:] 
         yval3 = uk3[i,j,:]  

#         print(i,j, uk[i,j,42]/uk2[i,j,42], uk2[i,j,42]/uk3[i,j,42])
  
         ymax  = 2.0 * yval.max()
         ymin2 = 0.8 * yval.min()
         ymin  = max(ymin2, 1.e-12*ymax)
   
         hvtitl   = "  T(eV)= %9.3f  RHO (gm/cc)= %9.7f" %(kt[i], rhotabl[j])
   
         fig = plt.figure(figsize=(9.0, 6.5))
         plt.title(hvtitl)
         plt.xscale('log')
         plt.yscale('log')
         plt.ylim(ymin, ymax)
         plt.xlabel('Photon Energy (eV)')
         plt.ylabel('Rosseland Mean Opacity (cm^2/gm)')

         plt.plot(hmid,  yval,  color='red', label="AESOP51")
         plt.plot(hmid2, yval2, color='blue', label="AESOP51 REV1")
         plt.plot(hmid3, yval3, color='green', label="MUDATHA51")

         plt.legend(loc='center left',bbox_to_anchor=(0.80,0.80), fontsize="small")

         plt.grid(True)
         pdf_pages.savefig(fig)
   
def plot_uk_vs_t(uk, kt, hnur, tevlabl, uklabl):
# Plot uk(nkt,nrho,nhv) versus kt  for various densities and frequency bands

   for m in range(7,50,15):                  
      dum1 = uk[:,0,m]
      dum2 = uk[:,2,m]
      dum3 = uk[:,4,m]
      dum4 = uk[:,6,m]
      
      ymax = 1.5 * max(dum4)
      ymin = 0.9 * min(dum1)
      ymin2 = max(ymin, 1.e-10*ymax)
   
      hvtitl = 'hv: ' + str(hnur[m]) + ' to ' + str(hnur[m+1]) + ' eV'
   
      fig = plt.figure(figsize=(9.0, 6.5))
      plt.title(hvtitl)
      plt.xscale('log')
      plt.yscale('log')
      plt.ylim(ymin2, ymax)
      plt.xlabel(tevlabl)
      plt.ylabel(uklabl)
      plt.plot(kt, dum1, kt, dum2, kt, dum3, kt, dum4)
      plt.grid(True)
      pdf_pages.savefig(fig)
      
def plothyd(icmax, rc, rho, tkel, sie, uc, dr, mass, lbltitl, figname, itype, rmax):

   fig, (ax0, ax1) = plt.subplots(nrows=2,figsize=(9.0, 6.5))

   rmin = 0.95*rc[0]
   
#  Plot RHO vs RC and UC vs RC
 
   ax0.set_title(lbltitl, fontsize = 10)
   ax0.set_xlim(rmin,rmax)
   ax0.set_xscale('log')
   ax0.set_yscale('log')
   ax0.set_ylabel('RHO (gm/cc)', color='blue', fontsize=8)
   ax0.plot(rc,rho)
   
   ax2 = ax0.twinx()
   ax2.set_ylabel('UC (km/s)', color='red', fontsize=8)
   uckm = 0.001 * uc
   ax2.plot(rc, uckm, color='red')
   
# Plot TK vs RC and DR vs RC

   ax1.set_xscale('log')
   ax1.set_xlim(rmin,rmax)
   ax1.set_yscale('log')
   ax1.set_xlabel('RC (m)')
   ax1.set_ylabel('T (K)', color='blue', fontsize=8)
   ax1.plot(rc,tkel)
   
   ax2 = ax1.twinx()
   ax2.set_yscale('log')
   ax2.set_ylabel('Cell Mass (kg)', color='red', fontsize=8)
   ax2.plot(rc, mass, color='red')  
   
   fig.tight_layout()

   if(itype > 0):   savefig(figname)
   if(itype < 1):   pdf_pages.savefig(fig) 
   
# Main Program Starts Here ***********************************

set_globals()

baseall   = '/home/e.symbalisty.adm/RFLOCHM29/'
basealts  = '/home/e.symbalisty.adm/RFLOCHM40/T3200/'
basealts2 = '/home/e.symbalisty.adm/RFLOCHM40/T3200alum/' 
basename  = '/home/e.symbalisty.adm/RFLOCode/TESTNOCHEM/'
basename2 = '/home/e.symbalisty.adm/RFLOCode/TESTCHEM/'
basename3 = '/home/e.symbalisty.adm/HYCHEM2013AF/src/KT320_KM08_NOCHEM/'
basename4 = '/home/e.symbalisty.adm/HYCHEM2013AF/src/KT320_KM08_CHEM/'
dataname  = '/home/e.symbalisty.adm/DataFiles/RFLOCHM/'
eosname   = 'eos.txt'
opacname  = 'Hai_230818b_f1_rosseland_78g.txt'
opacname2 = 'aesop51_ns_rev1.txt'
opacname3 = 'mudatha51.txt'
comptitl  = '60 kt, 1.63 km, CHEM, LLNL RP 78 (red) vs AESOP 51 (blue)'

plot_all         = False
plot_eos         = False
plot_opac        = False
plot_opac_nubins = False
plot_opac_cmp    = False
plot_hydro_movie = False
plot_hydro_pdf   = False
plot_power       = False
plot_rtvar_hydro = False
plot_planck      = False
plot_compare     = True
plot_compare4    = False
plot_alts        = False
plot_alts_cmp    = False
plot_xdep        = False

if(plot_power):
   pdf_pages = PdfPages(basename + 'power.pdf')
   fname     = basename + 'pwr_ave.plt'
   npts      = rdpowpts(fname)
   pl_title, tpwr, pow_band, totpow = rdpower(npts, fname)
   print(pl_title)
   plotpwr(tpwr, pow_band, totpow, npts, pl_title)
   print("plot_power done")
   pdf_pages.close()

if(plot_all):
   ylds   = ['T0010', 'T0018', 'T0032', 'T0056', 'T0100', 'T0180', \
             'T0320', 'T0560', 'T1000', 'T1800', 'T3200', 'T5600', \
             'KT0010', 'KT0018', 'KT0032', 'KT0056', 'KT0100', \
             'KT0180', 'KT0320', 'KT0560', 'KT1000', 'KT1800', \
             'KT3200', 'KT5600', 'MT10', 'MT18' ]
   alts   = ['KM00','KM04','KM08','KM12','KM16','KM20','KM24', \
             'KM28','KM32','KM36','KM40','KM44','KM48','KM52', \
             'KM56','KM60','KM64','KM68','KM72','KM76','KM80' ]
   for j in range(0,26,1):
     pdf_pages = PdfPages(baseall + 'PIX_v29/' + 'power_' + ylds[j] + '.pdf')

     for i in range(0,21,1):
        fname  = baseall + ylds[j] + '/' + alts[i] + '/pwr_ave.plt'
        npts      = rdpowpts(fname)
        pl_title, tpwr, pow_band, totpow = rdpower(npts, fname)
        plotpwr(tpwr, pow_band, totpow, npts, pl_title)
        print(fname + " power plot done")

     pdf_pages.close()
     print(" ")

if(plot_alts):
   pdf_pages = PdfPages(basealts + 'power_alts.pdf')
   alts   = ['KM00','KM04','KM08','KM12','KM16','KM20','KM24', \
             'KM28','KM32','KM36','KM40','KM44','KM48','KM52', \
             'KM56','KM60','KM64','KM68','KM72','KM76','KM80' ]
   for i in range(0,21,1):
      fname  = basealts + alts[i] + '/pwr_ave.plt'
      npts      = rdpowpts(fname)
      pl_title, tpwr, pow_band, totpow = rdpower(npts, fname)
      plotpwr(tpwr, pow_band, totpow, npts, pl_title)
      print(fname + " power plot done")

   pdf_pages.close()

if(plot_alts_cmp):
   pdf_pages = PdfPages(basealts + 'power_alts_cmp.pdf')
   alts   = ['KM00','KM04','KM08','KM12','KM16','KM20','KM24', \
             'KM28','KM32','KM36','KM40','KM44','KM48','KM52', \
             'KM56','KM60','KM64','KM68','KM72','KM76','KM80' ]
   for i in range(0,21,1):
      fname   = basealts  + alts[i] + '/pwr_ave.plt'
      fname2  = basealts2 + alts[i] + '/pwr_ave.plt'
      npts    = rdpowpts(fname)
      pl_title, tpwr, pow_band, totpow = rdpower(npts, fname)

      npts2    = rdpowpts(fname2)
      pl_title2, tpwr2, pow_band2, totpow2 = rdpower(npts2, fname2)

      lably  = 'Thermal Power (W)'
      hvtitl = alts[i] + '  ' + comptitl
      plotpwr2(tpwr, totpow, npts, tpwr2, totpow2, npts2, lably, hvtitl)

      print(alts[i] + " power compare done")

   pdf_pages.close()
 
if(plot_xdep):
   pdf_pages = PdfPages(basename + 'xdep.pdf')
   fname  = basename + 'xdep.owt'
   pl_title, rc, temp, npts, rc2, temp2, npts2 = rdxdep(fname)
   print(pl_title)
   plotxdep(rc, temp, npts, rc2, temp2, npts2, pl_title)
   print("xdep plot done")
   pdf_pages.close()

if(plot_compare):
   pdf_pages = PdfPages(basename + 'power_compare.pdf')
   fname  = basename + 'pwr_ave.plt'
   npts   = rdpowpts(fname)
   pl_title, tpwr, pow_band, totpow = rdpower(npts, fname)

   fname2  = basename2 + 'pwr_ave.plt'
   npts2   = rdpowpts(fname2)
   pl_title2, tpwr2, pow_band2, totpow2 = rdpower(npts2, fname2)
   
   lably  = 'Thermal Power (W)'
   hvtitl = comptitl
   plotpwr2(tpwr, totpow, npts, tpwr2, totpow2, npts2, lably, hvtitl)
   print("plot compare done")   
   pdf_pages.close()

if(plot_compare4):
   pdf_pages = PdfPages(basename2 + 'power_compare4.pdf')

   fname  = basename + 'pwr_ave.plt'
   npts      = rdpowpts(fname)
   pl_title, tpwr, pow_band, totpow = rdpower(npts, fname)

   fname2  = basename2 + 'pwr_ave.plt'
   npts2     = rdpowpts(fname2)
   pl_title2, tpwr2, pow_band2, totpow2 = rdpower(npts2, fname2)

   fname3  = basename3 + 'ptdata'
   fname4  = basename4 + 'ptdata'
   pl_title3, npts3, tpwr3, pow_band3, totpow3, pblue3 = rdptdata(fname3)
   pl_title4, npts4, tpwr4, pow_band4, totpow4, pblue4 = rdptdata(fname4)
   
   lably  = 'Thermal Power (W)'

#  1 and 3 are no chemistry; 2 and 4 are with chemistry

#  Compare HYCHEM2013AF, RFLOCHM Thermal Powers with No Chem
   hvtitl = 'KT320 KM08 RFLOCHM (red) and HYCHEM2013AF (blue) No Chemistry'
   plotpwr2(tpwr, totpow, npts, tpwr3, totpow3, npts3, lably, hvtitl)

#  Compare HYCHEM2013AF, RFLOCHM Thermal Powers with Chem
   hvtitl = 'KT320 KM08 RFLOCHM (red) and HYCHEM2013AF (blue) with Chemistry'
   plotpwr2(tpwr2, totpow2, npts2, tpwr4, totpow4, npts4, lably, hvtitl)

#  Compare RFLOCHM Thermal Powers with and without Chem
   hvtitl = 'KT320 KM08 RFLOCHM No Chem (red) and with Chem (blue)'
   plotpwr2(tpwr, totpow, npts, tpwr2, totpow2, npts2, lably, hvtitl)

#  Compare HYCHEM2013AF Thermal Powers with and without Chem
   hvtitl = 'KT320 KM08 HYCHEM2013AF No Chem (red) and with Chem (blue)'
   plotpwr2(tpwr3, totpow3, npts3, tpwr4, totpow4, npts4, lably, hvtitl)

   print(basename2 + 'power_compare4.pdf' + ' done')   
   pdf_pages.close()

if(plot_rtvar_hydro):
   pdf_pages = PdfPages(basename + 'rtvar.pdf')
   fname  = basename + 'rtvar.owt'
   npts   = rdpowpts(fname)
   tpwr   = np.zeros((npts))
   ptherm = np.zeros((npts))
   ke     = np.zeros((npts))
   ie     = np.zeros((npts))
   erad   = np.zeros((npts))
   etot   = np.zeros((npts))
   rshk   = np.zeros((npts))
   ritmx  = np.zeros((npts))
   titmxm = np.zeros((npts))
   titmx  = np.zeros((npts))
   titmxp = np.zeros((npts))
   rhorat = np.zeros((npts))
   flupow = np.zeros((npts))
   rdrtvar(tpwr, ptherm, ke, ie, erad, etot, rshk, ritmx, titmxm, 
        titmx, titmxp, rhorat, flupow, npts, fname)
   plotrthyd(tpwr, ptherm, ke, ie, erad, etot, rshk, ritmx, titmxm, 
        titmx, titmxp, rhorat, npts) 
   print("rtvar plots done")  
   pdf_pages.close()

if(plot_eos):
   pdf_pages = PdfPages(dataname + 'eos.pdf')
   fname   = dataname + eosname
   rhotabl = np.zeros((nrho))
   etable  = np.zeros((nkt))
   fp      = np.zeros((nkt,nrho))
   gt      = np.zeros((nkt,nrho))
   tev     = np.zeros((nkt,nrho))
   rhotabl, etable, fp, gt, tev = rdeos(fname)
   ploteos(etable, rhotabl, fp, tev)
   print("eos plot done")
   pdf_pages.close()
   
if(plot_opac):
   pdf_pages = PdfPages(dataname + 'opacity.pdf')
   fname   = dataname + opacname
   rhotabl = np.zeros((nrho))
   kt      = np.zeros((nkt))
   hnur    = np.zeros((nhvp))
   uk      = np.zeros((nkt,nrho,nhv))
   rhotabl, kt, hnur, uk = rdopac(fname,[nkt,nrho,nhv], [nhv+1])
   hmid    = np.zeros((nhv))
   for i in range(0,nhv,1):
     hmid[i] = 0.5*(hnur[i] + hnur[i+1])
   plot_uk(uk, kt, hmid, rhotabl, 'T (eV)', 'UK (cm^2/gm)', nhv, nrho, nkt, opacname)
   print("opacity plot done")
   pdf_pages.close()

if(plot_opac_nubins):
   pdf_pages = PdfPages(dataname + 'opacity.pdf')
   fname   = dataname + opacname
   nrho = 27
   nkt  = 100
   nhv  = 78
   nhvp = nhv + 1
   rhotabl = np.zeros((nrho))
   kt      = np.zeros((nkt))
   hnur    = np.zeros((nhvp))
   uk      = np.zeros((nkt,nrho,nhv))
   rhotabl, kt, hnur, uk = rdopac2(fname,[nkt,nrho,nhv], [nhv+1])
   hmid    = np.zeros((nhv))
   for i in range(0,nhv,1):
     hmid[i] = 0.5*(hnur[i] + hnur[i+1])
   plot_uk(uk, kt, hmid, rhotabl, 'T (eV)', 'UK (cm^2/gm)', nhv, nrho, nkt, opacname)
   print("opacity plot done")
   pdf_pages.close()

if(plot_opac_cmp):
   pdf_pages = PdfPages(dataname + 'opacity_compare.pdf')
   fname   = dataname + opacname
   fname2  = dataname + opacname2
   fname3  = dataname + opacname3
   rhotabl = np.zeros((nrho))
   kt      = np.zeros((nkt))
   rhotabl, kt, hnur,   uk   = rdopac(fname,  [nkt,nrho,nhv], [nhv+1])
   rhotabl, kt, hnur2,  uk2  = rdopac(fname2, [nkt,nrho,nhv], [nhv+1])
   rhotabl, kt, hnur3,  uk3  = rdopac(fname3, [nkt,nrho,nhv], [nhv+1])   
   hmid    = np.zeros((nhv))
   hmid2   = np.zeros((nhv))
   hmid3   = np.zeros((nhv))
   for i in range(0,nhv,1):
     hmid[i]  = 0.5*(hnur[i]  + hnur[i+1])
     hmid2[i] = 0.5*(hnur2[i] + hnur2[i+1])
     hmid3[i] = 0.5*(hnur3[i] + hnur3[i+1])
   plot_uk_compare(uk, hmid, uk2, hmid2, uk3, hmid3, kt, rhotabl, nhv, nrho, nkt)
   print("opacity compare plot done")
   pdf_pages.close()
   
if(plot_planck):
   pdf_pages = PdfPages(basename + 'planck.pdf') 
   fname = basename + 'planck_table.owt'
   kt    = np.zeros((nkt))
   hnur  = np.zeros((nhvp))
   b     = np.zeros((nkt,nhv))
   b, hnur, kt = rdplanck(fname)
   hmid    = np.zeros((nhv))
   for i in range(0,nhv,1):
     hmid[i] = 0.5*(hnur[i] + hnur[i+1])
   for i in range(0,nhv,1):
     dhvinv = 1./(hnur[i+1] - hnur[i])
     b[:,i] = b[:,i] * dhvinv
   plot_b(b, hnur, kt, hmid) 
   print("planck plot done")    
   pdf_pages.close()

if(plot_hydro_movie) or (plot_hydro_pdf):
#  create a list of the dmp file names
   dmpnames = [ ]
   runtitl  = ' '
   mydmp  = open(basename +'hydrodata',"r")
   mydmp.readline()
   runtitl = mydmp.readline()
   mydmp.readline()
   print(runtitl)

#  read in the dmp file names
   while True:
      dum = mydmp.readline()
      if len(dum) == 0:
         break
      dmpnames.append(dum[0:10])
   ndumps = len(dmpnames)
   print("ndumps = " + str(ndumps))
   
if(plot_hydro_movie):
#  loop over the dmp files and draw some pix and save as a GIF movie
   movie_name = basename + 'hydro.gif'
   images = []
   for i in range(0,ndumps,1):
     fname = basename + dmpnames[i]
     mydmp = open(fname,"r")    
     icmax, rc,dr,rho,uc,pr,tkel,mass,sie,cs,gamma, time,ncycle,runlbl = rd_rflo_dmp(fname)      
     
     lbltitl = runlbl + '     T(s) =' + str(time) + ' , ncycle = ' + str(ncycle)
     figname = 'last.png'
     itype   = 1		# for movie output
     rmax    = rc.max()		# reset for zooming in
   
     plothyd(icmax, rc, rho, tkel, sie, uc, dr, mass, lbltitl, figname, itype,rmax)

     images.append(imageio.imread(figname))
     
     print(dmpnames[i])

   imageio.mimsave(movie_name, images, duration=0.1, loop=3)

   print("hydro movie done")      


if(plot_hydro_pdf):
#  loop over the dmp files and draw some pix and save as a PDF file
#  rc is in meters
   pdf_pages = PdfPages(basename + 'hydro.pdf')
   nstrt = 0
   nend  = ndumps
	
   for i in range(nstrt,nend,1):
     fname = basename + dmpnames[i]
     mydmp = open(fname,"r")    
     icmax, rc,dr,rho,uc,pr,tkel,mass,sie,cs,gamma, time,ncycle,runlbl = rd_rflo_dmp(fname)      
     
     lbltitl = runlbl + '     T(s) =' + str(time) + ' , ncycle = ' + str(ncycle)
     figname = ' '
     itype   = 0		# for pdf output
     rmax = rc.max()		# reset for zooming in
   
     plothyd(icmax, rc, rho, tkel, sie, uc, dr, mass, lbltitl, figname, itype, rmax)
     
     print(dmpnames[i])

   print("hydro pdf done")      
   pdf_pages.close() 


