# from pylab import *
import matplotlib
import numpy as np
import math
import imageio
import matplotlib.animation as animation
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties

matplotlib.rcParams['figure.max_open_warning'] = 0

def set_globals(ciraname, eltname):
   global nrho, nkt, nhv, nhvp
   global ncira, zkm, rhokm, tk_km, pr_km, natm0, pi
   global strk_mf12, strk_mfnone, sxx_w12, sxx_none, sxx_w47
   global z_elt, nkm_elt, beta_ro78, beta_rop78, exptau
   global hevent, hsensor
   global bands_250_9000
   
   nrho = 20
   nkt  = 100
   nhv  = 78
   nhvp = nhv + 1
   pi   = np.arccos(-1.0)

   hevent  = 0.
   hsensor = 0.

#  read cira data
   zkm, rhokm, tk_km, pr_km, ncira = rd_cira(ciraname)

#  compute natm0 = number of atmospheres looking straight up, used for normalization
   natm0 = natm0_set()
#   print("natm0 = ", natm0)

#  set streak camera film weights from LLNL (Spriggs) - first 42 bins
   strk_mf12   = np.zeros([42])
   strk_mfnone = np.zeros([42])
   sxx_w12     = np.zeros([42])
   sxx_none    = np.zeros([42])
   sxx_w47     = np.zeros([42])

   strk_mf12[23:27]   = [ 0.0858, 0.2445, 0.3448, 0.3249 ]
   strk_mfnone[23:34] = [ 0.0146, 0.0416, 0.0586, 0.0552, 0.0584, 0.0131, 
                          0.0160, 0.0510, 0.1298, 0.1974, 0.3643]

   sxx_w12[23:27]  = [ 0.10107, 0.37690, 0.27189, 0.25015 ]
   sxx_none[23:34] = [ 0.02004, 0.07474, 0.05392, 0.04961, 0.05376, 0.10862,
                       0.12229, 0.13112, 0.12941, 0.12825, 0.12825]

#  Bands 31 & 32 which go to 30 & 31 in python indexing
   sxx_w47[30:32] = 1.0

   bands_250_9000 = np.zeros([42])
   set_bands(bands_250_9000)

#  z_elt = Elterman altitudes in km, 0 - 51, nkm_elt = 51
#  frequency dependent transmission factors akin to 0.87

   z_elt, beta_ro78, beta_rop_78, nkm_elt = rd_elterman78(eltname)
   exptau = np.zeros([42])

def set_bands(bands_250_9000):

# New photon bin edges, 79 edges = 78 bins

#                 ev           nm         width in nm
#    1         0.10000  12398.42041016   8505.66674920
#    2         0.31850   3892.75366096    686.70783707
#    3         0.38672   3206.04582389    565.61213913
#    4         0.46956   2640.43368476    465.80632829
#    5         0.57014   2174.62735647    383.62102790
#    6         0.69226   1791.00632857    151.00362881
#    7         0.75600   1640.00269976     80.00164157
#    8         0.79477   1560.00105819     84.94687752
#    9         0.84054   1475.05418066    185.05705126
#   10         0.96112   1289.99712941     19.99247792
#   11         0.97625   1270.00465149     55.17593377
#   12         1.02059   1214.82871772    114.82956854
#   13         1.12713   1099.99914918     49.99742068
#   14         1.18080   1050.00172850     20.00161135
#   15         1.20373   1030.00011715     29.48200787
#   16         1.23920   1000.51810928     80.51511071
#   17         1.34765    920.00299856     30.00239486
#   18         1.39308    890.00060371      9.99972288
#   19         1.40891    880.00088083     55.98298872
#   20         1.50463    824.01789212     34.01773855
#   21         1.56942    790.00015357     29.99963804
#   22         1.63137    760.00051553     19.99953609
#   23         1.67546    740.00097944     70.00079516
#   24         1.85051    670.00018428     20.00042491
#   25         1.90745    649.99975937     89.99883750
#   26         2.21400    560.00092187     10.00134728
#   27         2.25426    549.99957459      9.99868565
#   28         2.29600    540.00088894     74.00172296
#   29         2.66061    465.99916599      7.99938725
#   30         2.70708    457.99977873     15.99938622
#   31         2.80507    442.00039251     28.00073276
#   32         2.99479    413.99965975     10.00019381
#   33         3.06892    403.99946594      7.99892903
#   34         3.13091    396.00053691     16.88336835
#   35         3.27034    379.11716856      4.11768545
#   36         3.30625    374.99948310      5.99997387
#   37         3.36001    368.99950923      7.99998354
#   38         3.43447    360.99952570     10.99993756
#   39         3.54241    349.99958814     37.76287223
#   40         3.97084    312.23671591     55.08223536
#   41         4.82139    257.15448056     45.36555152
#   42         5.85414    211.78892903     37.36203784
#   43         7.10809    174.42689119     30.77119423
#   44         8.63065    143.65569697     25.34259479
#   45        10.47933    118.31310218     20.87188550
#   46        12.72400     97.44121668     17.18989181
#   47        15.44949     80.25132487     24.13149009
#   48        22.09276     56.11983478     16.87519905
#   49        31.59265     39.24463573     11.80083482
#   50        45.17749     27.44380091      8.25233153
#   51        64.60381     19.19146937      5.77086167
#   52        92.38345     13.42060771      4.03556497
#   53       132.10830      9.38504273      2.82207686
#   54       188.91490      6.56296587      1.97347913
#   55       270.14830      4.58948674      1.38005571
#   56       386.31210      3.20943103      0.96507365
#   57       552.42630      2.24435738      0.67487687
#   58       789.96970      1.56948050      0.47194201
#   59      1129.65700      1.09753849      0.33002881
#   60      1615.40900      0.76750968      0.23078965
#   61      2310.03500      0.53672002      0.16139133
#   62      3303.35000      0.37532869      0.11286105
#   63      4723.79000      0.26246765      0.07892385
#   64      6755.02000      0.18354380      0.05008395
#   65      9290.00000      0.13345985      0.02224331
#   66     11148.00000      0.11121654      0.01853609
#   67     13377.60000      0.09268045      0.01544674
#   68     16053.12000      0.07723371      0.01287227
#   69     19263.74000      0.06436144      0.01072691
#   70     23116.49000      0.05363453      0.00893909
#   71     27739.79000      0.04469544      0.00744924
#   72     33287.75000      0.03724620      0.00620770
#   73     39945.30000      0.03103850      0.00517308
#   74     47934.36000      0.02586541      0.00431090
#   75     57521.23000      0.02155451      0.00359242
#   76     69025.48000      0.01796209      0.00299368
#   77     82830.57000      0.01496841      0.00256999
#   78    100000.00000      0.01239842      0.00197958
#   79    119000.00000      0.01041884   

   bands_250_9000[0]    = 0.600
   bands_250_9000[1:40] = 1.

def rd_cira(fname):
   mydmp = open(fname)
   ncira = 101
   mydmp.readline()
   mydmp.readline()

   zkm   = np.zeros([ncira])
   rhokm = np.zeros([ncira])
   tk_km = np.zeros([ncira])
   pr_km = np.zeros([ncira])

   for iline in range(0,ncira,1):
     line = mydmp.readline()  
     values = [float(s) for s in line.split()]
     zkm[iline]   = values[0]
     rhokm[iline] = values[1]
     tk_km[iline] = values[3]
     pr_km[iline] = values[4]

   return zkm, rhokm, tk_km, pr_km, ncira

def natm0_set():

   # Compute the number of atmospheres from the ground to the end of the modeled atmosphere
   # Looking straight up, used for normalization

   ycutoff = 90.001
   dkm     = 0.25
   natm0   = 0.

   npts    = 1 + int(ycutoff/dkm)			# 250 m spacing should be good enough
   ykm     = np.zeros([npts])
   s_rho   = np.zeros([npts])

   for i in range(0,npts):
      ykm[i] = dkm * float(i)				# vertical path	
      qr     = rhoval(ykm[i])
      s_rho[i] = qr					# density along path of integration 
#      print(i, qr, ykm[i], s_rho[i])

   for i in range(0,npts-1):
      ds    = ykm[i+1] - ykm[i]
      natm0 = natm0 + 0.5*ds*(s_rho[i+1] + s_rho[i])	# trapezoidal integration

   return natm0

def natm_event_ground(hsensor, hevent, srange):

   # This is for ground based sensors (not satellites) on a flat earth

   ycutoff = 90.001
   dkm     = 0.25
   if(hsensor <= hevent):
     z1 = hsensor
     z2 = hevent
   if(hsensor > hevent):
     z1 = hevent
     z2 = hsensor 

   qcos = (z2 - z1) / srange

   npts    = 1 + int(srange /dkm)
   skm     = np.zeros([npts])		
   ykm     = np.zeros([npts])
   s_rho   = np.zeros([npts])
   natm    = 0.

   for i in range(0,npts):
      skm[i] = dkm * float(i)				
      ykm[i] = z1 + qcos * skm[i]
      qr     = rhoval(ykm[i])
      s_rho[i] = qr					# density along path of integration 

   for i in range(0,npts-1):
      ds   = skm[i+1] - skm[i]
      natm = natm + 0.5*ds*(s_rho[i+1] + s_rho[i])	# trapezoidal integration

   natm = natm / natm0
   return natm

def natm_event_sat(hsensor, lka):

   # This is for satellites - assumes spherical earth
   # lka = look angle in degrees (0 is for overhead)

   r_earth = 6378.
   deg2rad = np.arccos(-1.0) / 180.
   qlka    = min(lka,90.)

   # compute srange which will run from ~ 100 to 1100 km depending on lookangle
   # find srange such that the altitude above the spherical earth is 99 km
   
   qcos = cos(qlka*deg2rad)
   qsin = sin(qlka*deg2rad)
   z1   = r_earth + hevent
   z2   = r_earth + 99.	

   srange = np.sqrt(z2*z2 - z1*z1*qsin*qsin) - z1*qcos

   npts   = 401
   dkm    = srange / float(npts)
   skm    = np.zeros([npts])
   s_rho  = np.zeros([npts])
   
   for i in range(0,npts,1):
     skm[i] = dkm * float(i)
     ykm    = hevent + qcos * skm[i]
     qr     = rhoval(ykm[i])
     s_rho[i] = qr					# density along path of integration 

   natm = 0.
   for i in range(0,npts-1):
      ds   = skm[i+1] - skm[i]
      natm = natm + 0.5*ds*(s_rho[i+1] + s_rho[i])

   natm = natm / natm0
   return natm

def rhoval(zval):

   # find rho by log interpolation from an atmospheric density profile
   # this assumes zkm are in 1 km steps from 0

   ival = int(zval + 0.00001)
   qx   = (zval - zkm[ival]) / (zkm[ival+1] - zkm[ival])
   qr   = (1. - qx)*np.log(rhokm[ival]) + qx*np.log(rhokm[ival+1])
   qr   = np.exp(qr)

   return qr 

def rd_elterman78(fname):
   # ro implies Rayleigh  and Ozone
   # rop implies Rayleigh and Ozone and Aerosols
   # there are 51 Elterman altitudes (0 - 51 km)

   # non zero values correspond to photon bins 10 to 40 (31 values)
   #    which correspond to indices 9 to 39
   #    written out with 11 values per line

   # we are currently writing out the first 42 (of 78) photon bins

   nkm_elt    = 51
   beta_ro78  = np.zeros([nkm_elt,42])
   beta_rop78 = np.zeros([nkm_elt,42])
   z_elt      = np.zeros([nkm_elt])

   mydmp = open(fname)

   for i in range(0,40,1): mydmp.readline()

   for i in range(0,nkm_elt,1):
     z_elt[i] = float(mydmp.readline())
     line     = mydmp.readline()
     values   = [float(s) for s in line.split()]
     line2     = mydmp.readline()
     values2   = [float(s) for s in line2.split()]
     line3    = mydmp.readline()
     values3   = [float(s) for s in line3.split()]
     beta_ro78[i][9:20]   = values[0:11]
     beta_ro78[i][20:31]  = values2[0:11]
     beta_ro78[i][31:40]  = values3[0:9]
    
   for i in range(0,3,1): mydmp.readline()

   for i in range(0,nkm_elt,1):
     z_elt[i] = float(mydmp.readline())
     line     = mydmp.readline()
     values   = [float(s) for s in line.split()]
     line     = mydmp.readline()
     values2   = [float(s) for s in line.split()]
     line     = mydmp.readline()
     values3   = [float(s) for s in line.split()]
     beta_rop78[i][9:20]   = values[0:11]
     beta_rop78[i][20:31]  = values2[0:11]
     beta_rop78[i][31:40]  = values3[0:9]     
  
   mydmp.close()
   return z_elt, beta_ro78, beta_rop78, nkm_elt

def attenuate78_srange(pow_band, qatm, hsensor, hevent, srange, nsim, atten):

   # nhv = 78, 42 bins written out
   # beta_ro = beta_ro78  implies Rayleigh & Ozone
   # beta_ro = beta_ro78p implies Rayleigh & Ozone + Aerosols
   # hsensor = sensor height above sea level in km
   # hevent  = source height above sea level in km
   # srange  = slant range in km
   # exptau  = exp(tau[k])

   if(hsensor <= hevent):
     z1 = hsensor
     z2 = hevent
   if(hsensor > hevent):
     z1 = hevent
     z2 = hsensor

   powadj = np.zeros([nsim,42])
   powadj = pow_band

   if(atten < 3):		# Elterman
      if(z2 < 50.):
         if (atten == 1): exptau = elt_event78_srange( z1, z2, srange, beta_ro78)
         if (atten == 2): exptau = elt_event78_srange( z1, z2, srange, beta_rop78)
         for k in range(9,40):
            for itime in range(0,nsim):
               powadj[itime][k] = exptau[k] * pow_band[itime][k]

   if(atten == 3):		# Simple clear air everywhere
      qa     = 0.87**qatm
      powadj = qa * powadj

   return powadj, exptau

def attenuate78_lka(pow_band, qatm, hevent, lka, nsim, atten):

   # this is for satellites
   # nhv = 78, 42 bins written out
   # beta_ro = beta_ro78  implies Rayleigh & Ozone
   # beta_ro = beta_ro78p implies Rayleigh & Ozone + Aerosols
   # lka     = look angle (0 implies overhead)
   # exptau  = exp(tau[k])

   powadj = np.zeros([nsim,42])
   powadj[0:nsim,0:42] = pow_band[0:nsim, 0:42]

   if(atten < 3):		# Elterman
      if(atten == 1): exptau = elt_event78_lka(hevent, lka, beta_ro78)
      if(atten == 2): exptau = elt_event78_lka(hevent, lka, beta_rop78)
      for k in range(9,40):
         for itime in range(0,nsim):
            powadj[itime][k] = exptau[k] * pow_band[itime][k]

   if(atten == 3):		# Simple clear air everywhere
      qa     = 0.87**qatm
      powadj = qa * powadj

   return powadj, exptau


def elt_event78_srange(z1, z2, srange, beta):

   # This is for ground, or near ground sensors - not satellites
   # Assumes a flat earth
   # Elterman grid has 1 km vertical spacing from 0 to 50 km
   #     Therefore, no attenuation above 50 km
   # srange = slant range from sensor to source
   # exptau  = exp(tau[k])

   npts   = 101
   dkm    = srange/float(npts-1)
   skm    = np.zeros([npts])
   bval   = np.zeros([npts,42])		# beta values for each photon bin
   tau    = np.zeros([42])
   exptau = np.zeros([42])
   qcos   = (z2 - z1)/srange		# cosine of look angle
   
   for i in range(0,npts,1):
     skm[i] = dkm * float(i)
     ykm    = z1 + qcos*skm[i]
     if(ykm >= 50.): bval[i][0:42] = 0.
     if(ykm <  50.):
        iykm = int(ykm)
        qx   = ykm - float(iykm)
        qxm  = 1. - qx
        for k in range(9,40):
           bval[i][k] = qxm * beta[iykm][k] + qx * beta[iykm+1][k]

   for k in range(9,40):
      tau[k] = 0.
      for i in range(0,npts-1):
         ds = skm[i+1] - skm[i]
         tau[k] = tau[k] + 0.5*ds*(bval[i][k] + bval[i+1][k])

   exptau = exp(-tau)
   return exptau

def elt_event78_lka(hevent, lka, beta):

   # This is for satellites - assumes spherical earth
   # lka = look angle in degrees (0 is for overhead)
   # Elterman grid has 1 km vertical spacing from 0 to 50 km
   #     Therefore, no attenuation above 50 km
   # exptau  = exp(tau[k])

   r_earth = 6378.
   deg2rad = np.arccos(-1.0) / 180.
   qlka    = min(lka,90.)

   # compute srange which will run from ~ 50 to 800 km depending on lookangle
   # find srange such that the altitude above the spherical earth is 50 km
   
   qcos = cos(qlka*deg2rad)
   qsin = sin(qlka*deg2rad)
   z1   = r_earth + hevent
   z2   = r_earth + 50.			# Maximum Elterman altitude

   srange = np.sqrt(z2*z2 - z1*z1*qsin*qsin) - z1*qcos

   npts   = 201
   dkm    = srange / float(npts)
   skm    = np.zeros([npts])
   bval   = np.zeros([npts,42])		# beta values for each photon bin
   tau    = np.zeros([42])
   exptau = np.zeros([42])
   
   for i in range(0,npts,1):
     skm[i] = dkm * float(i)
     ykm    = hevent + qcos * skm[i]
     if(ykm >= 50.): bval[i][0:42] = 0.
     if(ykm <  50.):
        iykm = int(ykm)
        qx   = ykm - float(iykm)
        qxm  = 1. - qx
        for k in range(9,40):
           bval[i][k] = qxm * beta[iykm][k] + qx * beta[iykm+1][k]

   for k in range(9,40):
      tau[k] = 0.
      for i in range(0,npts-1):
         ds = skm[i+1] - skm[i]
         tau[k] = tau[k] + 0.5*ds*(bval[i][k] + bval[i+1][k])

   exptau = exp(-tau)
   return exptau

def rd_streak(fname):
   mydmp = open(fname)
#  relative power p/(p2 max)

   datalabl = mydmp.readline()
   for i in range(1,11): mydmp.readline()

   npts = int(mydmp.readline())

   tval = np.zeros([npts])
   pval = np.zeros([npts])
   for iline in range(0,npts,1):
     line = mydmp.readline()
     values = [float(s) for s in line.split()]
     tval[iline] = values[0]
     pval[iline] = values[1]

   return tval, pval, npts, datalabl

def rd_ha(fname):
   mydmp = open(fname)

   # header stuff
   for i in range(0,5,1): mydmp.readline()
   npts = int(mydmp.readline())
   mydmp.readline()
   mydmp.readline()

#  time in sec, irradiance in Cal/(sec * cm^2), power in Watts

   tval = np.zeros([npts])
   pval = np.zeros([npts])

   for iline in range(0,npts,1):
     line = mydmp.readline()
     values = [float(s) for s in line.split()]
     tval[iline] = values[0]
     pval[iline] = values[2]

   return tval, pval, npts
   
def norm_sim(tpwr, psim, nsim, tchk):
   pnorm = np.zeros([nsim])
   pmax = 0.
   for i in range(0,nsim,1):
     if (tpwr[i] > tchk and psim[i] > pmax): pmax = psim[i]

   qa = 1. / pmax
   for i in range(0,nsim,1): pnorm[i] = qa * psim[i]

   return pnorm

def rdpower78(fname):
   # We are currently writing out the first 42 bins
   mydmp = open(fname)

   mydmp.readline()
   run_titl = mydmp.readline()
   mydmp.readline()

   mydmp.readline()
   line = mydmp.readline()
   values = line.split()
   npts   = int(values[0]) 		# Number of time points

   mydmp.readline()
   line = mydmp.readline()
   values = line.split()
   nbins  = int(values[0]) 		# Number of frequency bins (= 42)

   # set arrays

   tpwr       = np.zeros([npts])
   pow_band78 = np.zeros([npts,nbins])
   totpow78   = np.zeros([npts])

   # determine how many lines to read with 12 values per line

   qa     = float(npts)/12.
   qb     = float(int(qa))
   addone = 0
   if(qa-qb > 0.05): addone = 1
   nmax = int(qa) + addone

   # read in the time array
   mydmp.readline()
   ibeg = 0
   for iline in range(0,nmax,1):
      iend   = min(ibeg+12, npts+1)
      ilast  = iend - ibeg
      line   = mydmp.readline()
      values = [float(s) for s in line.split()]
      tpwr[ibeg:iend] = values[0:ilast]
      ibeg += 12

   # read in the first nbins bands in a 78 bin simulation
   for m in range(0,nbins,1):
      mydmp.readline()
      ibeg = 0
      for iline in range(0,nmax,1):
         iend   = min(ibeg+12, npts+1)
         ilast  = iend - ibeg
         line   = mydmp.readline()
         values = [float(s) for s in line.split()]
         pow_band78[ibeg:iend,m] = values[0:ilast]
         ibeg += 12

   mydmp.close()

   # compute the thermal power - Sum of the first 42 bins
   for i in range(0,npts,1):
      totpow78[i] = 0.
      for m in range(0,nbins,1):
         totpow78[i] = totpow78[i] + pow_band78[i,m]

   return run_titl, tpwr, pow_band78, totpow78, npts, nbins      

   
def plotpwr(tpwr, pow_band, totpow, npts, pl_titl):

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
   ymax = 2.0 * ymax
   ymin = 0.8 * ymin
   
   thyld = 0.
   for i in range(0,npts-1):
     thyld = thyld + (tpwr[i+1]-tpwr[i])*(totpow[i+1]+totpow[i])/2.
   thyld = thyld / 4.185e12
   tyld_labl = ("Thermal Yield (kt) = %.3f" % (thyld))   
   hvtitl    = pl_titl
   
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
   plt.gcf().text(.14, .85, tyld_labl, fontsize=10, color='blue')
   pdf_pages.savefig(fig)
 
def plotpwr2(tpwr, sipow, npts, tpwr2, sipow2, npts2, lably, hvtitl, dname, 
       sim_titl, toff_labl, range_labl, natm_labl, labl_amp, sensor_labl, wgts):
   xmax1  = max(tpwr)
   xmax2  = max(tpwr2)
   xmax   = 3.0 * max(xmax1, xmax2)
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
   ymax = 2.0 * ymax
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
#   plt.grid(True)
   plt.gcf().text(.13, .86, sensor_labl, fontsize=8, color='blue')
   plt.gcf().text(.13, .84, dname,       fontsize=8, color='blue')
   plt.gcf().text(.13, .82, toff_labl,   fontsize=8, color='blue')
   qdiff       = hevent - hsensor
   labl_hob    = 'HOB    AGL (km) = %6.3f' %qdiff
   labl_sensor = 'Sensor ASL (km) = %6.3f '%hsensor
   labl_sim    = 'RFLOCHM: ' + sim_titl
   plt.gcf().text(.13, .80, labl_hob,    fontsize=8, color='black')
   plt.gcf().text(.13, .78, labl_sensor, fontsize=8, color='black')
   plt.gcf().text(.13, .76, range_labl,  fontsize=8, color='black')
   plt.gcf().text(.13, .74, natm_labl,   fontsize=8, color='black')
   plt.gcf().text(.13, .72, labl_amp,    fontsize=8, color='black')
   plt.gcf().text(.36, .86, labl_sim,    fontsize=8, color='red')

   yvl = 0.86
   plt.gcf().text(.825, yvl, 'Band     Wgt      T', fontsize=7, color='black')
   for m in range(0,42,1):
      if(wgts[m] > 0.):
          yvl = yvl - 0.015
          labl = '  %2.0f' %float(m+1) + '   %5.3f' %wgts[m] + '   %5.3f' %exptau[m]
          plt.gcf().text(.835, yvl, labl,fontsize=7, color='black')


   pdf_pages.savefig(fig)

def plot2sim(tsim1, psim1, nsim1, tsim2, psim2, nsim2, lably, pl_titl, sub_titl):

   xmax1 = tsim1.max()
   xmax2 = tsim2.max()
   xmax  = 1.1 * max(xmax1, xmax2)
   xmin  = 5.e-7

   ymin = 1.e30
   ymax = 0.
   for i in range(0,nsim1):
     if(tsim1[i] > xmin):
       if(psim1[i] > ymax): ymax = psim1[i]
       if(psim1[i] < ymin): ymin = psim1[i]

   for i in range(0,nsim2):
     if(tsim2[i] > xmin):
       if(psim2[i] > ymax): ymax = psim2[i]
       if(psim2[i] < ymin): ymin = psim2[i]

   yminchk = 1.e-4 * ymax
   if(ymin < yminchk): ymin = yminchk
   ymax = 1.2 * ymax
   ymin = 0.8 * ymin
   
   thyld1 = 0.
   for i in range(0,nsim1-1):
     thyld1 = thyld1 + (tsim1[i+1]-tsim1[i])*(psim1[i+1]+psim1[i])/2.
   thyld1 = thyld1 / 4.185e12
   tyld_labl1 = ("Thermal Yield (kt) = %.3f" % (thyld1))   

   thyld2 = 0.
   for i in range(0,nsim2-1):
     thyld2 = thyld2 + (tsim2[i+1]-tsim2[i])*(psim2[i+1]+psim2[i])/2.
   thyld2 = thyld2 / 4.185e12
   tyld_labl2 = ("Thermal Yield (kt) = %.3f" % (thyld2))  

   hvtitl    = pl_titl
   
   fig = plt.figure(figsize=(9.0, 6.5))
   plt.title(hvtitl)
   plt.xscale('log')
   plt.yscale('log')
   plt.xlim(xmin, xmax)
   plt.ylim(ymin, ymax)
   plt.xlabel('T (s)')
   plt.ylabel(lably)
   plt.plot(tsim1, psim1, color='red')
   plt.plot(tsim2, psim2, color='blue')
   # plt.grid(True)
   plt.gcf().text(.14, .85, sub_titl,   fontsize=8, color='black')
   plt.gcf().text(.14, .83, tyld_labl1, fontsize=8, color='red')
   plt.gcf().text(.14, .81, tyld_labl2, fontsize=8, color='blue')
   pdf_pages.savefig(fig)  
    
def t2match(tsim, psim, nsim, tdat, pdat, ndat, tchk):
   pamp = np.zeros([nsim])

   pdat_max = 0.
   for i in range(0,ndat,1):
      if(tdat[i] > tchk):
         if(pdat[i] > pdat_max): pdat_max = pdat[i]

   psim_max = 0.
   for i in range(0,nsim,1):
      if(tsim[i] > tchk):
         if(psim[i] > psim_max): psim_max = psim[i]

   qamp = pdat_max / psim_max
   pamp = qamp * psim

   return pamp, qamp

def ptot_calc(wgts, power, nsim, nbins):
   #  sum the weighted bands
   ptot = np.zeros([nsim])
   for i in range(0,nsim,1):
       for m in range(0,nbins,1):
          ptot[i] = ptot[i] + wgts[m] * power[i,m]

   return ptot

def print_sim_powers(tsim, pthermal, psensor, nsim, dirname, caselabl, datalabl, simname):

   fname = dirname +  '/' + caselabl + '_sim_data.txt'
   f = open(fname,"w+")

   f.write(datalabl)
   f.write("RFLOCHM Simulation Data of Event in line above \n")
   f.write(simname + " \n")
   f.write("NPTS below\n")
   f.write("%d \n" %nsim)
   f.write("Thermal Power (W) has NO sensor weights & NO atmospheric attenuation \n")
   f.write("Sensor Power (W) has sensor weights & atmospheric attneuation \n\n")
   f.write(" T (s)         Thermal (W)    Sensor (W) \n")
   for i in range(0,nsim,1):
     dum = '%12.5e' %tsim[i] + '  %12.5e' %pthermal[i] + '  %12.5e \n' %psensor[i]
     f.write(dum)
     
   f.close()

def toff_calc(tsim, psim, nsim, tdat, pdat, ndat, pfact):

   # this is called after the data and simulation have been normalized at t2max
   # if pfact < 1.  compute toffsdet
   #          >=1.  set toffset = 0, write label, and return

   toffset   = 0.

   if(pfact < 1.):
      itest     = 0
      ptest_sim = 0.
      ttest_sim = 0.

      pmax_sim = psim.max()
      pmax_dat = pdat.max()
      ptest    = pfact * pmax_sim			# This may be problem dependent

      for i in range(0,nsim,1):
         if(psim[i] > ptest):
            ttest_sim = tsim[i]
            ptest_sim = psim[i]
            break

      for i in range(0,ndat,1):
         if(pdat[i] > ptest_sim):
            itest = i
            break

      if(tdat[itest] < 0.003):
         qx      = (ptest_sim - pdat[itest-1])/(pdat[itest] - pdat[itest-1])
         ttest   = (1. - qx) * tdat[itest-1] + qx * tdat[itest]
         toffset = ttest_sim - ttest

#      print(pdat[itest-1], ptest_sim, pdat[itest])
#      print(tdat[itest-1]+toffset, ttest_sim, tdat[itest]+toffset)

   toff_labl = 'Data T offset (s) = %8.6f' %toffset
   tdat      = tdat + toffset
   return toffset, toff_labl, tdat

def rd_magpie_eos(fname):
   mydmp = open(fname,"r")
   mydmp.readline()
   nr  = int(mydmp.readline())
   mkt = int(mydmp.readline())

#  check dimensions
#   print(nr, nrho, mkt, nkt)

   if((nr-nrho) != 0):
     print("nr ne nrho in rd_magpie_eos")
     return
   if((mkt-nkt) != 0):
     print("mkt ne nkt in rd_magpie_eos")
     return
   
#  RHOTABL[nrho] with 7 values per line to read in

   mydmp.readline()
   ibeg = 0
   for iline in range(0,3,1):
      iend  = min(ibeg+7,nr)
      ilast = iend - ibeg
      line  = mydmp.readline()
      values = [float(s) for s in line.split()]
      rhotabl[ibeg:iend] = values[0:ilast]
      ibeg += 7
   
   print (rhotabl)
   
#  ETABLE[nkt] with 7 values per line 

   mydmp.readline()
   ibeg = 0
   for iline in range(0,15,1):
     iend  = min(ibeg+7, nkt)
     ilast = iend - ibeg
     line  = mydmp.readline()
     values = [float(s) for s in line.split()]
     etable[ibeg:iend] = values[0:ilast]
     ibeg += 7
     
   print(etable)
  
# GT[nkt,nrho] array with 7 values per line

   mydmp.readline()
   for ikt in range(0,nkt,1):
      ibeg = 0
      for iline in range(0,3,1):
        iend  = min(ibeg+7, nr)
        ilast = iend - ibeg
        line  = mydmp.readline()
        values = [float(s) for s in line.split()]
        gt[ikt][ibeg:iend] = values[0:ilast]
        ibeg += 7
      
   print("GT first and last: ", gt[0][0], gt[nkt-1][nrho-1])
       
# FP[nkt,nrho] array with 7 values per line

   mydmp.readline()
   for ikt in range(0,nkt,1):
      ibeg = 0
      for iline in range(0,3,1):
        iend  = min(ibeg+7, nr)
        ilast = iend - ibeg
        line  = mydmp.readline()
        values = [float(s) for s in line.split()]
        fp[ikt][ibeg:iend] = values[0:ilast]
        ibeg += 7
      
   print("Gamma-1 first and last: ", fp[0][0], fp[nkt-1][nrho-1]) 
   
# T[nkt,nrho] in eV array with 7 values per line 

   mydmp.readline()
   for ikt in range(0,nkt,1):
      ibeg = 0
      for iline in range(0,3,1):
        iend  = min(ibeg+7, nr)
        ilast = iend - ibeg
        line  = mydmp.readline()
        values = [float(s) for s in line.split()]
        tev[ikt][ibeg:iend] = values[0:ilast]
        ibeg += 7
      
   print("Tev first and last: ", tev[0][0], tev[nkt-1][nrho-1])      
   
   mydmp.close()
   return rhotabl, etable, fp, gt, tev

def rd_alum_eos(fname):
   mydmp = open(fname,"r")
   mydmp.readline()
   nr  = int(mydmp.readline())
   mkt = int(mydmp.readline())

#  check dimensions
#   print(nr, nrho, mkt, nkt)

   if((nr-nrho) != 0):
     print("nr ne nrho in rd_alum_eos")
     return
   if((mkt-nkt) != 0):
     print("mkt ne nkt in rd_alum_eos")
     return
   
#  RHOTABL[nrho] with 10 values per line to read in

   mydmp.readline()
   ibeg = 0
   for iline in range(0,2,1):
      iend  = min(ibeg+10,nr)
      ilast = iend - ibeg
      line  = mydmp.readline()
      values = [float(s) for s in line.split()]
      rhotabl[ibeg:iend] = values[0:ilast]
      ibeg += 10
   
   print (rhotabl)
   
#  ETABLE[nkt] with 7 values per line 

   mydmp.readline()
   ibeg = 0
   for iline in range(0,15,1):
     iend  = min(ibeg+7, nkt)
     ilast = iend - ibeg
     line  = mydmp.readline()
     values = [float(s) for s in line.split()]
     etable[ibeg:iend] = values[0:ilast]
     ibeg += 7
     
   print(etable)
  
# GT[nkt,nrho] array with 10 values per line

   mydmp.readline()
   for ikt in range(0,nkt,1):
      ibeg = 0
      for iline in range(0,2,1):
        iend  = min(ibeg+10, nr)
        ilast = iend - ibeg
        line  = mydmp.readline()
        values = [float(s) for s in line.split()]
        gt[ikt][ibeg:iend] = values[0:ilast]
        ibeg += 10
      
   print("GT first and last: ", gt[0][0], gt[nkt-1][nrho-1])
       
# FP[nkt,nrho] array with 10 values per line

   mydmp.readline()
   for ikt in range(0,nkt,1):
      ibeg = 0
      for iline in range(0,2,1):
        iend  = min(ibeg+10, nr)
        ilast = iend - ibeg
        line  = mydmp.readline()
        values = [float(s) for s in line.split()]
        fp[ikt][ibeg:iend] = values[0:ilast]
        ibeg += 10
      
   print("Gamma-1 first and last: ", fp[0][0], fp[nkt-1][nrho-1]) 
   
# cs[nkt,nrho] in eV array with 10 values per line 

   mydmp.readline()
   for ikt in range(0,nkt,1):
      ibeg = 0
      for iline in range(0,2,1):
        iend  = min(ibeg+10, nr)
        ilast = iend - ibeg
        line  = mydmp.readline()
        values = [float(s) for s in line.split()]
        cs[ikt][ibeg:iend] = values[0:ilast]
        ibeg += 10
      
   print("CS first and last: ", cs[0][0], cs[nkt-1][nrho-1])      
   
   mydmp.close()
   return rhotabl, etable, fp, gt, cs
   
def ploteos(etable, rhotabl, fp, tev, cs, titl):

   colorvar = cm.rainbow(np.linspace(0,1,nrho))
   
# plot gamma-1 vs sie for various densities

   ymax  = 0.7
   ymin  = 0.0
   
   fig = plt.figure(figsize=(9.0, 6.5))
   plt.title(titl)
   plt.xscale('log')
   plt.yscale('linear')
   plt.ylim(ymin, ymax)
   plt.xlabel('SIE (erg/gm)')
   plt.ylabel('Gamma - 1')

   for j in range(0,nrho,1):
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
   plt.title(titl)
   plt.xscale('log')
   plt.yscale('linear')
   plt.ylim(ymin, ymax)
   plt.xlabel('SIE (erg/gm)')
   plt.ylabel('(Gamma+1) / (Gamma-1)')

   for j in range(0,nrho,1):
      gmm  = fp[:,j]
      for i in range(0,nkt,1): 
         if (gmm[i] < 0.001): gmm[i] = 0.001
      yval = (fp[:,j] + 2.)/gmm
      if (j == 0): plt.plot(etable, yval, color=colorvar[j], label=str(rhotabl[j]) + " gm/cc")
      if (j > 0):  plt.plot(etable, yval, color=colorvar[j], label=rhotabl[j])

   plt.legend(loc='center left',bbox_to_anchor=(0.80,0.50), fontsize="small")

   plt.grid(True)
   pdf_pages.savefig(fig)
   
# Plot Tev versus Sie for various densities

   ymax  = 1.1 * tev.max()
   ymin  = 0.9 * tev.min()
   
   fig = plt.figure(figsize=(9.0, 6.5))
   plt.title(titl)
   plt.xscale('log')
   plt.yscale('log')
   plt.ylim(ymin, ymax)
   plt.xlabel('SIE (erg/gm)')
   plt.ylabel('T (eV)')

   for j in range(0,nrho,1):
      yval = tev[:,j]
      if (j == 0): plt.plot(etable, yval, color=colorvar[j], label=str(rhotabl[j]) + " gm/cc")
      if (j > 0):  plt.plot(etable, yval, color=colorvar[j], label=rhotabl[j])

   plt.legend(loc='center left',bbox_to_anchor=(0.80,0.40), fontsize="small")

   plt.grid(True)
   pdf_pages.savefig(fig) 
 
# Plot Sound Speed versus Sie for various densities

   ymax  = 5000.
   ymin  = 0.10
   
   fig = plt.figure(figsize=(9.0, 6.5))
   plt.title(titl)
   plt.xscale('log')
   plt.yscale('log')
   plt.ylim(ymin, ymax)
   plt.xlabel('SIE (erg/gm)')
   plt.ylabel('Sound Speed (km/s)')

   for j in range(0,nrho,1):
      yval = cs[:,j]
      if (j == 0): plt.plot(etable, yval, color=colorvar[j], label=str(rhotabl[j]) + " gm/cc")
      if (j > 0):  plt.plot(etable, yval, color=colorvar[j], label=rhotabl[j])

   plt.legend(loc='center left',bbox_to_anchor=(0.80,0.40), fontsize="small")

   plt.grid(True)
   pdf_pages.savefig(fig) 

def rdopac78(fname):

   # This reads (100,20,78) opacity tables

   print(fname)
   mydmp = open(fname,"r")
   mydmp.readline()
   ikt  = int(mydmp.readline())
   irho = int(mydmp.readline())
   ihv  = int(mydmp.readline())

   if(ikt != nkt):
      print("ikt ne nkt in rdopac78")
      return
   if(irho != nrho):
      print("irho ne nrho in rdopac78")
      return
   if(ihv != nhv):
      print("ihv ne nhv in rdopac78")
      return

   print(nkt,nrho,nhv)

   uk      = np.zeros([nkt,nrho,nhv])
   rhotabl = np.zeros([nrho])
   kt      = np.zeros([nkt])
   hnur    = np.zeros([nhv+1])

# 100 kt bins, 20 density bins, 78 frequency bins
   
# RHOTABL(nrho), 7 values per line

   mydmp.readline()
   ibeg = 0
   for iline in range(0,3,1):
     iend  = min(ibeg+7, nrho)
     ilast = iend - ibeg
     line  = mydmp.readline()
     values = [float(s) for s in line.split()]
     rhotabl[ibeg:iend] = values[0:ilast]
     ibeg += 7
   
   print(rhotabl)
   
# KT(nkt) table, 7 values per line 

   mydmp.readline()
   ibeg = 0
   for iline in range(0,15,1):
     iend  = min(ibeg+7, nkt+1)
     ilast = iend - ibeg
     line  = mydmp.readline()
     values = [float(s) for s in line.split()]
     kt[ibeg:iend] = values[0:ilast]
     ibeg += 7
     
   print(kt)      
   
# HNUR(nhvp) table, 7 values per line for etable

   mydmp.readline()
   ibeg = 0
   for iline in range(0,12,1):
     iend  = min(ibeg+7, nhvp)
     ilast = iend - ibeg
     line  = mydmp.readline()
     values = [float(s) for s in line.split()]
     hnur[ibeg:iend] = values[0:ilast]
     ibeg += 7
     
   print(hnur) 
   
# UK(nkt,nrho,nhv), 7 values per line

   mydmp.readline()			# uk label
   
   for m in range(0,78,1):
     mydmp.readline()   		# m label

     ukrho = np.zeros([nrho])     
     for iline in range(0,100,1):

       ibeg = 0
       for irho in range(0,3,1):
         iend  = min(ibeg+7, nrho)
         ilast = iend - ibeg
         line  = mydmp.readline()
         values = [float(s) for s in line.split()]
         ukrho[ibeg:iend] = values[0:ilast]
         ibeg += 7

       uk[iline,0:nrho,m] = ukrho[0:nrho]
            
   print("UK first and last: ", uk[0][0][0], uk[nkt-1][nrho-1][nhv-1])
   
   mydmp.close()
   return nkt, nrho, nhv, rhotabl, kt, hnur, uk

def rd_deb78(fname, jkt, jrho, jhv):

   # This reads (78,30,78) aluminum opacity tables

   print(fname)
   mydmp = open(fname,"r")
   mydmp.readline()
   ikt  = int(mydmp.readline())
   irho = int(mydmp.readline())
   ihv  = int(mydmp.readline())

   if(ikt != jkt):
      print("ikt ne jkt in rd_deb78")
      return
   if(irho != jrho):
      print("irho ne jrho in rd_deb78")
      return
   if(ihv != jhv):
      print("ihv ne jhv in rd_deb78")
      return

   print(jkt,jrho,jhv)

   uk      = np.zeros([jkt,jrho,jhv])
   rhotabl = np.zeros([jrho])
   kt      = np.zeros([jkt])
   hnur    = np.zeros([jhv+1])

# 78 kt bins, 30 density bins, 78 frequency bins
   
# RHOTABL(jrho), 10 values per line

   mydmp.readline()
   ibeg = 0
   for iline in range(0,3,1):
     iend  = min(ibeg+10, jrho)
     ilast = iend - ibeg
     line  = mydmp.readline()
     values = [float(s) for s in line.split()]
     rhotabl[ibeg:iend] = values[0:ilast]
     ibeg += 10
   
   print(rhotabl)
   
# KT(jkt) table, 7 values per line 

   mydmp.readline()
   ibeg = 0
   for iline in range(0,12,1):
     iend  = min(ibeg+7, jkt+1)
     ilast = iend - ibeg
     line  = mydmp.readline()
     values = [float(s) for s in line.split()]
     kt[ibeg:iend] = values[0:ilast]
     ibeg += 7
     
   print(kt)      
   
# HNUR(jhv+1) table, 7 values per line for etable

   mydmp.readline()
   ibeg = 0
   for iline in range(0,12,1):
     iend  = min(ibeg+7, jhv+1)
     ilast = iend - ibeg
     line  = mydmp.readline()
     values = [float(s) for s in line.split()]
     hnur[ibeg:iend] = values[0:ilast]
     ibeg += 7
     
   print(hnur) 
   
# UK(jkt,jrho,jhv), 7 values per line

   mydmp.readline()			# uk label
   
   for m in range(0,jhv,1):
     mydmp.readline()   		# m label

     ukrho = np.zeros([jrho])     
     for iline in range(0,jkt,1):

       ibeg = 0
       for irho in range(0,3,1):
         iend  = min(ibeg+10, jrho)
         ilast = iend - ibeg
         line  = mydmp.readline()			# 10 per line
         values = [float(s) for s in line.split()]
         ukrho[ibeg:iend] = values[0:ilast]
         ibeg += 10

       uk[iline,0:jrho,m] = ukrho[0:jrho]
            
   print("UK first and last: ", uk[0][0][0], uk[jkt-1][jrho-1][jhv-1])
   
   mydmp.close()
   return rhotabl, kt, hnur, uk

def plot_uk(uk, kt, hmid, rhotabl, tevlabl, uklabl, jhv, jrho, jkt, opactitl):

# Plot uk(nkt,nrho,nhv) versus hv  

   colorvar = cm.rainbow(np.linspace(0,1,jrho))
   dum = np.zeros((jrho,jhv))

   for i in range(0,jkt,1):
      ipr = i + 1
      for j in range(0,jrho,1): dum[j,:] = uk[i,j,:]                  
      
      ymax  = 2.0 * dum.max()
      ymin2 = 0.8 * dum.min()
      ymin  = max(ymin2, 1.e-12*ymax)
   
      hvtitl   = opactitl + ",     T(eV)= %9.3f" %(kt[i])
   
      fig = plt.figure(figsize=(9.0, 6.5))
      plt.title(hvtitl)
      plt.xscale('log')
      plt.yscale('log')
      plt.xlim(0.1,1.e6)
      plt.ylim(ymin, ymax)
      plt.xlabel('Photon Energy (eV)')
      plt.ylabel('Mean Opacity (cm^2/gm)')
      plt.gcf().text(.86, .89, 'gm/cc', fontsize=10, color='black')

      for j in range(0,jrho,1):
         yval = dum[j,:]
         plt.plot(hmid, yval, color=colorvar[j], label=rhotabl[j] )

      plt.legend(loc='center left',bbox_to_anchor=(0.90,0.50), fontsize="small")

      plt.grid(True)
      pdf_pages.savefig(fig)
      print("T[ ", ipr, "]  Done")

def plot_uk_hist(uk, kt, hnur, rhotabl, tevlabl, uklabl, jhv, jrho, jkt, opactitl):

# Plot uk(nkt,nrho,nhv) versus hv  

   colorvar = cm.rainbow(np.linspace(0,1,jrho))
   dum = np.zeros((jrho,jhv))

   xhist = np.zeros([2*jhv])
   ipt = -1
   for k in range(0,jhv):
     ipt = ipt + 1
     xhist[ipt] = hnur[k]
     ipt = ipt + 1
     xhist[ipt] = hnur[k+1]

   # for i in range(0,2,1):
   for i in range(0,jkt,1):
      ipr = i + 1
      for j in range(0,jrho,1): dum[j,:] = uk[i,j,:]                  
      
      ymax  = 2.0 * dum.max()
      ymin2 = 0.8 * dum.min()
      ymin  = max(ymin2, 1.e-12*ymax)
   
      hvtitl   = opactitl + ",     T(eV)= %9.3f" %(kt[i])
   
      fig = plt.figure(figsize=(9.0, 6.5))
      plt.title(hvtitl)
      plt.xscale('log')
      plt.yscale('log')
      plt.xlim(0.1,1.e6)
      plt.ylim(ymin, ymax)
      plt.xlabel('Photon Energy (eV)')
      plt.ylabel('Mean Opacity (cm^2/gm)')
      plt.gcf().text(.86, .89, 'gm/cc', fontsize=10, color='black')

      for j in range(0,jrho,1):
         yval  = dum[j,:]
         yhist = np.zeros([2*jhv])
         
         ipt = -1
         for k in range(0,jhv):
           ipt = ipt + 1
           yhist[ipt] = yval[k]
           ipt = ipt + 1
           yhist[ipt] = yval[k]
         
         plt.plot(xhist, yhist, color=colorvar[j], label=rhotabl[j] )

      plt.legend(loc='center left',bbox_to_anchor=(0.90,0.50), fontsize="small")

      plt.grid(True)
      pdf_pages.savefig(fig)
      print("T[ ", ipr, "]  Done")


def plot_uk_compare2(uk, hmid, uk2, hmid2, kt, rhotabl, nhv, nrho, nkt, labl1, labl2):

# Plot uk(nkt,nrho,nhv) versus hv  

   colorvar = cm.rainbow(np.linspace(0,1,nrho))
   dum = np.zeros((nrho,nhv))

   for i in range(0,nkt,1):
      ipr = i + 1
      for j in range(0,nrho,1):                 

         yval  = uk[i,j,:] 
         yval2 = uk2[i,j,:]   
  
         ymax1 = yval.max()
         ymax2 = yval2.max()
         ymax  = 2.0 * max(ymax1, ymax2)
         ymin1 = yval.min()
         ymin2 = yval2.min()
         ymin  = 0.8 * min(ymin1,ymin2)
         ymin  = max(ymin, 1.e-12*ymax)

         jpr = j + 1
   
         hvtitl   = "  T[%2.0f]= %9.3f eV,      RHO[%1.0f] = %5.3e  gm/cc" %(ipr, kt[i], jpr, rhotabl[j])
   
         fig = plt.figure(figsize=(9.0, 6.5))
         plt.title(hvtitl)
         plt.xscale('log')
         plt.yscale('log')
         plt.ylim(ymin, ymax)
         plt.xlabel('hv (eV)')
         plt.ylabel('Opacity (cm^2/gm)')

         plt.plot(hmid,  yval,  color='red',  label=labl1)
         plt.plot(hmid2, yval2, color='blue', label=labl2)

         plt.legend(loc='center left',bbox_to_anchor=(0.75,0.88), fontsize="small")

#         plt.grid(True)
         pdf_pages.savefig(fig)
      print("I =", ipr, "  Done")

def plot_uk_compare2_hist(uk, hnur, uk2, hnur2, kt, rhotabl, nhv, nrho, nkt, labl1, labl2):

# Plot uk(nkt,nrho,nhv) versus hv  

   colorvar = cm.rainbow(np.linspace(0,1,nrho))
   dum = np.zeros((nrho,nhv))

   xhist  = np.zeros([2*nhv])
   xhist2 = np.zeros([2*nhv])
   ipt = -1
   for k in range(0,nhv):
     ipt = ipt + 1
     xhist[ipt]  = hnur[k]
     xhist2[ipt] = hnur2[k]
     ipt = ipt + 1
     xhist[ipt]  = hnur[k+1]
     xhist2[ipt] = hnur2[k+1]

   indx= [3,9,15,18,19]

   # for ij in range(0,5,1):
   for i in range(0,nkt,1):
      ipr = i + 1
      for j in range(0,nrho,1):                 

         yval  = uk[i,j,:] 
         yval2 = uk2[i,j,:]   
  
         ymax1 = yval.max()
         ymax2 = yval2.max()
         ymax  = 2.0 * max(ymax1, ymax2)
         ymin1 = yval.min()
         ymin2 = yval2.min()
         ymin  = 0.8 * min(ymin1,ymin2)
         ymin  = max(ymin, 1.e-12*ymax)

         yhist  = np.zeros([2*nhv])
         yhist2 = np.zeros([2*nhv])
         
         ipt = -1
         for k in range(0,nhv):
           ipt = ipt + 1
           yhist[ipt]  = yval[k]
           yhist2[ipt] = yval2[k]
           ipt = ipt + 1
           yhist[ipt]  = yval[k]
           yhist2[ipt] = yval2[k]

         jpr = j + 1
   
         hvtitl   = "  T[%2.0f]= %9.3f eV,      RHO[%1.0f] = %5.3e  gm/cc" %(ipr, kt[i], jpr, rhotabl[j])
   
         fig = plt.figure(figsize=(9.0, 6.5))
         plt.title(hvtitl)
         plt.xscale('log')
         plt.yscale('log')
         plt.ylim(ymin, ymax)
         plt.xlabel('hv (eV)')
         plt.ylabel('Opacity (cm^2/gm)')

         plt.plot(xhist,  yhist,  color='red',  label=labl1)
         plt.plot(xhist2, yhist2, color='blue', label=labl2)

         plt.legend(loc='center left',bbox_to_anchor=(0.75,0.88), fontsize="small")

#         plt.grid(True)
         pdf_pages.savefig(fig)
      print("I =", ipr, "  Done")

def plotalts(tval, pval, pmax, pmin, ncurves, yldname, alts):

   colorvar = plt.cm.rainbow(np.linspace(0,1,ncurves))

   xmax  = 10
   xmin  = 1.e-6

   ymin = 0.8 * pmin
   ymax = 1.5 * pmax
   
   fig = plt.figure(figsize=(9.0, 6.5))
   plt.title(yldname)
   plt.xscale('log')
   plt.yscale('log')
   plt.xlim(xmin, xmax)
   plt.ylim(ymin, ymax)
   plt.xlabel('T (s)')
   plt.ylabel('Thermal Power (W)')
   xpl  = 0.92
   ypl  = 0.15
   dely = 0.035
   font = FontProperties()
   font.set_family("monospace")
   for i in range(0,ncurves):
     xp = tval[i]
     yp = pval[i]
     plt.plot(xp, yp, color=colorvar[i])
     plt.gcf().text(xpl, ypl, alts[i], fontproperties= font, fontsize=10, color=colorvar[i])
     ypl = ypl + dely

   pdf_pages.savefig(fig)

def plotylds(tval, pval, pmax, pmin, ncurves, altname, ylds):

   colorvar = plt.cm.rainbow(np.linspace(0,1,ncurves))

   xmax  = 10
   xmin  = 1.e-6

   ymin = 0.8 * pmin
   ymax = 1.5 * pmax
   
   fig = plt.figure(figsize=(9.0, 6.5))
   plt.title(altname)
   plt.xscale('log')
   plt.yscale('log')
   plt.xlim(xmin, xmax)
   plt.ylim(ymin, ymax)
   plt.xlabel('T (s)')
   plt.ylabel('Thermal Power (W)')
   xpl  = 0.92
   ypl  = 0.15
   dely = 0.035
   font = FontProperties()
   font.set_family("monospace")
   for i in range(0,ncurves):
     xp = tval[i]
     yp = pval[i]
     plt.plot(xp, yp, color=colorvar[i])
     plt.gcf().text(xpl, ypl, ylds[i], fontproperties= font, fontsize=10, color=colorvar[i])
     ypl = ypl + dely

   pdf_pages.savefig(fig)    

def rdrtvar(fname):

   mydmp = open(fname)
   mydmp.readline()
   sim_titl = mydmp.readline()
   mydmp.readline()
   line   = mydmp.readline()
   values = line.split()
   npts   = int(values[0]) 		# Number of time points

   tpwr   = np.zeros([npts])
   ptherm = np.zeros([npts])
   ke     = np.zeros([npts])
   ie     = np.zeros([npts])
   erad   = np.zeros([npts])
   etot   = np.zeros([npts])
   rshk   = np.zeros([npts])
   ritmx  = np.zeros([npts])
   titmxm = np.zeros([npts])
   titmx  = np.zeros([npts])
   titmxp = np.zeros([npts])
   rhorat = np.zeros([npts])
   flupow = np.zeros([npts])

   mydmp.readline()
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

   return sim_titl, npts, tpwr, ptherm, ke, ie, erad, etot, rshk, ritmx

def plotrtvar(tpwr, ptherm, ke, ie, erad, etot, rshk, ritmx, npts, sim_titl):

   subtitl = 'Energies: Total (blue), Internal (green), Kinetic (red)'

   ymax  = 1.05 * max(etot)
   ymin  = 0.9 * min(ke)
   ymin2 = max(ymin, 1.e-6*ymax)
   xmax  = 1.1 * max(tpwr)
   xmin  = tpwr[0]   
   fig   = plt.figure(figsize=(9.0, 6.5))
   plt.title(sim_titl)
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
   plt.gcf().text(0.20, 0.89, subtitl, fontsize=11, color='black')
   pdf_pages.savefig(fig)

   subtitl = 'Thermal Power With No Averaging'

   ymax  = 1.5*ptherm.max()
   ymin  = 0.8*ptherm.min()
   if(ymin < 1.e-5*ymax): ymin = 1.e-5*ymax
   fig   = plt.figure(figsize=(9.0, 6.5))
   plt.title(sim_titl)
   plt.xscale('log')
   plt.yscale('log')
   plt.ylim(ymin, ymax)
   plt.xlabel('T (s)')
   plt.ylabel('Power (W)')
   plt.plot(tpwr, ptherm,  color='blue')
   plt.grid(True)
   plt.gcf().text(0.20, 0.89, subtitl, fontsize=11, color='black')
   pdf_pages.savefig(fig)
   
   subtitl = 'R Shock in blue,  R Fireball in red'

   ritmx = 0.01 * ritmx		# convert to m
   rshk  = 0.01 * rshk		# convert to m   
   ymax  = 1.05 * max(rshk)
   ymin  = 0.0
   fig   = plt.figure(figsize=(9.0, 6.5))
   plt.title(sim_titl)
   plt.xscale('log')
   plt.yscale('log')
   plt.ylim(0.10, ymax)
   plt.xlabel('T (s)')
   plt.ylabel('Radii (m)')
   plt.plot(tpwr, rshk,  color='blue')
   plt.plot(tpwr, ritmx, color='red')
   plt.grid(True)
   plt.gcf().text(0.20, 0.89, subtitl, fontsize=11, color='black')
   pdf_pages.savefig(fig)
   
def rd_rflo_dmp(fname):

   mydmp  = open(fname,"r")  
   sim_titl = mydmp.readline()
   sim_titl = sim_titl.strip()
   mydmp.readline()

   time   = float(mydmp.readline())		# time in seconds
   ncycle = int(mydmp.readline())		
   ndc    = int(mydmp.readline())		# number of debris cells
   itmx   = int(mydmp.readline())		# Fireball radius index
   irmx   = int(mydmp.readline())		# Shock radius index
   icmax  = int(mydmp.readline())		# number of cells

   rc    = np.zeros((icmax))  			# cell center in m
   dr    = np.zeros((icmax))			# cell width in cm
   rho   = np.zeros((icmax))			# mass density (rho/cc)
   uc    = np.zeros((icmax))			# velocity (m/s)
   pr    = np.zeros((icmax))			# Pressure cgs
   tkel  = np.zeros((icmax))			# T in Kelvin
   mass  = np.zeros((icmax))			# cell mass in kg
   sie   = np.zeros((icmax))			# internal energy erg/gm
   cs    = np.zeros((icmax)) 			# Sound Speed (m/s)  
   gamma = np.zeros((icmax))			# cell effective gamma

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
   return icmax, rc, dr, rho, uc, pr, tkel, mass, sie, cs, gamma, time, ncycle, sim_titl

def plothyd(icmax, rc, rho, tkel, sie, uc, dr, mass, gamma, sim_titl, subtitl, figname, itype, rmax):

   fig, (ax0, ax1) = plt.subplots(nrows=2,figsize=(9.0, 6.5))

   rmin = 0.95*rc[0]
   
#  Plot RHO vs RC and UC vs RC

   fulltitl = sim_titl[0:45] + '  ' + subtitl
 
   ax0.set_title(fulltitl, fontsize = 10)
   ax0.set_xlim(rmin,rmax)
   ax0.set_xscale('log')
   ax0.set_yscale('log')
   ax0.set_ylabel('RHO (gm/cc)', color='blue', fontsize=8)
   ax0.plot(rc,rho)
   
   ax2 = ax0.twinx()
   ax2.set_ylabel('UC (km/s)', color='red', fontsize=8)
   uckm = 0.001 * uc
   ax2.plot(rc, uckm, color='red')
   
# Plot TK vs RC and GAMMA - 1 vs RC

   gmm = gamma - 1.

   ax1.set_xscale('log')
   ax1.set_xlim(rmin,rmax)
   ax1.set_yscale('log')
   ax1.set_xlabel('RC (m)')
   ax1.set_ylabel('T (K)', color='blue', fontsize=8)
   ax1.plot(rc,tkel)
   
   ax2 = ax1.twinx()
   ax2.set_yscale('linear')
   ax2.set_ylim(0., 1.)
   ax2.set_ylabel('GAMMA - 1', color='red', fontsize=8)
   ax2.plot(rc, gmm, color='red')  
   
   fig.tight_layout()

   if(itype > 0):   plt.savefig(figname)		# gif output
   if(itype < 1):   pdf_pages.savefig(fig) 		# pdf output


def rdxdep(fname):
   mydmp = open(fname)

#  xdep.owt is written twice - coarse grid and full grid

#  coarse grid

   mydmp.readline()
   sim_titl = mydmp.readline()
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
   sim_titl = mydmp.readline()
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
   return sim_titl, rc, temp, npts, rc2, temp2, npts2

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


# Main Program Starts Here ***********************************

basename  = '/home/e.symbalisty.adm/RFLOCode/USTESTS/'
basename2 = '/home/e.symbalisty.adm/RFLOCode/USTESTS/'
dataname  = '/home/e.symbalisty.adm/DataFiles/'
dataname2 = '/home/e.symbalisty.adm/DataFiles/USTESTS/'
baseall   = '/home/e.symbalisty.adm/RFLOCode/'
basealts  = '/home/e.symbalisty.adm/RFLOCode/'
basesim   = '/home/e.symbalisty.adm/RFLOCode/BookRuns/KT1KM12/'

eosname   = 'RFLOCHM/magpie_eos.txt'
eosname2  = 'RFLOCHM/eos_3720.txt'

opacname1 = 'RFLOCHM/Hai_78_rosseland.txt'
opacname2 = 'RFLOCHM/LANL_AbsRosseland_78_v1.txt'
opac_pix  = dataname + 'RFLOCHM/LLNL_LANL_Rosseland.pdf' 	# Used in plot_opac78
# opacname2 = 'RFLOCHM/alum_ross_78.txt'
# opacname1 = 'RFLOCHM/alum_planck_78.txt'
# opac_pix  = dataname + 'RFLOCHM/Aluminum_Planck.pdf'		# Used in plot_deb78
	
ciraname = dataname + 'RFLOCHM/cira_100km.out'
eltname  = dataname + 'RFLOCHM/elt_hychem_beta_78.owt'
set_globals(ciraname, eltname)

plot_magpie_eos = False
plot_alum_eos   = False
plot_opac78     = False
plot_deb78      = False
plot_opac_cmp2  = False

plot_all        = False
plot_alts       = False
plot_ylds       = False
plot_power      = False
plot_compare    = False		# US Tests Comparisons
plot_special    = False
plot_sim        = False
plot_2sims      = False
plot_rtvar      = False
plot_hydro_pdf  = False
plot_hydro_gif  = True
plot_xdep       = False
atten = 1			# Elterman clear air attenuation
# atten = 2			# Elterman aerosol attenuation
# atten = 3			# 0.87^Natm clear air attenuation

ncases   = 10
casename = ['ABLE_TS', 'BAKER', 'CHARLIE_BJ', 'CLIMAX', 'DIXIE', 'ENCORE', 'GRABLE', 
            'HA', 'WASP', 'WASP_PRIME']
complabl = ':   LANL AbsRosseland (red) '			# Label for the simulations
simpwr   = '/pwr78_ave_iopac1_abs_rosseland.plt'			# Pick the simulation output file
cmpname  = basename + 'power_lanl_compare_data.pdf'			# Name for the comparison PDF file
sim1     ='pwr78_ave_iopac9.plt'
sim2     ='pwr78_ave_iopac1_abs_rosseland.plt'

if(plot_power):
   pdf_pages = PdfPages(basename + 'power.pdf')
   for icase in range(ncases-1,ncases,1):

     #  read simulation power versus time
     fname  = basename +  casename[icase] + '/pwr78_ave.plt'
     sim_titl, tpwr, pow_band78, totpow78, nsim, nbins = rdpower78(fname)
     pl_titl =  casename[i] + '  ' + sim_titl[0:80]

     plotpwr(tpwr, pow_band78, totpow78, nsim, pl_titl)

     print(casename[icase] + " done")
   pdf_pages.close()

if(plot_sim):
   pdf_pages = PdfPages(basesim + 'power.pdf')

   #  read simulation 
   fname  = basesim + 'pwr78_ave.plt'
   sim_titl, tpwr, pow_band78, totpow78, nsim, nbins = rdpower78(fname)

   plotpwr(tpwr, pow_band78, totpow78, nsim, sim_titl)

   print("plot_sim done")
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
     pdf_pages = PdfPages(baseall + 'PIX_v71/' + 'power_' + ylds[j] + '.pdf')

     for i in range(0,21,1):
        fname  = baseall + ylds[j] + '/' + alts[i] + '/pwr78_ave.plt'
        sim_titl, tpwr, pow_band78, totpow78, nsim, nbins = rdpower78(fname)
        pl_titl = sim_titl[0:34]

        plotpwr(tpwr, pow_band78, totpow78, nsim, pl_titl)

        print(fname + " power plot done")

     pdf_pages.close()
     print(" ")

if(plot_alts):
   owtname = 'PIX_v71/' + 'power_llnl_rp_alts.pdf'
   pdf_pages = PdfPages(basealts + owtname)
   ylds   = ['T0010', 'T0032', 'T0100', 'T0320', 'T1000', 'T3200', 
             'KT0010', 'KT0100', 'MT10']
   alts   = ['KM00','KM04','KM08','KM12','KM16','KM20','KM24', \
             'KM28','KM32','KM36','KM40','KM44','KM48','KM52', \
             'KM56','KM60','KM64','KM68','KM72','KM76','KM80' ]
   nylds = 9
   nalts = 21

   for iyld in range(0,nylds,1):
     yldname = ylds[iyld]
     tval = [ ]
     pval = [ ]
     pmax = 0.
     pmin = 1.e20
     for ialt in range(0,nalts,1):
       fname = basealts + yldname + '/' + alts[ialt] + '/pwr78_ave.plt'
       sim_titl, tpwr, pow_band78, totpow78, nsim, nbins = rdpower78(fname)
       ymax = totpow78.max()
       if(ymax > pmax): pmax = ymax
       tval.append(tpwr)
       pval.append(totpow78)
       for k in range(0,nsim):
          if(tpwr[k] > 1.e-6):
             if(totpow78[k] < pmin): pmin = totpow78[k]

     plotalts(tval, pval, pmax, pmin, nalts, yldname, alts)
     print(yldname + ' done')

   print(owtname)
   print(fname + " plots alts done \n")
   pdf_pages.close()

if(plot_ylds):
   owtname = 'PIX_v71/' + 'power_llnl_rp_ylds.pdf'
   pdf_pages = PdfPages(basealts + owtname)
   ylds   = ['T0010', 'T0032', 'T0100', 'T0320', 'T1000', 'T3200', 
             'KT0010', 'KT0100', 'MT10']
   alts   = ['KM00','KM04','KM08','KM12','KM16','KM20','KM24', \
             'KM28','KM32','KM36','KM40','KM44','KM48','KM52', \
             'KM56','KM60','KM64','KM68','KM72','KM76','KM80' ]
   nylds = 9
   nalts = 21

   for ialt in range(0,nalts,1):
     altname = alts[ialt]
     tval = [ ]
     pval = [ ]
     pmax = 0.
     pmin = 1.e20
     for iyld in range(0,nylds,1):
       fname = basealts + ylds[iyld] + '/' + altname + '/pwr78_ave.plt'
       sim_titl, tpwr, pow_band78, totpow78, nsim, nbins = rdpower78(fname)
       ymax = totpow78.max()
       if(ymax > pmax): pmax = ymax
       tval.append(tpwr)
       pval.append(totpow78)
       for k in range(0,nsim):
          if(tpwr[k] > 1.e-6):
             if(totpow78[k] < pmin): pmin = totpow78[k]

     plotylds(tval, pval, pmax, pmin, nylds, altname, ylds)
     print(altname + ' done')

   print(owtname)
   print(fname + " plots alts done \n")
   pdf_pages.close()

if(plot_2sims):
   pdf_pages = PdfPages(basename + 'power_WASP_PRIME_Compare_2.pdf')

   #  read simulation # 1
   fname  = basename + sim1
   sim_titl, tpwr, pow_band78, totpow78, nsim, nbins = rdpower78(fname)

   #  read simulation # 2
   fname  = basename + sim2
   sim_titl2, tpwr2, pow_band78_2, totpow78_2, nsim2, nbins = rdpower78(fname)

   lably    = 'Thermal Power (W)'
   comptitl = 'Wasp Prime: LLNL R&P (red)  and LANL AbsRosseland (blue) - (100,20,78) binning '
   sim_titl = 'USTESTS/WASP_PRIME'
   plot2sim(tpwr, totpow78, nsim, tpwr2, totpow78_2, nsim2, lably, comptitl, sim_titl)

   print("plot_2sims done")
   pdf_pages.close()

if(plot_special):

   pdf_pages = PdfPages(basename2 + 'power_kt320_compare.pdf')

   #  read simulation in base simulation
   fname  = basename2 + 'KM08/pwr78_ave.plt'
   sim_titl, tpwr, pow_band78, totpow78, nsim, nbins = rdpower78(fname)

   #  read NOCHEM simulation
   fname  = basename2 + 'KM08_NOCHEM/pwr78_ave.plt'
   sim_titl2, tpwr2, pow_band78_2, totpow78_2, nsim2, nbins = rdpower78(fname)

   #  read Air & Aluminum simulation
   fname  = basename2 + 'KM08_AL/pwr78_ave.plt'
   sim_titl3, tpwr3, pow_band78_3, totpow78_3, nsim3, nbins = rdpower78(fname)

   # compare chem vs no chem sims
   lably    = 'Thermal Power (W)'
   comptitl = 'KT320 at 8 KM with CHEM (red) and without CHEM (blue)'
   sim_titl = 'KT0320/KM08 and KT0320/KM08_NOCHEM'
   plot2sim(tpwr, totpow78, nsim, tpwr2, totpow78_2, nsim2, lably, comptitl, sim_titl)

   # compare air vs air & aluminum sims
   lably    = 'Thermal Power (W)'
   comptitl = 'KT320 at 8 KM: Air only (red) and Air + Aluminum (blue)'
   sim_titl = 'KT0320/KM08 and KT0320/KM08_AL'
   plot2sim(tpwr, totpow78, nsim, tpwr3, totpow78_3, nsim3, lably, comptitl, sim_titl)

   print("plot_special done")
   pdf_pages.close()

if(plot_compare):
   pdf_pages = PdfPages(cmpname)
   for icase in range(ncases-1,ncases,1):
     tchk = 0.004		# test for max power after time = tchk in seconds

     if(casename[icase] == 'ABLE_TS'):

       comptitl  = casename[icase] + complabl + ' vs  & Streak Camera Data (blue)'

       #  read streak camera data
       dname = 'Able_TS_13027_exp_vs_time0.txt'
       fname  = dataname2 + casename[icase] + '/' + dname
       tdat, pdat, ndat, datalabl = rd_streak(fname)

       #  estimate number of atmospheres between source and sensor
       hevent  = 1.180		# 793 feet above ground level
       hsensor = 0.938		# 3077 feet 
       srange  = 3.097		# From Spriggs
       roamb   = 1.077		# gm/liter
       pamb    = 888.5		# mbar
       qatm = natm_event_ground(hsensor, hevent, srange) 
       range_labl = 'Slant Range (km) = %6.2f' %srange 
       natm_labl  = 'NATM = %5.2f' %qatm

       #  read simulation
       fname  = basename + casename[icase] + simpwr
       run_titl, tpwr, pow_band78, totpow78, nsim, nbins = rdpower78(fname)
       sim_titl = run_titl[0:15] + ',   rho0 (g/l) = %5.3f' % roamb + ',   p0 (mbar) =%5.1f' %pamb

       #  adjust for atmospheric attenuation
       powadj78, exptau = attenuate78_srange(pow_band78, qatm, hsensor, hevent, srange, nsim, atten)

       #  sum the weighted bands
       ptot = ptot_calc(sxx_w12, powadj78, nsim, nbins)

       #  normalize simulation to power at second maximum
       pnorm = norm_sim(tpwr, ptot, nsim, tchk)

       #  add time offset to data if desired
       toffset, toff_labl, tdat = toff_calc(tpwr, pnorm, nsim, tdat, pdat, ndat, 0.10)

       #  draw the plot  
       lably  = 'Normalized SXX W12 Power'
       plotpwr2(tpwr, pnorm, nsim, tdat, pdat, ndat, lably, comptitl, dname, sim_titl,
          toff_labl, range_labl, natm_labl, ' ', 'SXX W12', sxx_w12)

       # text file outputs
       print_sim_powers(tpwr, totpow78, ptot, nsim, basename+casename[icase], casename[icase], datalabl, fname)

       print(casename[icase] + " done")

     if(casename[icase] == 'BAKER'):

       comptitl  = casename[icase] + complabl + ' vs  & Streak Camera Data (blue)'

       #  read streak camera data
       dname = 'Baker_BJ_10651_exp_vs_time0.txt'
       fname  = dataname2 + casename[icase] + '/' + dname
       tdat, pdat, ndat, datalabl = rd_streak(fname)

       #  estimate number of atmospheres between source and sensor
       hevent  = 1.619		# 1118 feet above ground level
       hsensor = 1.278		# 4193 feet 
       srange  = 3.616		# From Spriggs
       roamb   = 1.020		# gm/liter
       pamb    = 840.		# mbar
       qatm = natm_event_ground(hsensor, hevent, srange) 
       range_labl = 'Slant Range (km) = %6.2f' %srange 
       natm_labl = 'NATM = %5.2f' %qatm

       #  read simulation
       fname  = basename + casename[icase] + simpwr
       run_titl, tpwr, pow_band78, totpow78, nsim, nbins = rdpower78(fname)
       sim_titl = run_titl[0:15] + ',   rho0 (g/l) = %5.3f' % roamb + ',   p0 (mbar) =%5.1f' %pamb

       #  adjust for atmospheric attenuation
       powadj78, exptau = attenuate78_srange(pow_band78, qatm, hsensor, hevent, srange, nsim, atten)

       #  sum the weighted bands
       ptot = ptot_calc(sxx_w47, powadj78, nsim, nbins)

       #  normalize simulation to power at second maximum
       pnorm = norm_sim(tpwr, ptot, nsim, tchk)

       #  add time offset to data if desired
       toffset, toff_labl, tdat = toff_calc(tpwr, pnorm, nsim, tdat, pdat, ndat, 0.10)

       #  draw the plot  
       lably  = 'Normalized SXX W47 Power'
       plotpwr2(tpwr, pnorm, nsim, tdat, pdat, ndat, lably, comptitl, dname, sim_titl,
          toff_labl, range_labl, natm_labl, ' ', 'SXX W47', sxx_w47)

       # text file outputs
       print_sim_powers(tpwr, totpow78, ptot, nsim, basename+casename[icase], casename[icase], datalabl, fname)

       print(casename[icase] + " done")      

     if(casename[icase] == 'CHARLIE_BJ'):

       comptitl  = casename[icase] + complabl + ' vs  & Streak Camera Data (blue)'

       #  read streak camera data
       dname = 'Charlie_BJ_10732_exp_vs_time0.txt'
       fname  = dataname2 + casename[icase] + '/' + dname
       tdat, pdat, ndat, datalabl = rd_streak(fname)

       #  estimate number of atmospheres between source and sensor
       hevent  = 1.623		# 1132 feet above ground level
       hsensor = 1.278		# 4193 feet 
       srange  = 3.626		# From Spriggs
       roamb   = 1.015		# gm/liter
       pamb    = 835.		# mbar
       qatm = natm_event_ground(hsensor, hevent, srange) 
       range_labl = 'Slant Range (km) = %6.2f' %srange 
       natm_labl = 'NATM = %5.2f' %qatm

       #  read simulation
       fname  = basename + casename[icase] + simpwr
       run_titl, tpwr, pow_band78, totpow78, nsim, nbins = rdpower78(fname)
       sim_titl = run_titl[0:15] + ',   rho0 (g/l) = %5.3f' % roamb + ',   p0 (mbar) =%5.1f' %pamb

       #  adjust for atmospheric attenuation
       powadj78, exptau = attenuate78_srange(pow_band78, qatm, hsensor, hevent, srange, nsim, atten)

       #  sum the weighted bands
       ptot = ptot_calc(sxx_w12, powadj78, nsim, nbins)

       #  normalize simulation to power at second maximum
       pnorm = norm_sim(tpwr, ptot, nsim, tchk)

       #  add time offset to data if desired
       toffset, toff_labl, tdat = toff_calc(tpwr, pnorm, nsim, tdat, pdat, ndat, 0.10)

       #  draw the plot  
       lably  = 'Normalized SXX W12 Power'
       plotpwr2(tpwr, pnorm, nsim, tdat, pdat, ndat, lably, comptitl, dname, sim_titl,
          toff_labl, range_labl, natm_labl, ' ', 'SXX W12', sxx_w12)

       # text file outputs
       print_sim_powers(tpwr, totpow78, ptot, nsim, basename+casename[icase], casename[icase], datalabl, fname)

       print(casename[icase] + " done")      

     if(casename[icase] == 'CLIMAX'):

       comptitl  = casename[icase] + complabl + ' vs  & Streak Camera Data (blue)'

       #  read streak camera MF_W12 data
       dname = 'Climax_171028_exp_vs_time0.txt'
       fname = dataname2 + casename[icase] +'/' + dname
       tdat, pdat, ndat, datalabl = rd_streak(fname)

       #  read streak camera MF_NONE data
       dname2 = 'Climax_171027_exp_vs_time0.txt'
       fname  = dataname2 + casename[icase] + '/' + dname2
       tdat2, pdat2, ndat2, datalabl2 = rd_streak(fname)

       #  estimate number of atmospheres between source and sensor
       hevent  = 1.6395		# 1334 feet = 0.406 km  above ground level
       hsensor = 1.2268		# 4025 feet = 1.227
       srange  = 5.176 		# From Spriggs
       roamb   = 1.004		# gm/liter
       pamb    = 824.		# mbar
       qatm = natm_event_ground(hsensor, hevent, srange) 
       range_labl = 'Slant Range (km) = %6.2f' %srange 
       natm_labl = 'NATM = %5.2f' %qatm

       #  read base simulation
       fname  = basename + casename[icase] + simpwr
       run_titl, tpwr, pow_band78, totpow78, nsim, nbins = rdpower78(fname)
       sim_titl = run_titl[0:15] + ',   rho0 (g/l) = %5.3f' % roamb + ',   p0 (mbar) =%5.1f' %pamb

       #  read no chem simulation
       fname  = basename + 'CLIMAX_NOCHEM/pwr78_ave.plt'
       run_titl2, tpwr2, pow_band78_2, totpow78_2, nsim2, nbins = rdpower78(fname)
       sim_titl2 = run_titl2[0:15] + ',   rho0 (g/l) = %5.3f' % roamb + ',   p0 (mbar) =%5.1f' %pamb

       # read air + aluminum simulation
       fname  = basename + 'CLIMAX_AL/pwr78_ave.plt'
       run_titl3, tpwr3, pow_band78_3, totpow78_3, nsim3, nbins = rdpower78(fname)
       sim_titl3 = run_titl3[0:15] + ',   rho0 (g/l) = %5.3f' % roamb + ',   p0 (mbar) =%5.1f' %pamb

       #  adjust for atmospheric attenuation
       powadj78,   exptau = attenuate78_srange(pow_band78,   qatm, hsensor, hevent, srange, nsim,  atten)
       powadj78_2, exptau = attenuate78_srange(pow_band78_2, qatm, hsensor, hevent, srange, nsim2, atten)
       powadj78_3, exptau = attenuate78_srange(pow_band78_3, qatm, hsensor, hevent, srange, nsim3, atten)

       #  sum the weighted bands to get total simulated sensor power for MF12 data
       ptot  = ptot_calc(strk_mf12, powadj78,   nsim, nbins)

       # normalize simulation to power at second maximum 
       pnorm  = norm_sim(tpwr,  ptot,  nsim, tchk)

       # text file outputs
       print_sim_powers(tpwr, totpow78, ptot, nsim, basename+casename[icase], casename[icase]+'_MF12', datalabl, fname)

       #  add time offset to data if desired
       toffset, toff_labl, tdat = toff_calc(tpwr, pnorm, nsim, tdat, pdat, ndat, 0.10)

       #  draw the plot  
       lably  = 'Normalized MF12 Power'
       plotpwr2(tpwr, pnorm, nsim, tdat, pdat, ndat, lably, comptitl, dname, sim_titl,
         toff_labl, range_labl, natm_labl,'  ', 'MF12', strk_mf12)

       #  sum the weighted bands to get total simulated sensor power for MF NONE data
       ptot  = ptot_calc(strk_mfnone, powadj78,   nsim,  nbins)
       ptot2 = ptot_calc(strk_mfnone, powadj78_2, nsim2, nbins)
       ptot3 = ptot_calc(strk_mfnone, powadj78_3, nsim3, nbins)

       #  normalize simulation to power at second maximum
       pnorm  = norm_sim(tpwr,  ptot,  nsim,  tchk)
       pnorm2 = norm_sim(tpwr2, ptot2, nsim2, tchk)

       #  add time offset to data if desired
       toffset, toff_labl, tdat2 = toff_calc(tpwr, pnorm, nsim, tdat2, pdat2, ndat2, 0.10)

       #  draw the plot  
       lably  = 'Normalized MF NONE Power'
       plotpwr2(tpwr, pnorm, nsim, tdat2, pdat2, ndat2, lably, comptitl, dname2, sim_titl,
          toff_labl, range_labl, natm_labl, ' ', 'MF NONE', strk_mfnone)

       # compare chem vs no chem for the MF NONE case
       lably    = 'MFNONE Power (W)'
       comptitl = 'CLIMAX with CHEM (red) and without CHEM (blue)'
       sim_titl = 'CLIMAX and CLIMAX_NOCHEM, 60 KT, 1.64 km ASL'
       plotpwr2(tpwr, ptot, nsim, tpwr2, ptot2, nsim2, lably, comptitl, ' ', sim_titl,
          '  ', range_labl, natm_labl, ' ', 'MF NONE', strk_mfnone)

       # compare air vs air & aluminum for the MF NONE case
       lably    = 'MFNONE Power (W)'
       comptitl = 'CLIMAX: Air only (red) and Air + Aluminum (blue)'
       sim_titl = 'CLIMAX and CLIMAX_AL, 60 KT, 1.64 km ASL'
       plotpwr2(tpwr, ptot, nsim, tpwr3, ptot3, nsim3, lably, comptitl, ' ', sim_titl,
          '  ', range_labl, natm_labl, ' ', 'MF NONE', strk_mfnone)

       # text file outputs
       print_sim_powers(tpwr, totpow78, ptot, nsim, basename+casename[icase], casename[icase]+'_MFNONE', datalabl2, fname)

       print(casename[icase] + " done") 

     if(casename[icase] == 'DIXIE'):

       comptitl  = casename[icase] + complabl + ' vs  & Streak Camera Data (blue)'

       #  read MF_NONE streak camera data
       dname = 'Dixie_17328_exp_vs_time0.txt'
       fname = dataname2 + casename[icase] + '/' + dname
       tdat, pdat, ndat, datalabl = rd_streak(fname)

       #  estimate number of atmospheres between source and sensor
       hevent  = 3.22		# 10211 feet
       hsensor = 1.277		# 4191 feet
       srange  = 7.95		# From Spriggs
       roamb   = 0.8764		# gm/liter
       pamb    = 686.		# mbar
       qatm = natm_event_ground(hsensor, hevent, srange) 
       range_labl = 'Slant Range (km) = %6.2f' %srange 
       natm_labl  = 'NATM = %5.2f' %qatm

       #  read simulation
       fname  = basename + casename[icase] + simpwr
       run_titl, tpwr, pow_band78, totpow78, nsim, nbins = rdpower78(fname)
       sim_titl = run_titl[0:15] + ',   rho0 (g/l) = %5.3f' % roamb + ',   p0 (mbar) =%5.1f' %pamb

       #  adjust for atmospheric attenuation
       powadj78, exptau = attenuate78_srange(pow_band78, qatm, hsensor, hevent, srange, nsim, atten)

       #  sum the weighted bands to get total simulated sensor power for MF NONE data
       ptot = ptot_calc(strk_mfnone, powadj78, nsim, nbins)

       #  normalize simulation to power at second maximum
       pnorm = norm_sim(tpwr, ptot, nsim, tchk)

       #  add time offset to data if desired
       toffset, toff_labl, tdat = toff_calc(tpwr, pnorm, nsim, tdat, pdat, ndat, 0.10)

       #  draw the plot  
       lably  = 'Normalized MF NONE Power'
       plotpwr2(tpwr, pnorm, nsim, tdat, pdat, ndat, lably, comptitl, dname, sim_titl,
          toff_labl, range_labl, natm_labl, ' ', 'MF NONE', strk_mfnone)

       # text file outputs
       print_sim_powers(tpwr, totpow78, ptot, nsim, basename+casename[icase], casename[icase], datalabl, fname)

       print(casename[icase] + " done")   

     if(casename[icase] == 'ENCORE'):

       comptitl  = casename[icase] + complabl + ' vs  & Streak Camera Data (blue)'

       #  read MF_NONE streak camera data
       dname = 'Encore_17727_exp_vs_time0.txt'
       fname = dataname2 + casename[icase] + '/' + dname
       tdat, pdat, ndat, datalabl = rd_streak(fname)

       #  estimate number of atmospheres between source and sensor
       hevent  = 1.678		# 2433 feet above ground
       hsensor = 0.938		# 3077 feet
       srange  = 2.637		# From Spriggs
       roamb   = 1.0212		# gm/liter
       pamb    = 825.		# mbar
       qatm = natm_event_ground(hsensor, hevent, srange) 
       range_labl = 'Slant Range (km) = %6.2f' %srange 
       natm_labl  = 'NATM = %5.2f' %qatm

       #  read simulation
       fname  = basename + casename[icase] + simpwr
       run_titl, tpwr, pow_band78, totpow78, nsim, nbins = rdpower78(fname)
       sim_titl = run_titl[0:15] + ',   rho0 (g/l) = %5.3f' % roamb + ',   p0 (mbar) =%5.1f' %pamb

       #  adjust for atmospheric attenuation
       powadj78, exptau = attenuate78_srange(pow_band78, qatm, hsensor, hevent, srange, nsim, atten)

       #  sum the weighted bands to get total simulated sensor power for MF NONE data
       ptot = ptot_calc(strk_mfnone, powadj78, nsim, nbins)

       #  normalize simulation to power at second maximum
       pnorm = norm_sim(tpwr, ptot, nsim, tchk)

       #  add time offset to data if desired
       toffset, toff_labl, tdat = toff_calc(tpwr, pnorm, nsim, tdat, pdat, ndat, 0.10)

       #  draw the plot  
       lably  = 'Normalized MF NONE Power'
       plotpwr2(tpwr, pnorm, nsim, tdat, pdat, ndat, lably, comptitl, dname, sim_titl,
          toff_labl, range_labl, natm_labl, ' ', 'MF NONE', strk_mfnone)

       # text file outputs
       print_sim_powers(tpwr, totpow78, ptot, nsim, basename+casename[icase], casename[icase], datalabl, fname)

       print(casename[icase] + " done")  

     if(casename[icase] == 'GRABLE'):

       comptitl  = casename[icase] + complabl + ' vs  & Streak Camera Data (blue)'

       #  read MF_NONE streak camera data
       dname = 'Grable_17927_exp_vs_time0.txt'
       fname = dataname2 + casename[icase] + '/' + dname
       tdat, pdat, ndat, datalabl = rd_streak(fname)

       #  estimate number of atmospheres between source and sensor
       hevent  = 1.1005		# 162.5 m above ground
       hsensor = 0.938		# 3077 feet
       srange  = 2.802		# From Spriggs
       roamb   = 1.0742		# gm/liter
       pamb    = 884.		# mbar
       qatm = natm_event_ground(hsensor, hevent, srange) 
       range_labl = 'Slant Range (km) = %6.2f' %srange 
       natm_labl  = 'NATM = %5.2f' %qatm

       #  read simulation
       fname  = basename + casename[icase] + simpwr
       run_titl, tpwr, pow_band78, totpow78, nsim, nbins = rdpower78(fname)
       sim_titl = run_titl[0:15] + ',   rho0 (g/l) = %5.3f' % roamb + ',   p0 (mbar) =%5.1f' %pamb

       #  adjust for atmospheric attenuation
       powadj78, exptau = attenuate78_srange(pow_band78, qatm, hsensor, hevent, srange, nsim, atten)

       #  sum the weighted bands to get total simulated sensor power for MF NONE data
       ptot = ptot_calc(strk_mfnone, powadj78, nsim, nbins)

       #  normalize simulation to power at second maximum
       pnorm = norm_sim(tpwr, ptot, nsim, tchk)

       #  add time offset to data if desired
       toffset, toff_labl, tdat = toff_calc(tpwr, pnorm, nsim, tdat, pdat, ndat, 0.10)

       #  draw the plot  
       lably  = 'Normalized MF NONE Power'
       plotpwr2(tpwr, pnorm, nsim, tdat, pdat, ndat, lably, comptitl, dname, sim_titl,
          toff_labl, range_labl, natm_labl, ' ', 'MF NONE', strk_mfnone)

       # text file outputs
       print_sim_powers(tpwr, totpow78, ptot, nsim, basename+casename[icase], casename[icase], datalabl, fname)

       print(casename[icase] + " done")     
   

     if(casename[icase] == 'HA'):

       comptitl  = casename[icase] + complabl + '& 250-9000 nm Data (blue)'

       #  read HA data
       dname = 'ha.txt'
       fname = dataname2 + casename[icase] + '/' + dname
       tdat, pdat, ndat = rd_ha(fname)

       #  estimate number of atmospheres between source and sensor
       hevent  = 11.16		# 32582 + 4025 = 36607 feet = 11.16 km
       hsensor = 1.2268		# 4025 feet = 1.227
       srange  = 19.33 		# sqrt(16.58^2 + 9.931^2]
       roamb   = 0.343		# gm/liter
       pamb    = 222.		# mbar
       qatm = natm_event_ground(hsensor, hevent, srange) 
       range_labl = 'Slant Range (km) = %6.2f' %srange 
       natm_labl  = 'NATM = %5.2f' %qatm

       #  read simulation
       fname  = basename + casename[icase] + simpwr
       run_titl, tpwr, pow_band78, totpow78, nsim, nbins = rdpower78(fname)
       sim_titl = run_titl[0:15] + ',   rho0 (g/l) = %5.3f' % roamb + ',   p0 (mbar) =%5.1f' %pamb

       #  adjust for atmospheric attenuation
       powadj78, exptau = attenuate78_srange(pow_band78, qatm, hsensor, hevent, srange, nsim, atten)

       #  sum the weighted bands
       ptot = ptot_calc(bands_250_9000, powadj78, nsim, nbins)

       #  match powers at second maximum
       pamp, qamp = t2match(tpwr, ptot, nsim, tdat, pdat, ndat, tchk)
       labl_amp = 'AMP = %5.3f' %qamp

       #  add time offset to data if desired
       toffset, toff_labl, tdat = toff_calc(tpwr, pamp, nsim, tdat, pdat, ndat, 0.03)

       #  draw the plot  
       plotpwr2(tpwr, pamp, nsim, tdat, pdat, ndat,'Thermal Power (W)',comptitl, dname, sim_titl,
         toff_labl, range_labl, natm_labl, labl_amp, '250 - 9000 nm', bands_250_9000)

       print(casename[icase] + " done")

     if(casename[icase] == 'WASP'):

       comptitl  = casename[icase] + complabl + '& 250-9000 nm Data (blue)'

       #  read data
       dname = 'wasp.txt'
       fname = dataname2 + casename[icase] + '/' + dname
       tdat, pdat, ndat = rd_ha(fname)

       #  estimate number of atmospheres between source and sensor
       hevent  = 1.510		# 762 + 4195 = 4957 feet =  1.510 km
       hsensor = 1.279		# 4195 feet 
       srange  = 16.58 		# = 10.3 miles, From Hilbun XCEL file, ? range or slant range
       roamb   = 1.101		# gm/liter
       pamb    = 847.		# mbar
       qatm = natm_event_ground(hsensor, hevent, srange) 
       range_labl = 'Slant Range (km) = %6.2f' %srange 
       natm_labl  = 'NATM = %5.2f' %qatm 

       #  read simulation
       fname  = basename + casename[icase] + simpwr
       run_titl, tpwr, pow_band78, totpow78, nsim, nbins = rdpower78(fname)
       sim_titl = run_titl[0:15] + ',   rho0 (g/l) = %5.3f' % roamb + ',   p0 (mbar) =%5.1f' %pamb

       #  adjust for atmospheric attenuation
       powadj78, exptau = attenuate78_srange(pow_band78, qatm, hsensor, hevent, srange, nsim, atten)

       #  sum the weighted bands
       ptot = ptot_calc(bands_250_9000, powadj78, nsim, nbins)

       #  match powers at second maximum
       pamp, qamp = t2match(tpwr, ptot, nsim, tdat, pdat, ndat, tchk)
       labl_amp = 'AMP = %5.3f' %qamp

       #  add time offset to data if desired
       toffset, toff_labl, tdat = toff_calc(tpwr, pamp, nsim, tdat, pdat, ndat, 10.)

       #  draw the plot  
       plotpwr2(tpwr, pamp, nsim, tdat, pdat, ndat,'Thermal Power (W)', comptitl, dname, sim_titl,
         toff_labl, range_labl, natm_labl, labl_amp, '250 - 9000 nm', bands_250_9000)

       print(casename[icase] + " done")

     if(casename[icase] == 'WASP_PRIME'):

       comptitl  = casename[icase] + complabl + '& 250-9000 nm Data (blue)'

       #  read data
       dname = 'wasp_prime.txt'
       fname = dataname2 + casename[icase] + '/' + dname
       tdat, pdat, ndat = rd_ha(fname)

       #  estimate number of atmospheres between source and sensor
       hevent  = 1.504		# 739 + 4194 = 4933 feet = 1.504 km
       hsensor = 1.278		# 4194 feet = 1.227
       srange  = 16.58 		# = 10.3 miles, From Hilbun XCEL file
       roamb   = 1.032		# gm/liter
       pamb    = 847.		# mbar
       qatm = natm_event_ground(hsensor, hevent, srange) 
       range_labl = 'Slant Range (km) = %6.2f' %srange 
       natm_labl  = 'NATM = %5.2f' %qatm

       #  read simulation
       fname  = basename + casename[icase] + simpwr
       run_titl, tpwr, pow_band78, totpow78, nsim, nbins = rdpower78(fname)
       sim_titl = run_titl[0:15] + ',   rho0 (g/l) = %5.3f' % roamb + ',   p0 (mbar) =%5.1f' %pamb

       #  adjust for atmospheric attenuation
       powadj78, exptau = attenuate78_srange(pow_band78, qatm, hsensor, hevent, srange, nsim, atten)

       #  sum the weighted bands
       ptot = ptot_calc(bands_250_9000, powadj78, nsim, nbins)

       #  match powers at second maximum
       pamp, qamp = t2match(tpwr, ptot, nsim, tdat, pdat, ndat, tchk)
       labl_amp = 'AMP = %5.3f' %qamp

       #  add time offset to data if desired
       toffset, toff_labl, tdat = toff_calc(tpwr, pamp, nsim, tdat, pdat, ndat, 10.)

       #  draw the plot  
       plotpwr2(tpwr, pamp, nsim, tdat, pdat, ndat, '250 - 9000 nm Power (W)', comptitl, dname, sim_titl,
         toff_labl, range_labl, natm_labl, labl_amp, '250 - 9000 nm', bands_250_9000)

       print(casename[icase] + " done")
   
   pdf_pages.close()

if(plot_magpie_eos):
   pdf_pages = PdfPages(dataname + 'RFLOCHM/magpie_eos.pdf')
   fname   = dataname + eosname

   rhotabl = np.zeros([nrho])
   etable  = np.zeros([nkt])
   fp      = np.zeros([nkt,nrho])
   gt      = np.zeros([nkt,nrho])
   tev     = np.zeros([nkt,nrho])
   cs      = np.zeros([nkt,nrho])

   rhotabl, etable, fp, gt, tev = rd_magpie_eos(fname)

   for i in range(0,nkt,1):
      for j in range(0,nrho,1):
         cs[i,j] = 1.e-5 * np.sqrt(fp[i,j] * (1. + fp[i,j]) * etable[i])

   ploteos(etable, rhotabl, fp, tev, cs, "LANL MAGPIE AIR")
   print("eos plot done")
   pdf_pages.close()

if(plot_alum_eos):
   pdf_pages = PdfPages(dataname + 'RFLOCHM/aluminum_eos.pdf')
   fname   = dataname + eosname2

   rhotabl = np.zeros([nrho])
   etable  = np.zeros([nkt])
   fp      = np.zeros([nkt,nrho])
   gt      = np.zeros([nkt,nrho])
   tev     = np.zeros([nkt,nrho])
   cs	   = np.zeros([nkt,nrho])

   rhotabl, etable, fp, gt, cs = rd_alum_eos(fname)

   cs      = 1.e-5 * cs			# cm/s to km/s
   for i in range(0,nkt,1):
      for j in range(0,nrho,1):
         tev[i,j] = etable[i] * gt[i,j]

   ploteos(etable, rhotabl, fp, tev, cs, "LANL SESAME 3720 ALUMINUM")
   print("eos plot done")
   pdf_pages.close()

if(plot_opac78):
   pdf_pages = PdfPages(opac_pix)
   fname = dataname + opacname1

   nkt, nrho, nhv, rhotabl, kt, hnur, uk = rdopac78(fname)

   hmid  = np.zeros([nhv])
   for i in range(0,nhv,1):
     hmid[i] = 0.5*(hnur[i] + hnur[i+1])

   # plot_uk(uk, kt, hmid, rhotabl, 'T (eV)', 'UK (cm^2/gm)', nhv, nrho, nkt, opacname1)
   plot_uk_hist(uk, kt, hnur, rhotabl, 'T (eV)', 'UK (cm^2/gm)', nhv, nrho, nkt, opacname1)

   print("air opacity plot done")
   pdf_pages.close()

if(plot_deb78):
   pdf_pages = PdfPages(opac_pix)
   fname = dataname + opacname1

   jkt  = 78
   jrho = 30
   jhv  = 78
   rhotabl, kt, hnur, uk = rd_deb78(fname, jkt, jrho, jhv)

   plot_uk_hist(uk, kt, hnur, rhotabl, 'T (eV)', 'UK (cm^2/gm)', jhv, jrho, jkt, opacname1)

   print("debris opacity plot done")
   pdf_pages.close()

if(plot_opac_cmp2):
   pdf_pages = PdfPages(opac_pix)
   fname   = dataname + opacname1
   fname2  = dataname + opacname2

   nkt, nrho, nhv, rhotabl, kt, hnur,  uk  = rdopac78(fname)
   nkt, nrho, nhv, rhotabl, kt, hnur2, uk2 = rdopac78(fname2)

   labl1 = "LLNL Rosseland"
   labl2 = "LANL AbsRosseland_78_v1"
   plot_uk_compare2_hist(uk, hnur, uk2, hnur2, kt, rhotabl, nhv, nrho, nkt, labl1, labl2)

   print("plot_opac_cmp2 done")
   pdf_pages.close()

if(plot_rtvar):
   pdf_pages = PdfPages(basesim + 'rtvar.pdf')
   fname  = basesim + 'rtvar.owt'
   
   # read rtvar.owt 
   sim_titl, npts, tpwr, ptherm, ke, ie, erad, etot, rshk, ritmx = rdrtvar(fname)

   # draw pix
   plotrtvar(tpwr, ptherm, ke, ie, erad, etot, rshk, ritmx, npts, sim_titl) 

   print("rtvar plots done")  
   pdf_pages.close()

if(plot_hydro_gif) or (plot_hydro_pdf):
#  create a list of the dmp file names
   dmpnames = [ ]
   
   mydmp  = open(basesim +'hydrodata',"r")
   mydmp.readline()
   sim_titl = mydmp.readline()
   mydmp.readline()
   print(sim_titl)

#  read in the dmp file names
   while True:
      dum = mydmp.readline()
      if len(dum) == 0:
         break
      dmpnames.append(dum[0:10])
   ndumps = len(dmpnames)
   print("ndumps = " + str(ndumps))

if(plot_hydro_pdf):
   #  loop over the dmp files and draw some pix and save as a PDF file
   #  rc is in meters
   pdf_pages = PdfPages(basesim + 'hydro.pdf')
   nstrt = 0
   nend  = ndumps
	
   for i in range(nstrt,nend,1):
     fname = basesim + dmpnames[i]
     mydmp = open(fname,"r")    
     icmax, rc, dr, rho, uc, pr, tkel, mass, sie, cs, gamma, time, ncycle, sim_titl = rd_rflo_dmp(fname)      
     
     subtitl = ' T(s) =' + str(time) + ' , ncycle = ' + str(ncycle)
     figname = ' '
     itype   = 0		# for pdf output
     rmax = rc.max()		# full grid, reset to desired value for zooming in
   
     plothyd(icmax, rc, rho, tkel, sie, uc, dr, mass, gamma, sim_titl, subtitl, figname, itype, rmax)
     
     print(dmpnames[i])

   print("hydro pdf done")      
   pdf_pages.close() 

if(plot_hydro_gif):
#  loop over the dmp files and draw some pix and save as a GIF movie
   movie_name = basesim + 'hydro.gif'
   images = []
   for i in range(0,ndumps,1):
     fname = basesim + dmpnames[i]
     mydmp = open(fname,"r") 
     icmax, rc, dr, rho, uc, pr, tkel, mass, sie, cs, gamma, time, ncycle, sim_titl = rd_rflo_dmp(fname)        

     subtitl = ' T(s) =' + str(time) + ' , ncycle = ' + str(ncycle)     
     figname = 'last.png'
     itype   = 1		# for movie gif output
     rmax    = rc.max()		# full grid, reset to desired value for zooming in
   
     plothyd(icmax, rc, rho, tkel, sie, uc, dr, mass, gamma, sim_titl, subtitl, figname, itype, rmax)

     images.append(imageio.imread(figname))
     
     print(dmpnames[i])

   imageio.mimsave(movie_name, images, duration=0.1, loop=3)

   print("hydro movie done")

if(plot_xdep):
   pdf_pages = PdfPages(basesim + 'xdep.pdf')

   fname  = basesim + 'xdep.owt'
   sim_titl, rc, temp, npts, rc2, temp2, npts2 = rdxdep(fname)

   plotxdep(rc, temp, npts, rc2, temp2, npts2, sim_titl)

   print("xdep plot done")
   pdf_pages.close()      

