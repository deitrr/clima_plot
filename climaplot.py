import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplcol
import vplot
import sys
import os
import subprocess
import pdb
import re
import scipy.signal as sig
from imp import reload
from textwrap import wrap
import pdb
import matplotlib.animation as animation
import ebm_analytical as ebm
import h5py
from scipy import interpolate
from matplotlib import ticker

vplred = '#C91111'
vplorg = '#E09401'
vpllbl = '#13AED5'
vpldbl = '#1321D8'
vplpur = '#642197'

def OLRhm16(T,pCO2):
  tmpk = np.log10(T)
  phi = np.log10(pCO2)
  f = 9.12805643869791438760*(tmpk*tmpk*tmpk*tmpk)+4.58408794768168803557*(tmpk*tmpk*tmpk)*phi- \
      8.47261075643147449910e+01*(tmpk*tmpk*tmpk)+4.35517381112690282752e-01*(tmpk*phi*tmpk*phi)-\
      2.86355036260417961103e+01*(tmpk*tmpk)*phi+2.96626642498045896446e+02*(tmpk*tmpk)-\
      6.01082900358299240806e-02*tmpk*(phi*phi*phi)-2.60414691486954641420*tmpk*(phi*phi)+\
      5.69812976563675661623e+01*tmpk*phi-4.62596100127381816947e+02*tmpk+\
      2.18159373001564722491e-03*(phi*phi*phi*phi)+1.61456772400726950023e-01*(phi*phi*phi)+\
      3.75623788187470086797*(phi*phi)-3.53347289223180354156e+01*phi+\
      2.75011005409836684521e+02
  Int = pow(10.0,f)/1000.
  return Int
 
def dOLRdThm16(T,pCO2):
  tmpk = np.log10(T)
  phi = np.log10(pCO2)
  f = 4*9.12805643869791438760*(tmpk*tmpk*tmpk)+3*4.58408794768168803557*(tmpk*tmpk)*phi- \
      3*8.47261075643147449910e+01*(tmpk*tmpk)+2*4.35517381112690282752e-01*tmpk*(phi*phi)-\
      2*2.86355036260417961103e+01*tmpk*phi+2*2.96626642498045896446e+02*tmpk-\
      6.01082900358299240806e-02*(phi*phi*phi)-2.60414691486954641420*(phi*phi)+\
      5.69812976563675661623e+01*phi-4.62596100127381816947e+02;
  dI = OLRhm16(T,pCO2) * f / (T); 
  return dI
 
def OLRwk97(T,phi): 
  Int = 9.468980 - 7.714727e-5*phi - 2.794778*T - 3.244753e-3*phi*T-3.547406e-4*(phi*phi) \
      + 2.212108e-2*(T*T) + 2.229142e-3*(phi*phi)*T + 3.088497e-5*phi*(T*T) \
      - 2.789815e-5*(phi*T*phi*T) - 3.442973e-3*(phi*phi*phi) - 3.361939e-5*(T*T*T) \
      + 9.173169e-3*(phi*phi*phi)*T - 7.775195e-5*(phi*phi*phi)*(T*T) \
      - 1.679112e-7*phi*(T*T*T) + 6.590999e-8*(phi*phi)*(T*T*T) \
      + 1.528125e-7*(phi*phi*phi)*(T*T*T) - 3.367567e-2*(phi*phi*phi*phi) - 1.631909e-4*(phi*phi*phi*phi)*T \
      + 3.663871e-6*(phi*phi*phi*phi)*(T*T) - 9.255646e-9*(phi*phi*phi*phi)*(T*T*T)
  return Int
  
def StefBoltz(T):
  sigma = 5.67e-8
  return sigma*T**4

def examine_olr(plname, dir='.', show = True):
  dirf = dir+'/SeasonalClimateFiles'
  if not os.path.exists(dirf):
    raise StandardError('SeasonalClimateFiles directory does not exist')
  else:
    check = 0
    for f in subprocess.check_output('echo '+dirf+'/*.DailyInsol.*',shell=True).split():
      f1 = re.split('\.',re.split('/',f.decode('ascii'))[-1])  #split apart output file
    
      if len(f1) == 4:
        timestamp = f1[3]      
      elif len(f1) == 5:
        timestamp = f1[3]+'.'+f1[4]
        
      time0 = np.float(timestamp)  
      
      if time0 == 0:      
        #get system and planet names
        sysname = f1[0]
        plname = f1[1]
        insolf = dirf+'/'+sysname+'.'+plname+'.DailyInsol.'+timestamp
        tempf = dirf+'/'+sysname+'.'+plname+'.SeasonalTemp.'+timestamp
#         icef = dirf+'/'+sysname+'.'+plname+'.SeasonalIceBalance.'+timestamp
        planckbf = dirf+'/'+sysname+'.'+plname+'.PlanckB.'+timestamp
        check = 1
        
    if check == 0:
      raise StandardError('Climate data not found for time %f'%time)
    
    insol = np.loadtxt(insolf,unpack=True)
    temp = np.loadtxt(tempf,unpack=True)
#     ice = np.loadtxt(icef,unpack=True)
    planckb = np.loadtxt(planckbf,unpack=True)
    output = vplot.GetOutput(dir)
    ctmp = 0
    for p in range(len(output.bodies)):
      if output.bodies[p].name == plname:
        body = output.bodies[p]
        ctmp = 1
      else:
        if p == len(output.bodies)-1 and ctmp == 0:
          raise Exception("Planet %s not found in folder %s"%(plname,dir))
  
    lats = body.Latitude[0]
    try:
      obl = body.Obliquity[np.where(body.Time==time)[0]]
    except:
      obl = getattr(output.log.initial,plname).Obliquity
      if obl.unit == 'rad':
        obl *= 180/np.pi
    
    try:
      ecc = body.Eccentricity[np.where(body.Time==time)[0]]
    except:
      ecc = getattr(output.log.initial,plname).Eccentricity
      
    try:
      longp = (body.LongP+body.PrecA)[np.where(body.Time==time)[0]]
    except:
      try:
        longp = getattr(output.log.initial,plname).LongP
      except:
        try:
          longp = (getattr(output.log.initial,plname).LongA+getattr(out.log.initial,plname).ArgP)
          if longp.unit == 'rad':
            longp *= 180/np.pi
          longp = longp%360
        except:
          longp = 0  
    
    f = open(dir+'/'+plname+'.in','r')
    lines = f.readlines()
    f.close()
    #pdb.set_trace()
    for i in range(len(lines)):
      if lines[i].split() != []:
        if lines[i].split()[0] == 'dRotPeriod':
          P = -1*np.float(lines[i].split()[1]) 
        if lines[i].split()[0] == 'dSemi':
          semi = np.float(lines[i].split()[1]) 
          if semi < 0:
            semi *= -1
        if lines[i].split()[0] == 'dpCO2':
          pco2 = np.float(lines[i].split()[1])
  
    fig = plt.figure(figsize=(10,10))
    fig.suptitle('Obl = %f, Ecc = %f, LongP = %f'%(obl,ecc,longp),fontsize=20)
    
    scale = 4*np.shape(insol)[1]/np.shape(temp)[1]
    Bavg = (np.zeros((240,151))+body.PlanckBAvg[0]).T
    plt.subplot(2,2,1)
    c2=plt.contourf(np.arange(np.shape(planckb)[1])*scale,lats,planckb/Bavg,cmap='plasma')
    plt.title(r'Planck Coeff B [W m$^{-2}$ C$^{-1}$] (4 orbits)',fontsize=12)
    plt.colorbar(c2)
    plt.ylim(lats[0],lats[-1])
    plt.ylabel('Latitude (degrees)')
    
    plt.subplot(2,2,2)
    olr = np.reshape(body.FluxOut,(len(body.Time),len(lats)))
    c = plt.contourf(body.Time,lats,olr.T,cmap='plasma')
    plt.ylabel('Latitude')
    plt.title('Outgoing flux [W/m$^2$]')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    plt.colorbar(c)
    
    plt.subplot(2,2,3)
    temp0 = body.TempLat[0]
    ttmp = np.linspace(np.min(temp0),np.max(temp0),100)
    olrwk = OLRhm16(ttmp+273.15,pco2)
    olr0 = body.FluxOut[0]
    plt.plot(temp0,olr0,'k.')
    plt.plot(ttmp,olrwk,'r-')
    plt.xlabel('Temp')
    plt.ylabel('Outgoing flux [W/m$^2$]')

    tmin = body.TempMinLat[0]-273.15
    tmax = body.TempMaxLat[0]-273.15
    
    olrmin = OLRhm16(tmin+273.15,pco2)
    olrmax = OLRhm16(tmax+273.15,pco2)
    
    Aline = OLRhm16(temp0+273.15,pco2)-body.PlanckBAvg[0]*(temp0+273.15)
    olrminL = Aline + body.PlanckBAvg[0]*(tmin+273.15)
    olrmaxL = Aline + body.PlanckBAvg[0]*(tmax+273.15)
    
    pdb.set_trace()
    plt.subplot(2,2,4)
    plt.plot(tmin,(olrmin-olrminL),'bo')
#     plt.plot(tmin,olrminL,'b--')
    plt.plot(tmax,(olrmax-olrmaxL),'ro')
#     plt.plot(tmax,olrmaxL,'r--')
    plt.xlabel('Temp')
    plt.ylabel('Outgoing flux [W/m$^2$]')

    plt.savefig('olr_detail.pdf')
    if show:
      plt.show()
    else:
      plt.close()
    
def seasonal_maps(time, dir = '.', show = True):
  """
  Creates plots of insolation, temperature, and ice balance
  over the course of an orbit (4 orbits for temp)
  
  Parameters
  ----------
  time : float
    The time of the data you want to plot (see SeasonalClimateFiles directory)
    
  Keyword Arguments
  -----------------
  dir : string
    Directory of vplanet simulation (default = '.')
  show : bool
    Show plot in Python (default = True)
  
  Output
  ------
  PDF format plot with name 'surf_seas_<time>.pdf'
  
  """
  dirf = dir+'/SeasonalClimateFiles'
  if not os.path.exists(dirf):
    raise StandardError('SeasonalClimateFiles directory does not exist')
  else:
    check = 0
    for f in subprocess.check_output('echo '+dirf+'/*.DailyInsol.*',shell=True).split():
      f1 = re.split('\.',re.split('/',f.decode('ascii'))[-1])  #split apart output file
    
      if len(f1) == 4:
        timestamp = f1[3]      
      elif len(f1) == 5:
        timestamp = f1[3]+'.'+f1[4]
        
      time0 = np.float(timestamp)  
      
      if time0 == time:      
        #get system and planet names
        sysname = f1[0]
        plname = f1[1]
        insolf = dirf+'/'+sysname+'.'+plname+'.DailyInsol.'+timestamp
        tempf = dirf+'/'+sysname+'.'+plname+'.SeasonalTemp.'+timestamp
        icef = dirf+'/'+sysname+'.'+plname+'.SeasonalIceBalance.'+timestamp
        planckbf = dirf+'/'+sysname+'.'+plname+'.PlanckB.'+timestamp
        check = 1
        
    if check == 0:
      raise StandardError('Climate data not found for time %f'%time)
    
    insol = np.loadtxt(insolf,unpack=True)
    temp = np.loadtxt(tempf,unpack=True)
    ice = np.loadtxt(icef,unpack=True)
    planckb = np.loadtxt(planckbf,unpack=True)
    output = vplot.GetOutput(dir)
    ctmp = 0
    for p in range(len(output.bodies)):
      if output.bodies[p].name == plname:
        body = output.bodies[p]
        ctmp = 1
      else:
        if p == len(output.bodies)-1 and ctmp == 0:
          raise Exception("Planet %s not found in folder %s"%(plname,dir))
        
    lats = body.Latitude[0]
    try:
      obl = body.Obliquity[np.where(body.Time==time)[0]]
    except:
      obl = getattr(output.log.initial,plname).Obliquity
      if obl.unit == 'rad':
        obl *= 180/np.pi
    
    try:
      ecc = body.Eccentricity[np.where(body.Time==time)[0]]
    except:
      ecc = getattr(output.log.initial,plname).Eccentricity
      
    try:
      longp = (body.LongP+body.PrecA)[np.where(body.Time==time)[0]]
    except:
      try:
        longp = getattr(output.log.initial,plname).LongP
      except:
        try:
          longp = (getattr(output.log.initial,plname).LongA+getattr(out.log.initial,plname).ArgP)
          if longp.unit == 'rad':
            longp *= 180/np.pi
          longp = longp%360
        except:
          longp = 0

    fig = plt.figure(figsize=(10,10))
    fig.suptitle('Time = %f, Obl = %f, Ecc = %f, LongP = %f'%(time,obl,ecc,longp),fontsize=20)
    plt.subplot(2,2,1)
    plt.title(r'Insolation [W/m$^2$] (1 orbit)',fontsize=12)
    c1=plt.contourf(np.arange(np.shape(insol)[1]),lats,insol,cmap='plasma')
    plt.colorbar(c1)
    plt.ylim(lats[0],lats[-1])
    plt.ylabel('Latitude (degrees)')
    
    scale = 4*np.shape(insol)[1]/np.shape(temp)[1]
    plt.subplot(2,2,2)
    c2=plt.contourf(np.arange(np.shape(temp)[1])*scale,lats,temp,cmap='plasma')
    plt.title(r'Surface Temp [$^{\circ}$C] (1 orbit)',fontsize=12)
    plt.colorbar(c2)
    plt.ylim(lats[0],lats[-1])
    plt.ylabel('Latitude (degrees)')
    plt.xlim(0,np.shape(temp)[1]*scale/4.)
    
    scale = np.shape(insol)[1]/np.shape(ice)[1]
    plt.subplot(2,2,3)
    c3=plt.contourf(np.arange(np.shape(ice)[1])*scale,lats,ice,cmap='Blues_r')
    plt.title(r'Ice balance [kg/m$^2$/s] (1 orbit)',fontsize=12)
    plt.colorbar(c3)
    plt.ylim(lats[0],lats[-1])
    plt.ylabel('Latitude (degrees)')
    
    scale = 4*np.shape(insol)[1]/np.shape(temp)[1]
    plt.subplot(2,2,4)
    c2=plt.contourf(np.arange(np.shape(planckb)[1])*scale,lats,planckb,cmap='plasma')
    plt.title(r'Planck Coeff B [W m$^{-2}$ C$^{-1}$] (4 orbits)',fontsize=12)
    plt.colorbar(c2)
    plt.ylim(lats[0],lats[-1])
    plt.ylabel('Latitude (degrees)')
    
    plt.savefig('surf_seas_%.0f.pdf'%time)
    if show:
      plt.show()
    else:
      plt.close()

def clim_spec(plname,dir='.',xrange=False,orbit=True,show=True):
  """
  Creates plots of insolation, temperature, albedo, ice mass,
  and bed rock height over the length of the simulation
  
  Parameters
  ----------
  plname : string
    The name of the planet with .Climate data
    
  Keyword Arguments
  -----------------
  dir : string
    Directory of vplanet simulation (default = '.')
  xrange : float tuple, list, or numpy array
    Range of x-values (time) to restrict plot
    (default = False (no restriction))
  orbit : bool
    Plot orbital data (obliquity, eccentricity, COPP)
    (default = False)
  show : bool
    Show plot in Python (default = True)
  
  Output
  ------
  PDF format plot with name 'evol_<dir>.pdf'
  
  """
  if not isinstance(dir,(list,tuple)):
    dir = [dir]
  
  nfiles = len(dir)

  if nfiles > 1 and orbit == True:
    raise Exception("Error: cannot plot multiple files when orbit = True")
  
  if orbit == True:
    fig = plt.figure(figsize=(15,7))
  else:
    fig = plt.figure(figsize=(10*nfiles,15))

  fig.subplots_adjust(wspace=0.4,top=0.9,right=0.98)

  hdffile = 0
  if os.path.exists('climatedata.hdf5'):
    hdff = h5py.File('climatedata.hdf5')
    hdffile = 1

  for ii in np.arange(nfiles):
    if hdffile:
      time = hdff[dir[ii]]['Time'][...]
      try:
        ecc = hdff[dir[ii]]['Eccentricity'][...]
      except:
        ecc = np.zeros_like(time)+hdff[dir[ii]]['EccInit'][...]
        
      try:
        inc = (hdff[dir[ii]]['Inc'][...])
      except:
        inc = np.zeros_like(time)
        
      try:  
        obl = (hdff[dir[ii]]['Obliquity'][...])
      except:
        obl = np.zeros_like(time)+hdff[dir[ii]]['OblInit'][...]
        
      preca = hdff[dir[ii]]['PrecA'][...]
      try:
        argp = (hdff[dir[ii]]['ArgP'][...])
        longa = hdff[dir[ii]]['LongA'][...]
        longp = (argp+longa+preca)*np.pi/180.0
      except:
        longp = preca*np.pi/180.0
        
      lats = np.unique(hdff[dir[ii]]['Latitude'][...])
      TempLat = hdff[dir[ii]]['TempLat'][...]
      AlbedoLat = hdff[dir[ii]]['AlbedoLat'][...]
      IceHeight = hdff[dir[ii]]['IceHeight'][...]
      BedrockH = hdff[dir[ii]]['BedrockH'][...]
      AnnInsol = hdff[dir[ii]]['AnnInsol'][...]
      
    else:
      out = vplot.GetOutput(dir[ii])
    
    #pdb.set_trace()
  
      ctmp = 0
      for p in range(len(out.bodies)):
        if out.bodies[p].name == plname:
          body = out.bodies[p]
          ctmp = 1
        else:
          if p == len(out.bodies)-1 and ctmp == 0:
            raise Exception("Planet %s not found in folder %s"%(plname,dir[ii]))
  
      try:
        ecc = body.Eccentricity
      except:
        ecc = np.zeros_like(body.Time)+getattr(out.log.initial,plname).Eccentricity
    
      try:
        inc = body.Inc
      except:
        inc = np.zeros_like(body.Time)
    
      try:
        obl = body.Obliquity
      except:
        obltmp = getattr(out.log.initial,plname).Obliquity
        if obltmp.unit == 'rad':
          obltmp *= 180/np.pi
        obl = np.zeros_like(body.Time)+obltmp

      try:
        longp = (body.ArgP + body.LongA + body.PrecA)*np.pi/180.0
      except:
        longp = body.PrecA*np.pi/180.0
    
      lats = np.unique(body.Latitude)
      time = body.Time
      TempLat = body.TempLat
      AlbedoLat = body.AlbedoLat
      IceHeight = body.IceHeight
      BedrockH = body.BedrockH
      AnnInsol = body.AnnInsol
    
    f = open(dir[ii]+'/'+plname+'.in','r')
    lines = f.readlines()
    f.close()
    pco2 = 0
    #pdb.set_trace()
    for i in range(len(lines)):
      if lines[i].split() != []:
        if lines[i].split()[0] == 'dRotPeriod':
          P = -1*np.float(lines[i].split()[1]) 
        if lines[i].split()[0] == 'dSemi':
          semi = np.float(lines[i].split()[1]) 
          if semi < 0:
            semi *= -1
        if lines[i].split()[0] == 'dpCO2':
          pco2 = np.float(lines[i].split()[1])
  
    esinv = ecc*np.sin(longp)*np.sin(obl*np.pi/180.)
    
    titlestr = []
    titlestr.append(r'$a = %f, pCO_2 = %f$'%(semi,pco2))
    titlestr.append(r'$e_0 = %f, i_0 = %f^{\circ}, \psi_0 = %f^{\circ}, P_{rot} = %f$ d'%(ecc[0],inc[0],obl[0],P))
    fig.subplots_adjust(wspace=0.3)

    nlats = len(lats)
    ntimes = len(time)
    
    # plot temperature
    temp = np.reshape(TempLat,(ntimes,nlats))
    if orbit == True:
      ax1 = plt.subplot(2,2,1)
    else:
      ax1 = plt.subplot(5,nfiles,ii+1)
    pos = ax1.figbox.get_points()
    c = plt.contourf(time/1e6,lats,temp.T,cmap='plasma')
    plt.ylabel('Latitude ($^{\circ}$)',fontsize=24,fontweight='bold')
    plt.title(r'Surface Temp ($^{\circ}$C)')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60],fontsize=22,fontweight='bold')
    plt.xticks([0,0.5,1,1.5,2],fontsize=22,fontweight='bold')
    if xrange == False:
      left = 0
    else:
      left = xrange[0]
#     plt.text(left,140,'\n'.join(titlestr),fontsize=20) 
    if xrange:
      plt.xlim(xrange)
    clb=plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
    clb.ax.set_yticklabels(['%d'%val for val in np.linspace(-50,40,10)],weight='bold',fontsize=18)

    # plot albedo
    # alb = np.reshape(AlbedoLat,(ntimes,nlats))
#     if orbit == True:
#       ax2 = plt.subplot(4,2,3)
#     else:
#       ax2 = plt.subplot(5,nfiles,ii+2*nfiles+1)
#     pos = ax2.figbox.get_points()
#     c = plt.contourf(time,lats,alb.T,cmap='Blues_r',rasterized=True)
#     plt.ylabel('Latitude')
#     plt.title('Albedo (TOA)')
#     plt.ylim(-90,90)
#     plt.yticks([-60,-30,0,30,60])
#     if xrange:
#       plt.xlim(xrange)
#     plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
#   

    # plot ice height
    ice = np.reshape(IceHeight,(ntimes,nlats))
    if orbit == True:
      ax3 = plt.subplot(2,2,3)
    else:
      ax3 = plt.subplot(5,nfiles,ii+3*nfiles+1)
    pos = ax3.figbox.get_points()
    c = plt.contourf(time/1e6,lats,ice.T/1000,cmap='Blues_r')
    plt.ylabel('Latitude ($^{\circ}$)',fontsize=24,fontweight='bold')
    plt.title('Ice sheet height (km)')
    plt.ylim(-90,90)
  #   plt.xlim(0,2e6)
    plt.yticks([-60,-30,0,30,60],fontsize=22,fontweight='bold')
    plt.xticks([0,0.5,1,1.5,2],fontsize=22,fontweight='bold')
    plt.xlabel('Time (Myr)',fontsize=24,fontweight='bold')
    if xrange:
      plt.xlim(xrange)
    clb=plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
    clb.ax.set_yticklabels(np.linspace(0,3.5,8),weight='bold',fontsize=18)
    # ax3p = ax3.twinx()
  #   plt.plot(body.Time,esinv,linestyle = 'solid',marker='None',color='salmon',linewidth=2)
  

    # plot bedrock
   #  brock = np.reshape(BedrockH,(ntimes,nlats))
#     if orbit == True:
#       ax4 = plt.subplot(4,2,7)
#     else:
#       ax4 = plt.subplot(5,nfiles,ii+4*nfiles+1)
#     pos = ax4.figbox.get_points()
#     c = plt.contourf(time,lats,brock.T,cmap='Reds_r',rasterized=True)
#     plt.ylabel('Latitude')
#     plt.title('Bedrock height [m]')
#     plt.ylim(-90,90)
#     plt.yticks([-60,-30,0,30,60])
#     plt.xlabel('Time [years]')
#     if xrange:
#       plt.xlim(xrange)
#     plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
#   

    # plot insolation
   #  insol = np.reshape(AnnInsol,(ntimes,nlats))
#     if orbit == True:
#       ax5 = plt.subplot(4,2,2)
#     else:
#       ax5 = plt.subplot(5,nfiles,ii+nfiles+1)
#     pos = ax5.figbox.get_points()
#     c = plt.contourf(time,lats,insol.T,cmap='plasma',rasterized=True)
#     plt.ylabel('Latitude')
#     plt.title(r'Annual average insolation [W m$^{-2}$]')
#     plt.ylim(-90,90)
#     plt.yticks([-60,-30,0,30,60])
#     if xrange:
#       plt.xlim(xrange)
#     plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))

    if orbit == True:
      #obliquity
      plt.subplot(2,2,2)
      plt.plot(time/1e6,obl,linestyle = 'solid',marker='None',color='darkblue',linewidth =2)
      plt.ylabel('Obliquity ($^{\circ}$)',fontsize=24,fontweight='bold')
      plt.xticks([0,0.5,1,1.5,2],fontsize=22,fontweight='bold')
      plt.yticks(fontsize=22,fontweight='bold')
      if xrange:
        plt.xlim(xrange)

      #eccentricity
      plt.subplot(2,2,4)
      plt.plot(time/1e6,ecc,linestyle = 'solid',marker='None',color='darkorchid',linewidth =2)
      plt.ylabel('Eccentricity',fontsize=24,fontweight='bold')
      plt.xlabel('Time (Myr)',fontsize=24,fontweight='bold')
      plt.xticks([0,0.5,1,1.5,2],fontsize=22,fontweight='bold')
      plt.yticks(fontsize=20,fontweight='bold')
      if xrange:
        plt.xlim(xrange)
# 
#       #e sin(obl) sin varpi
#       plt.subplot(4,2,8)
#       plt.plot(time,esinv,linestyle = 'solid',marker='None',color='salmon',linewidth=2,rasterized=True)
#       plt.ylabel('COPP')
#       plt.xlabel('Time [years]')
#       if xrange:
#         plt.xlim(xrange)

    if dir[ii] == '.':
      dir[ii] = 'cwd'
  
  #fig.suptitle('\n'.join(titlestr),fontsize=20) 
  
  if xrange:
    sfile = 'spec_'+'_'.join(dir)+'_%d_%d.pdf'%(xrange[0],xrange[1])
  else:
    sfile = 'spec_'+'_'.join(dir)+'.pdf'
  plt.savefig(sfile)

  if show:
    plt.show()
  else:
    plt.close()

def clim_evol(plname,dir='.',xrange=False,orbit=False,show=True):
  """
  Creates plots of insolation, temperature, albedo, ice mass,
  and bed rock height over the length of the simulation
  
  Parameters
  ----------
  plname : string
    The name of the planet with .Climate data
    
  Keyword Arguments
  -----------------
  dir : string
    Directory of vplanet simulation (default = '.')
  xrange : float tuple, list, or numpy array
    Range of x-values (time) to restrict plot
    (default = False (no restriction))
  orbit : bool
    Plot orbital data (obliquity, eccentricity, COPP)
    (default = False)
  show : bool
    Show plot in Python (default = True)
  
  Output
  ------
  PDF format plot with name 'evol_<dir>.pdf'
  
  """
  if not isinstance(dir,(list,tuple)):
    dir = [dir]
  
  nfiles = len(dir)

  if nfiles > 1 and orbit == True:
    raise Exception("Error: cannot plot multiple files when orbit = True")
  
  if orbit == True:
    fig = plt.figure(figsize=(16,13))
  else:
    fig = plt.figure(figsize=(10*nfiles,15))

  fig.subplots_adjust(wspace=0.3,top=0.9)

  hdffile = 0
  if os.path.exists('climatedata.hdf5'):
    hdff = h5py.File('climatedata.hdf5')
    hdffile = 1

  for ii in np.arange(nfiles):
    if hdffile:
      time = hdff[dir[ii]]['Time'][...]
      try:
        ecc = hdff[dir[ii]]['Eccentricity'][...]
      except:
        ecc = np.zeros_like(time)+hdff[dir[ii]]['EccInit'][...]
        
      try:
        inc = (hdff[dir[ii]]['Inc'][...])
      except:
        inc = np.zeros_like(time)
        
      try:  
        obl = (hdff[dir[ii]]['Obliquity'][...])
      except:
        obl = np.zeros_like(time)+hdff[dir[ii]]['OblInit'][...]
        
      preca = hdff[dir[ii]]['PrecA'][...]
      try:
        argp = (hdff[dir[ii]]['ArgP'][...])
        longa = hdff[dir[ii]]['LongA'][...]
        longp = (argp+longa+preca)*np.pi/180.0
      except:
        longp = preca*np.pi/180.0
        
      lats = np.unique(hdff[dir[ii]]['Latitude'][...])
      TempLat = hdff[dir[ii]]['TempLat'][...]
      AlbedoLat = hdff[dir[ii]]['AlbedoLat'][...]
      IceHeight = hdff[dir[ii]]['IceHeight'][...]
      BedrockH = hdff[dir[ii]]['BedrockH'][...]
      AnnInsol = hdff[dir[ii]]['AnnInsol'][...]
      
    else:
      out = vplot.GetOutput(dir[ii])
    
    #pdb.set_trace()
  
      ctmp = 0
      for p in range(len(out.bodies)):
        if out.bodies[p].name == plname:
          body = out.bodies[p]
          ctmp = 1
        else:
          if p == len(out.bodies)-1 and ctmp == 0:
            raise Exception("Planet %s not found in folder %s"%(plname,dir[ii]))
  
      try:
        ecc = body.Eccentricity
      except:
        ecc = np.zeros_like(body.Time)+getattr(out.log.initial,plname).Eccentricity
    
      try:
        inc = body.Inc
      except:
        inc = np.zeros_like(body.Time)
    
      try:
        obl = body.Obliquity
      except:
        obltmp = getattr(out.log.initial,plname).Obliquity
        if obltmp.unit == 'rad':
          obltmp *= 180/np.pi
        obl = np.zeros_like(body.Time)+obltmp
    
      try:
        longp = (body.ArgP + body.LongA + body.PrecA)*np.pi/180.0
      except:
        longp = body.PrecA*np.pi/180.0
    
      lats = np.unique(body.Latitude)
      time = body.Time
      TempLat = body.TempLat
      AlbedoLat = body.AlbedoLat
      IceHeight = body.IceHeight
      BedrockH = body.BedrockH
      AnnInsol = body.AnnInsol
    
    f = open(dir[ii]+'/'+plname+'.in','r')
    lines = f.readlines()
    f.close()
    pco2 = 0
    #pdb.set_trace()
    for i in range(len(lines)):
      if lines[i].split() != []:
        if lines[i].split()[0] == 'dRotPeriod':
          P = -1*np.float(lines[i].split()[1]) 
        if lines[i].split()[0] == 'dSemi':
          semi = np.float(lines[i].split()[1]) 
          if semi < 0:
            semi *= -1
        if lines[i].split()[0] == 'dpCO2':
          pco2 = np.float(lines[i].split()[1])
  
    esinv = ecc*np.sin(longp)*np.sin(obl*np.pi/180.)
    
    titlestr = []
    titlestr.append(r'$a = %f, pCO_2 = %f$'%(semi,pco2))
    titlestr.append(r'$e_0 = %f, i_0 = %f^{\circ}, \psi_0 = %f^{\circ}, P_{rot} = %f$ d'%(ecc[0],inc[0],obl[0],P))
    fig.subplots_adjust(wspace=0.3)

    nlats = len(lats)
    ntimes = len(time)
    
    if time[-1] >= 1e6:
      time /= 1e6
      tunit = ' (Myr)'
      if xrange:
        xrange /= 1e6
    else:
      tunit = ' (years)'
      
    # plot temperature
    temp = np.reshape(TempLat,(ntimes,nlats))
    if orbit == True:
      ax1 = plt.subplot(4,2,1)
    else:
      ax1 = plt.subplot(5,nfiles,ii+1)
    pos = ax1.figbox.get_points()
    c = plt.contourf(time,lats,temp.T,cmap='plasma',rasterized=True)
    plt.ylabel('Latitude ($^{\circ}$)')
    plt.title(r'Surface Temp ($^{\circ}$C)')
    plt.ylim(np.min(lats),np.max(lats))
    plt.yticks([-60,-30,0,30,60])
    if xrange == False:
      left = 0
    else:
      left = xrange[0]
#     plt.text(left,140,'\n'.join(titlestr),fontsize=20) 
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    # plot albedo
    alb = np.reshape(AlbedoLat,(ntimes,nlats))
    if orbit == True:
      ax2 = plt.subplot(4,2,3)
    else:
      ax2 = plt.subplot(5,nfiles,ii+2*nfiles+1)
    pos = ax2.figbox.get_points()
    c = plt.contourf(time,lats,alb.T,cmap='Blues_r',rasterized=True)
    plt.ylabel('Latitude ($^{\circ}$)')
    plt.title('Albedo')
    plt.ylim(np.min(lats),np.max(lats))
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  

    # plot ice height
    ice = np.reshape(IceHeight,(ntimes,nlats))
    if orbit == True:
      ax3 = plt.subplot(4,2,5)
    else:
      ax3 = plt.subplot(5,nfiles,ii+3*nfiles+1)
    pos = ax3.figbox.get_points()
    c = plt.contourf(time,lats,ice.T,cmap='Blues_r',rasterized=True)
    plt.ylabel('Latitude ($^{\circ}$)')
    plt.title('Ice sheet height (m)')
    plt.ylim(np.min(lats),np.max(lats))
  #   plt.xlim(0,2e6)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
    # ax3p = ax3.twinx()
  #   plt.plot(body.Time,esinv,linestyle = 'solid',marker='None',color='salmon',linewidth=2)
  

    # plot bedrock
    brock = np.reshape(BedrockH,(ntimes,nlats))
    if orbit == True:
      ax4 = plt.subplot(4,2,7)
    else:
      ax4 = plt.subplot(5,nfiles,ii+4*nfiles+1)
    pos = ax4.figbox.get_points()
    c = plt.contourf(time,lats,brock.T,cmap='Reds_r',rasterized=True)
    plt.ylabel('Latitude ($^{\circ}$)')
    plt.title('Bedrock height (m)')
    plt.ylim(np.min(lats),np.max(lats))
    plt.yticks([-60,-30,0,30,60])
    plt.xlabel('Time'+tunit)
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  

    # plot insolation
    insol = np.reshape(AnnInsol,(ntimes,nlats))
    if orbit == True:
      ax5 = plt.subplot(4,2,2)
    else:
      ax5 = plt.subplot(5,nfiles,ii+nfiles+1)
    pos = ax5.figbox.get_points()
    c = plt.contourf(time,lats,insol.T,cmap='plasma',rasterized=True)
    plt.ylabel('Latitude ($^{\circ}$)')
    plt.title(r'Annual average insolation (W m$^{-2}$)')
    plt.ylim(np.min(lats),np.max(lats))
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))

    if orbit == True:
      #obliquity
      plt.subplot(4,2,4)
      plt.plot(time,obl,linestyle = 'solid',marker='None',color='darkblue',linewidth =2,rasterized=True)
      plt.ylabel('Obliquity ($^{\circ}$)')
      if xrange:
        plt.xlim(xrange)

      #eccentricity
      plt.subplot(4,2,6)
      plt.plot(time,ecc,linestyle = 'solid',marker='None',color='darkorchid',linewidth =2,rasterized=True)
      plt.ylabel('Eccentricity')
      if xrange:
        plt.xlim(xrange)

      #e sin(obl) sin varpi
      plt.subplot(4,2,8)
      plt.plot(time,esinv,linestyle = 'solid',marker='None',color='salmon',linewidth=2,rasterized=True)
      plt.ylabel('COPP')
      plt.xlabel('Time'+tunit)
      if xrange:
        plt.xlim(xrange)

    if dir[ii] == '.':
      dir[ii] = 'cwd'
  
  #fig.suptitle('\n'.join(titlestr),fontsize=20) 
  
  if xrange:
    sfile = 'evol_'+'_'.join(dir)+'_%d_%d.pdf'%(xrange[0],xrange[1])
  else:
    sfile = 'evol_'+'_'.join(dir)+'.pdf'
  plt.savefig(sfile,dpi =300)

  if show:
    plt.show()
  else:
    plt.close()

def albLW(plname,dir='.',xrange=False,orbit=False,show=True):
  """
  Creates plots of average temperature, min temp, and max temp
  over the length of the simulation
  
  Parameters
  ----------
  plname : string
    The name of the planet with .Climate data
    
  Keyword Arguments
  -----------------
  dir : string
    Directory of vplanet simulation (default = '.')
  xrange : float tuple, list, or numpy array
    Range of x-values (time) to restrict plot
    (default = False (no restriction))
  orbit : bool
    Plot orbital data (obliquity, eccentricity, COPP)
    (default = False)
  show : bool
    Show plot in Python (default = True)
  
  Output
  ------
  PDF format plot with name 'evol_<dir>.pdf'
  
  """
  if not isinstance(dir,(list,tuple)):
    dir = [dir]
  
  nfiles = len(dir)

  if nfiles > 1 and orbit == True:
    raise Exception("Error: cannot plot multiple files when orbit = True")
  
  titlestr = []
  
  if orbit == True:
    fig = plt.figure(figsize=(16,12))
  else:
    fig = plt.figure(figsize=(10*nfiles,12))

  for ii in np.arange(nfiles):
    out = vplot.GetOutput(dir[ii])
  
    for p in range(len(out.bodies)):
      if out.bodies[p].name == plname:
        body = out.bodies[p]
        ctmp = 1
      else:
        if p == len(out.bodies)-1 and ctmp == 0:
          raise Exception("Planet %s not found in folder %s"%(plname,folders[j]))
  
    try:
      ecc = body.Eccentricity
    except:
      ecc = np.zeros_like(body.Time)+getattr(out.log.initial,plname).Eccentricity
    
    try:
      inc = body.Inc
    except:
      inc = np.zeros_like(body.Time)
    
    try:
      obl = body.Obliquity
    except:
      obltmp = getattr(out.log.initial,plname).Obliquity
      if obltmp.unit == 'rad':
        obltmp *= 180/np.pi
      obl = np.zeros_like(body.Time)+obltmp

    f = open(dir[ii]+'/'+plname+'.in','r')
    lines = f.readlines()
    f.close()
    for i in range(len(lines)):
      if lines[i].split() != []:
        if lines[i].split()[0] == 'dRotPeriod':
          P = -1*np.float(lines[i].split()[1])  

    try:
      longp = (body.ArgP + body.LongA + body.PrecA)*np.pi/180.0
    except:
      longp = body.PrecA*np.pi/180.0
    
    esinv = ecc*np.sin(longp)*np.sin(obl*np.pi/180.)

  
    titlestr.append(r'$e_0 = %f, i_0 = %f, \psi_0 = %f, P_{rot} = %f$'%(ecc[0],inc[0],obl[0],P)) 
    fig.subplots_adjust(wspace=0.3)

    lats = np.unique(body.Latitude)
    nlats = len(lats)
    ntimes = len(body.Time)

    # plot temperature
    maxscl = np.max(body.AlbedoLandLat)
    minscl = np.min(body.AlbedoLandLat)
    norm = mplcol.Normalize(vmin=minscl,vmax=maxscl)
  
  
    temp = np.reshape(body.AlbedoLat,(ntimes,nlats))
    if orbit == True:
      ax1 = plt.subplot(3,2,1)
    else:
      ax1 = plt.subplot(3,nfiles,ii+1)
    pos = ax1.figbox.get_points()
    c = plt.contourf(body.Time,lats,temp.T,cmap='Blues_r',norm = norm)
    plt.ylabel('Latitude')
    plt.title('Albedo')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    tempmin = np.reshape(body.AlbedoLandLat,(ntimes,nlats))
    if orbit == True:
      ax3 = plt.subplot(3,2,3)
    else:
      ax3 = plt.subplot(3,nfiles,ii+nfiles+1)
    pos = ax3.figbox.get_points()
    c = plt.contourf(body.Time,lats,tempmin.T,cmap='Blues_r',norm = norm)
    plt.ylabel('Latitude')
    plt.title('Land Albedo')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    tempmax = np.reshape(body.AlbedoWaterLat,(ntimes,nlats))
    if orbit == True:
      ax5 = plt.subplot(3,2,5)
    else:
      ax5 = plt.subplot(3,nfiles,ii+2*nfiles+1)
    pos = ax5.figbox.get_points()
    c = plt.contourf(body.Time,lats,tempmax.T,cmap='Blues_r',norm = norm)
    plt.ylabel('Latitude')
    plt.title('Water albedo')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    if orbit == True:
      #obliquity
      plt.subplot(3,2,2)
      plt.plot(body.Time,obl,linestyle = 'solid',marker='None',color='darkblue',linewidth =2)
      plt.ylabel('Obliquity')
      if xrange:
        plt.xlim(xrange)

      #eccentricity
      plt.subplot(3,2,4)
      plt.plot(body.Time,ecc,linestyle = 'solid',marker='None',color='darkorchid',linewidth =2)
      plt.ylabel('Eccentricity')
      if xrange:
        plt.xlim(xrange)

      #e sin(obl) sin varpi
      plt.subplot(3,2,6)
      plt.plot(body.Time,esinv,linestyle = 'solid',marker='None',color='salmon',linewidth=2)
      plt.ylabel('COPP')
      plt.xlabel('Time [years]')
      if xrange:
        plt.xlim(xrange)

    if dir[ii] == '.':
      dir[ii] = 'cwd'

  fig.suptitle('\n'.join(titlestr),fontsize=20) 

  if xrange:
    sfile = 'tempLW_'+'_'.join(dir)+'_%d_%d.pdf'%(xrange[0],xrange[1])
  else:
    sfile = 'tempLW_'+'_'.join(dir)+'.pdf'
  plt.savefig(sfile)
  if show:
    plt.show()
  else:
    plt.close()


def tempLW(plname,dir='.',xrange=False,orbit=False,show=True):
  """
  Creates plots of average temperature, min temp, and max temp
  over the length of the simulation
  
  Parameters
  ----------
  plname : string
    The name of the planet with .Climate data
    
  Keyword Arguments
  -----------------
  dir : string
    Directory of vplanet simulation (default = '.')
  xrange : float tuple, list, or numpy array
    Range of x-values (time) to restrict plot
    (default = False (no restriction))
  orbit : bool
    Plot orbital data (obliquity, eccentricity, COPP)
    (default = False)
  show : bool
    Show plot in Python (default = True)
  
  Output
  ------
  PDF format plot with name 'evol_<dir>.pdf'
  
  """
  if not isinstance(dir,(list,tuple)):
    dir = [dir]
  
  nfiles = len(dir)

  if nfiles > 1 and orbit == True:
    raise Exception("Error: cannot plot multiple files when orbit = True")
  
  titlestr = []
  
  if orbit == True:
    fig = plt.figure(figsize=(16,13))
  else:
    fig = plt.figure(figsize=(10*nfiles,15))

  for ii in np.arange(nfiles):
    out = vplot.GetOutput(dir[ii])
  
    for p in range(len(out.bodies)):
      if out.bodies[p].name == plname:
        body = out.bodies[p]
        ctmp = 1
      else:
        if p == len(out.bodies)-1 and ctmp == 0:
          raise Exception("Planet %s not found in folder %s"%(plname,folders[j]))
  
    try:
      ecc = body.Eccentricity
    except:
      ecc = np.zeros_like(body.Time)+getattr(out.log.initial,plname).Eccentricity
    
    try:
      inc = body.Inc
    except:
      inc = np.zeros_like(body.Time)
    
    try:
      obl = body.Obliquity
    except:
      obltmp = getattr(out.log.initial,plname).Obliquity
      if obltmp.unit == 'rad':
        obltmp *= 180/np.pi
      obl = np.zeros_like(body.Time)+obltmp

    f = open(dir[ii]+'/'+plname+'.in','r')
    lines = f.readlines()
    f.close()
    for i in range(len(lines)):
      if lines[i].split() != []:
        if lines[i].split()[0] == 'dRotPeriod':
          P = -1*np.float(lines[i].split()[1])  

    try:
      longp = (body.ArgP + body.LongA + body.PrecA)*np.pi/180.0
    except:
      longp = body.PrecA*np.pi/180.0
    
    esinv = ecc*np.sin(longp)*np.sin(obl*np.pi/180.)

  
    titlestr.append(r'$e_0 = %f, i_0 = %f, \psi_0 = %f, P_{rot} = %f$'%(ecc[0],inc[0],obl[0],P)) 
    fig.subplots_adjust(wspace=0.3)

    lats = np.unique(body.Latitude)
    nlats = len(lats)
    ntimes = len(body.Time)

    # plot temperature
    maxscl = np.max(body.TempLandLat)
    minscl = np.min(body.TempLandLat)
    norm = mplcol.Normalize(vmin=minscl,vmax=maxscl)
  
  
    temp = np.reshape(body.TempLat,(ntimes,nlats))
    if orbit == True:
      ax1 = plt.subplot(4,2,1)
    else:
      ax1 = plt.subplot(5,nfiles,ii+1)
    pos = ax1.figbox.get_points()
    c = plt.contourf(body.Time,lats,temp.T,cmap='plasma',norm = norm)
    plt.ylabel('Latitude')
    plt.title('Surface Temp [$^{\circ}$C]')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    tempmin = np.reshape(body.TempLandLat,(ntimes,nlats))
    if orbit == True:
      ax3 = plt.subplot(4,2,3)
    else:
      ax3 = plt.subplot(5,nfiles,ii+2*nfiles+1)
    pos = ax3.figbox.get_points()
    c = plt.contourf(body.Time,lats,tempmin.T,cmap='plasma',norm = norm)
    plt.ylabel('Latitude')
    plt.title('Land annual surface Temp [$^{\circ}$C]')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
    
    tempmin = np.reshape(body.TempMaxLand,(ntimes,nlats))
    if orbit == True:
      ax3 = plt.subplot(4,2,5)
    else:
      ax3 = plt.subplot(5,nfiles,ii+3*nfiles+1)
    pos = ax3.figbox.get_points()
    c = plt.contourf(body.Time,lats,tempmin.T,cmap='plasma')
    plt.ylabel('Latitude')
    plt.title('Land maximum surface Temp [$^{\circ}$C]')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    tempmax = np.reshape(body.TempWaterLat,(ntimes,nlats))
    if orbit == True:
      ax5 = plt.subplot(4,2,7)
    else:
      ax5 = plt.subplot(5,nfiles,ii+4*nfiles+1)
    pos = ax5.figbox.get_points()
    c = plt.contourf(body.Time,lats,tempmax.T,cmap='plasma',norm = norm)
    plt.ylabel('Latitude')
    plt.title('Water annual surface Temp [$^{\circ}$C]')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    tempmax = np.reshape(body.TempMaxWater,(ntimes,nlats))
    if orbit == True:
      ax5 = plt.subplot(4,2,2)
    else:
      ax5 = plt.subplot(5,nfiles,ii+nfiles+1)
    pos = ax5.figbox.get_points()
    c = plt.contourf(body.Time,lats,tempmax.T,cmap='plasma')
    plt.ylabel('Latitude')
    plt.title('Water maximum surface Temp [$^{\circ}$C]')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    if orbit == True:
      #obliquity
      plt.subplot(4,2,4)
      plt.plot(body.Time,obl,linestyle = 'solid',marker='None',color='darkblue',linewidth =2)
      plt.ylabel('Obliquity')
      if xrange:
        plt.xlim(xrange)

      #eccentricity
      plt.subplot(4,2,6)
      plt.plot(body.Time,ecc,linestyle = 'solid',marker='None',color='darkorchid',linewidth =2)
      plt.ylabel('Eccentricity')
      if xrange:
        plt.xlim(xrange)

      #e sin(obl) sin varpi
      plt.subplot(4,2,8)
      plt.plot(body.Time,esinv,linestyle = 'solid',marker='None',color='salmon',linewidth=2)
      plt.ylabel('COPP')
      plt.xlabel('Time [years]')
      if xrange:
        plt.xlim(xrange)

    if dir[ii] == '.':
      dir[ii] = 'cwd'

  fig.suptitle('\n'.join(titlestr),fontsize=20) 

  if xrange:
    sfile = 'tempLW_'+'_'.join(dir)+'_%d_%d.pdf'%(xrange[0],xrange[1])
  else:
    sfile = 'tempLW_'+'_'.join(dir)+'.pdf'
  plt.savefig(sfile)
  if show:
    plt.show()
  else:
    plt.close()

def temp_spec(plname,dir='.',xrange=False,orbit=True,show=True):
  """
  Creates plots of average temperature, min temp, and max temp
  over the length of the simulation
  
  Parameters
  ----------
  plname : string
    The name of the planet with .Climate data
    
  Keyword Arguments
  -----------------
  dir : string
    Directory of vplanet simulation (default = '.')
  xrange : float tuple, list, or numpy array
    Range of x-values (time) to restrict plot
    (default = False (no restriction))
  orbit : bool
    Plot orbital data (obliquity, eccentricity, COPP)
    (default = False)
  show : bool
    Show plot in Python (default = True)
  
  Output
  ------
  PDF format plot with name 'evol_<dir>.pdf'
  
  """
  if not isinstance(dir,(list,tuple)):
    dir = [dir]
  
  nfiles = len(dir)

  if nfiles > 1 and orbit == True:
    raise Exception("Error: cannot plot multiple files when orbit = True")
  
  titlestr = []
  
  if orbit == True:
    fig = plt.figure(figsize=(14,7))
  else:
    fig = plt.figure(figsize=(10*nfiles,10))

  hdffile = 0
  if os.path.exists('climatedata.hdf5'):
    hdff = h5py.File('climatedata.hdf5')
    hdffile = 1

  for ii in np.arange(nfiles):
    if hdffile:
      time = hdff[dir[ii]]['Time'][...]
      try:
        ecc = hdff[dir[ii]]['Eccentricity'][...]
      except:
        ecc = np.zeros_like(time)+hdff[dir[ii]]['EccInit'][...]
        
      try:
        inc = (hdff[dir[ii]]['Inc'][...])
      except:
        inc = np.zeros_like(time)
        
      try:  
        obl = (hdff[dir[ii]]['Obliquity'][...])
      except:
        obl = np.zeros_like(time)+hdff[dir[ii]]['OblInit'][...]
        
      preca = hdff[dir[ii]]['PrecA'][...]
      try:
        argp = (hdff[dir[ii]]['ArgP'][...])
        longa = hdff[dir[ii]]['LongA'][...]
        longp = (argp+longa+preca)*np.pi/180.0
      except:
        longp = preca*np.pi/180.0
        
      lats = np.unique(hdff[dir[ii]]['Latitude'][...])
      TempLat = hdff[dir[ii]]['TempLat'][...]
      TempMinLat = hdff[dir[ii]]['TempMinLat'][...]
      TempMaxLat = hdff[dir[ii]]['TempMaxLat'][...]
      
    else:
      out = vplot.GetOutput(dir[ii])
  
      for p in range(len(out.bodies)):
        if out.bodies[p].name == plname:
          body = out.bodies[p]
          ctmp = 1
        else:
          if p == len(out.bodies)-1 and ctmp == 0:
            raise Exception("Planet %s not found in folder %s"%(plname,folders[j]))
  
      try:
        ecc = body.Eccentricity
      except:
        ecc = np.zeros_like(body.Time)+getattr(out.log.initial,plname).Eccentricity
    
      try:
        inc = body.Inc
      except:
        inc = np.zeros_like(body.Time)
    
      try:
        obl = body.Obliquity
      except:
        obltmp = getattr(out.log.initial,plname).Obliquity
        if obltmp.unit == 'rad':
          obltmp *= 180/np.pi
        obl = np.zeros_like(body.Time)+obltmp
        
      try:
        longp = (body.ArgP + body.LongA + body.PrecA)*np.pi/180.0
      except:
        longp = body.PrecA*np.pi/180.0
        
      lats = np.unique(body.Latitude)
      time = body.Time
      TempLat = body.TempLat
      TempMinLat = body.TempMinLat
      TempMaxLat = body.TempMaxLat
    
    f = open(dir[ii]+'/'+plname+'.in','r')
    lines = f.readlines()
    f.close()
    for i in range(len(lines)):
      if lines[i].split() != []:
        if lines[i].split()[0] == 'dRotPeriod':
          P = -1*np.float(lines[i].split()[1])  


    
    esinv = ecc*np.sin(longp)*np.sin(obl*np.pi/180.)

  
    titlestr.append(r'$e_0 = %f, i_0 = %f, \psi_0 = %f, P_{rot} = %f$'%(ecc[0],inc[0],obl[0],P)) 
    fig.subplots_adjust(wspace=0.3,hspace=0.15)

    #lats = np.unique(body.Latitude)
    nlats = len(lats)
    ntimes = len(time)

    # plot temperature
    maxscl = np.max(TempMaxLat)
    minscl = np.min(TempMinLat)
    norm = mplcol.Normalize(vmin=minscl,vmax=maxscl)
  
  
    temp = np.reshape(TempLat,(ntimes,nlats))
    # if orbit == True:
#       ax1 = plt.subplot(3,2,1)
#     else:
#       ax1 = plt.subplot(3,nfiles,ii+1)
#     pos = ax1.figbox.get_points()
#     c = plt.contourf(time,lats,temp.T,cmap='plasma',norm = norm)
#     plt.ylabel('Latitude')
#     plt.title('Average Surface Temp [$^{\circ}$C]')
#     plt.ylim(-90,90)
#     plt.yticks([-60,-30,0,30,60])
#     if xrange:
#       plt.xlim(xrange)
#     plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
#   
#     tempmin = np.reshape(TempMinLat,(ntimes,nlats))
#     if orbit == True:
#       ax3 = plt.subplot(3,2,3)
#     else:
#       ax3 = plt.subplot(3,nfiles,ii+nfiles+1)
#     pos = ax3.figbox.get_points()
#     c = plt.contourf(time,lats,tempmin.T,cmap='plasma',norm = norm)
#     plt.ylabel('Latitude')
#     plt.title('Minimum Surface Temp [$^{\circ}$C]')
#     plt.ylim(-90,90)
#     plt.yticks([-60,-30,0,30,60])
#     if xrange:
#       plt.xlim(xrange)
#     plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    tempmax = np.reshape(TempMaxLat,(ntimes,nlats))
    if orbit == True:
      ax5 = plt.subplot(2,2,1)
    else:
      ax5 = plt.subplot(3,nfiles,ii+2*nfiles+1)
    pos = ax5.figbox.get_points()
    c = plt.contourf(time,lats,tempmax.T,cmap='plasma',norm = norm)
    plt.ylabel('Latitude ($^{\circ}$)')
    plt.title('Maximum Surface Temp [$^{\circ}$C]')
    plt.ylim(np.min(lats),np.max(lats))
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    # tempmax = np.reshape(TempMaxLat,(ntimes,nlats))
#     if orbit == True:
#       ax5 = plt.subplot(2,2,3)
#     else:
#       ax5 = plt.subplot(3,nfiles,ii+2*nfiles+1)
#     pos = ax5.figbox.get_points()
#     c = plt.contourf(time,lats,tempmax.T,cmap='plasma',norm = norm)
#     plt.ylabel('Latitude')
#     plt.title('Maximum Surface Temp [$^{\circ}$C]')
#     plt.ylim(-90,90)
#     plt.yticks([-60,-30,0,30,60])
#     if xrange:
#       plt.xlim(2.5e5,5e5)
#     plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    if orbit == True:
      #obliquity
      plt.subplot(2,2,3)
      plt.plot(time,obl,linestyle = 'solid',marker='None',color='darkblue',linewidth =2)
      plt.ylabel('Obliquity ($^{\circ}$)')
      if xrange:
        plt.xlim(xrange)
      plt.xlabel('Time (years)')

      #eccentricity
      plt.subplot(2,2,2)
      plt.plot(time,ecc,linestyle = 'solid',marker='None',color='darkorchid',linewidth =2)
      plt.ylabel('Eccentricity')
      if xrange:
        plt.xlim(xrange)
      

      #e sin(obl) sin varpi
      plt.subplot(2,2,4)
      plt.plot(time,esinv,linestyle = 'solid',marker='None',color='salmon',linewidth=2)
      plt.ylabel('COPP')
      plt.xlabel('Time (years)')
      if xrange:
        plt.xlim(xrange)

    if dir[ii] == '.':
      dir[ii] = 'cwd'

#   fig.suptitle('\n'.join(titlestr),fontsize=20) 

  if xrange:
    sfile = 'tspec_'+'_'.join(dir)+'_%d_%d.pdf'%(xrange[0],xrange[1])
  else:
    sfile = 'tspec_'+'_'.join(dir)+'.pdf'
  plt.savefig(sfile)
  if show:
    plt.show()
  else:
    plt.close()

def tempminmax(plname,dir='.',xrange=False,orbit=False,show=True):
  """
  Creates plots of average temperature, min temp, and max temp
  over the length of the simulation
  
  Parameters
  ----------
  plname : string
    The name of the planet with .Climate data
    
  Keyword Arguments
  -----------------
  dir : string
    Directory of vplanet simulation (default = '.')
  xrange : float tuple, list, or numpy array
    Range of x-values (time) to restrict plot
    (default = False (no restriction))
  orbit : bool
    Plot orbital data (obliquity, eccentricity, COPP)
    (default = False)
  show : bool
    Show plot in Python (default = True)
  
  Output
  ------
  PDF format plot with name 'evol_<dir>.pdf'
  
  """
  if not isinstance(dir,(list,tuple)):
    dir = [dir]
  
  nfiles = len(dir)

  if nfiles > 1 and orbit == True:
    raise Exception("Error: cannot plot multiple files when orbit = True")
  
  titlestr = []
  
  if orbit == True:
    fig = plt.figure(figsize=(16,10))
  else:
    fig = plt.figure(figsize=(10*nfiles,10))

  hdffile = 0
  if os.path.exists('climatedata.hdf5'):
    hdff = h5py.File('climatedata.hdf5')
    hdffile = 1

  for ii in np.arange(nfiles):
    if hdffile:
      time = hdff[dir[ii]]['Time'][...]
      try:
        ecc = hdff[dir[ii]]['Eccentricity'][...]
      except:
        ecc = np.zeros_like(time)+hdff[dir[ii]]['EccInit'][...]
        
      try:
        inc = (hdff[dir[ii]]['Inc'][...])
      except:
        inc = np.zeros_like(time)
        
      try:  
        obl = (hdff[dir[ii]]['Obliquity'][...])
      except:
        obl = np.zeros_like(time)+hdff[dir[ii]]['OblInit'][...]
        
      preca = hdff[dir[ii]]['PrecA'][...]
      try:
        argp = (hdff[dir[ii]]['ArgP'][...])
        longa = hdff[dir[ii]]['LongA'][...]
        longp = (argp+longa+preca)*np.pi/180.0
      except:
        longp = preca*np.pi/180.0
        
      lats = np.unique(hdff[dir[ii]]['Latitude'][...])
      TempLat = hdff[dir[ii]]['TempLat'][...]
      TempMinLat = hdff[dir[ii]]['TempMinLat'][...]
      TempMaxLat = hdff[dir[ii]]['TempMaxLat'][...]
      
    else:
      out = vplot.GetOutput(dir[ii])
  
      for p in range(len(out.bodies)):
        if out.bodies[p].name == plname:
          body = out.bodies[p]
          ctmp = 1
        else:
          if p == len(out.bodies)-1 and ctmp == 0:
            raise Exception("Planet %s not found in folder %s"%(plname,folders[j]))
  
      try:
        ecc = body.Eccentricity
      except:
        ecc = np.zeros_like(body.Time)+getattr(out.log.initial,plname).Eccentricity
    
      try:
        inc = body.Inc
      except:
        inc = np.zeros_like(body.Time)
    
      try:
        obl = body.Obliquity
      except:
        obltmp = getattr(out.log.initial,plname).Obliquity
        if obltmp.unit == 'rad':
          obltmp *= 180/np.pi
        obl = np.zeros_like(body.Time)+obltmp
        
      try:
        longp = (body.ArgP + body.LongA + body.PrecA)*np.pi/180.0
      except:
        longp = body.PrecA*np.pi/180.0
        
      lats = np.unique(body.Latitude)
      time = body.Time
      TempLat = body.TempLat
      TempMinLat = body.TempMinLat
      TempMaxLat = body.TempMaxLat
    
    f = open(dir[ii]+'/'+plname+'.in','r')
    lines = f.readlines()
    f.close()
    for i in range(len(lines)):
      if lines[i].split() != []:
        if lines[i].split()[0] == 'dRotPeriod':
          P = -1*np.float(lines[i].split()[1])  


    
    esinv = ecc*np.sin(longp)*np.sin(obl*np.pi/180.)

  
    titlestr.append(r'$e_0 = %f, i_0 = %f, \psi_0 = %f, P_{rot} = %f$'%(ecc[0],inc[0],obl[0],P)) 
    fig.subplots_adjust(wspace=0.3)

    #lats = np.unique(body.Latitude)
    nlats = len(lats)
    ntimes = len(time)

    # plot temperature
    maxscl = np.max(TempMaxLat)
    minscl = np.min(TempMinLat)
    norm = mplcol.Normalize(vmin=minscl,vmax=maxscl)
  
  
    temp = np.reshape(TempLat,(ntimes,nlats))
    if orbit == True:
      ax1 = plt.subplot(3,2,1)
    else:
      ax1 = plt.subplot(3,nfiles,ii+1)
    pos = ax1.figbox.get_points()
    c = plt.contourf(time,lats,temp.T,cmap='plasma',norm = norm)
    plt.ylabel('Latitude')
    plt.title('Average Surface Temp [$^{\circ}$C]')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    tempmin = np.reshape(TempMinLat,(ntimes,nlats))
    if orbit == True:
      ax3 = plt.subplot(3,2,3)
    else:
      ax3 = plt.subplot(3,nfiles,ii+nfiles+1)
    pos = ax3.figbox.get_points()
    c = plt.contourf(time,lats,tempmin.T,cmap='plasma',norm = norm)
    plt.ylabel('Latitude')
    plt.title('Minimum Surface Temp [$^{\circ}$C]')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    tempmax = np.reshape(TempMaxLat,(ntimes,nlats))
    if orbit == True:
      ax5 = plt.subplot(3,2,5)
    else:
      ax5 = plt.subplot(3,nfiles,ii+2*nfiles+1)
    pos = ax5.figbox.get_points()
    c = plt.contourf(time,lats,tempmax.T,cmap='plasma',norm = norm)
    plt.ylabel('Latitude')
    plt.title('Maximum Surface Temp [$^{\circ}$C]')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    if orbit == True:
      #obliquity
      plt.subplot(3,2,2)
      plt.plot(time,obl,linestyle = 'solid',marker='None',color='darkblue',linewidth =2)
      plt.ylabel('Obliquity')
      if xrange:
        plt.xlim(xrange)

      #eccentricity
      plt.subplot(3,2,4)
      plt.plot(time,ecc,linestyle = 'solid',marker='None',color='darkorchid',linewidth =2)
      plt.ylabel('Eccentricity')
      if xrange:
        plt.xlim(xrange)

      #e sin(obl) sin varpi
      plt.subplot(3,2,6)
      plt.plot(time,esinv,linestyle = 'solid',marker='None',color='salmon',linewidth=2)
      plt.ylabel('COPP')
      plt.xlabel('Time [years]')
      if xrange:
        plt.xlim(xrange)

    if dir[ii] == '.':
      dir[ii] = 'cwd'

#   fig.suptitle('\n'.join(titlestr),fontsize=20) 

  if xrange:
    sfile = 'temp_'+'_'.join(dir)+'_%d_%d.pdf'%(xrange[0],xrange[1])
  else:
    sfile = 'temp_'+'_'.join(dir)+'.pdf'
  plt.savefig(sfile)
  if show:
    plt.show()
  else:
    plt.close()
    
def flux_evol(plname,dir='.',xrange=False,orbit=False,show=True):
  """
  Creates plots of average temperature, min temp, and max temp
  over the length of the simulation
  
  Parameters
  ----------
  plname : string
    The name of the planet with .Climate data
    
  Keyword Arguments
  -----------------
  dir : string
    Directory of vplanet simulation (default = '.')
  xrange : float tuple, list, or numpy array
    Range of x-values (time) to restrict plot
    (default = False (no restriction))
  orbit : bool
    Plot orbital data (obliquity, eccentricity, COPP)
    (default = False)
  show : bool
    Show plot in Python (default = True)
  
  Output
  ------
  PDF format plot with name 'evol_<dir>.pdf'
  
  """
  if not isinstance(dir,(list,tuple)):
    dir = [dir]
  
  nfiles = len(dir)

  if nfiles > 1 and orbit == True:
    raise Exception("Error: cannot plot multiple files when orbit = True")
  
  titlestr = []
  
  if orbit == True:
    fig = plt.figure(figsize=(16,12))
  else:
    fig = plt.figure(figsize=(10*nfiles,12))

  for ii in np.arange(nfiles):
    out = vplot.GetOutput(dir[ii])
  
    for p in range(len(out.bodies)):
      if out.bodies[p].name == plname:
        body = out.bodies[p]
        ctmp = 1
      else:
        if p == len(out.bodies)-1 and ctmp == 0:
          raise Exception("Planet %s not found in folder %s"%(plname,folders[j]))
  
    try:
      ecc = body.Eccentricity
    except:
      ecc = np.zeros_like(body.Time)+getattr(out.log.initial,plname).Eccentricity
    
    try:
      inc = body.Inc
    except:
      inc = np.zeros_like(body.Time)
    
    try:
      obl = body.Obliquity
    except:
      obltmp = getattr(out.log.initial,plname).Obliquity
      if obltmp.unit == 'rad':
        obltmp *= 180/np.pi
      obl = np.zeros_like(body.Time)+obltmp

    f = open(dir[ii]+'/'+plname+'.in','r')
    lines = f.readlines()
    f.close()
    for i in range(len(lines)):
      if lines[i].split() != []:
        if lines[i].split()[0] == 'dRotPeriod':
          P = -1*np.float(lines[i].split()[1])  

    try:
      longp = (body.ArgP + body.LongA + body.PrecA)*np.pi/180.0
    except:
      longp = body.PrecA*np.pi/180.0
    
    esinv = ecc*np.sin(longp)*np.sin(obl*np.pi/180.)

  
    titlestr.append(r'$e_0 = %f, i_0 = %f, \psi_0 = %f, P_{rot} = %f$'%(ecc[0],inc[0],obl[0],P)) 
    fig.subplots_adjust(wspace=0.3)

    lats = np.unique(body.Latitude)
    nlats = len(lats)
    ntimes = len(body.Time)

    # plot temperature
#     maxscl = np.max(body.TempMaxLat)
#     minscl = np.min(body.TempMinLat)
#     norm = mplcol.Normalize(vmin=minscl,vmax=maxscl)
  
  
    temp = np.reshape(body.TempLat,(ntimes,nlats))
    if orbit == True:
      ax1 = plt.subplot(4,2,1)
    else:
      ax1 = plt.subplot(4,nfiles,ii+1)
    pos = ax1.figbox.get_points()
    c = plt.contourf(body.Time,lats,temp.T,cmap='plasma')
    plt.ylabel('Latitude')
    plt.title('Surface Temp [$^{\circ}$C]')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    insol = np.reshape(body.AnnInsol,(ntimes,nlats))
    if orbit == True:
      ax3 = plt.subplot(4,2,3)
    else:
      ax3 = plt.subplot(4,nfiles,ii+nfiles+1)
    pos = ax3.figbox.get_points()
    c = plt.contourf(body.Time,lats,insol.T,cmap='plasma')
    plt.ylabel('Latitude')
    plt.title('Insolation [W/m$^2$]')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    olr = np.reshape(body.FluxOut,(ntimes,nlats))
    if orbit == True:
      ax5 = plt.subplot(4,2,5)
    else:
      ax5 = plt.subplot(4,nfiles,ii+2*nfiles+1)
    pos = ax5.figbox.get_points()
    c = plt.contourf(body.Time,lats,olr.T,cmap='plasma')
    plt.ylabel('Latitude')
    plt.title('Outgoing flux [W/m$^2$]')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
    
    plb = np.reshape(body.PlanckBAvg,(ntimes,nlats))
    if orbit == True:
      ax7 = plt.subplot(4,2,7)
    else:
      ax7 = plt.subplot(4,nfiles,ii+3*nfiles+1)
    pos = ax7.figbox.get_points()
    c = plt.contourf(body.Time,lats,plb.T,cmap='plasma')
    plt.ylabel('Latitude')
    plt.title('Planck B coeff [W/m$^2$/$^{\circ}$C]')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    if orbit == True:
      #obliquity
      plt.subplot(4,2,2)
      plt.plot(body.Time,obl,linestyle = 'solid',marker='None',color='darkblue',linewidth =2)
      plt.ylabel('Obliquity')
      if xrange:
        plt.xlim(xrange)

      #eccentricity
      plt.subplot(4,2,4)
      plt.plot(body.Time,ecc,linestyle = 'solid',marker='None',color='darkorchid',linewidth =2)
      plt.ylabel('Eccentricity')
      if xrange:
        plt.xlim(xrange)

      #e sin(obl) sin varpi
      plt.subplot(4,2,6)
      plt.plot(body.Time,esinv,linestyle = 'solid',marker='None',color='salmon',linewidth=2)
      plt.ylabel('COPP')
      plt.xlabel('Time [years]')
      if xrange:
        plt.xlim(xrange)

    if dir[ii] == '.':
      dir[ii] = 'cwd'

  fig.suptitle('\n'.join(titlestr),fontsize=20) 

  if xrange:
    sfile = 'temp_'+'_'.join(dir)+'_%d_%d.pdf'%(xrange[0],xrange[1])
  else:
    sfile = 'temp_'+'_'.join(dir)+'.pdf'
  plt.savefig(sfile)
  if show:
    plt.show()
  else:
    plt.close()

def ice_lats(time, lat, iceh):
  icelatsouth = np.zeros(len(time))
  icelatnorth = np.zeros(len(time))
  
  for i in range(len(time)):
    ice = lat[i][iceh[i]>0]
    if len(ice[ice<0]) != 0:
      icelatsouth[i] = np.max(ice[ice<0])
    else:
      icelatsouth[i] = -90
    if len(ice[ice>0]) != 0:   
      icelatnorth[i] = np.min(ice[ice>0])
    else:
      icelatnorth[i] = 90
  
  return icelatsouth, icelatnorth

def ice_fft(plname,dir='.',log = False):
  out = vplot.GetOutput(dir)

  ctmp = 0
  for p in range(len(out.bodies)):
      if out.bodies[p].name == plname:
        body = out.bodies[p]
        ctmp = 1
      else:
        if p == len(out.bodies)-1 and ctmp == 0:
          raise Exception("Planet %s not found in folder %s"%(plname,dir[ii]))

  icelatsouth, icelatnorth = ice_lats(body.Time, body.Latitude, body.IceHeight)
  n65 = np.where(np.abs(body.Latitude[0]-65)==np.min(np.abs(body.Latitude[0]-65)))[0]
  s65 = np.where(np.abs(body.Latitude[0]+65)==np.min(np.abs(body.Latitude[0]+65)))[0]
  norths = np.where(body.Latitude[0]>=0)[0]
  souths = np.where(body.Latitude[0]<0)[0]

  icehsouth = body.IceHeight[:,s65[0]]
  icehnorth = body.IceHeight[:,n65[0]]

  icemsouth = np.sum(body.IceMass[:,souths],1)
  icemnorth = np.sum(body.IceMass[:,norths],1)

  #----ice data periodogramaphones----------------  
  #datasouth = icelatsouth[50:]-np.mean(icelatsouth[50:])
  #datanorth = icelatnorth[50:]-np.mean(icelatnorth[50:])
#   datasouth = icehsouth[50:]-np.mean(icehsouth[50:])# 
#   datanorth = icehnorth[50:]-np.mean(icehnorth[50:])
  datasouth = icemsouth[50:]-np.mean(icemsouth[50:])
  datanorth = icemnorth[50:]-np.mean(icemnorth[50:])
  datatotal = body.TotIceMass[50:] - np.mean(body.TotIceMass[50:])
  
  freqs, powsouth = sig.periodogram(datasouth,fs=0.001,window='bartlett')
  freqs, pownorth = sig.periodogram(datanorth,fs=0.001,window='bartlett')
  freqs, powtotal = sig.periodogram(datatotal,fs=0.001,window='bartlett')

  powsouth *= 1./np.max(powsouth)
  pownorth *= 1./np.max(pownorth)
  powtotal *= 1./np.max(powtotal)
  
  #----milankovitch data periodododoodoo-----------
  try:
    dataobliq = body.Obliquity - np.mean(body.Obliquity)
    dataeccen = body.Eccentricity - np.mean(body.Eccentricity)
    COPP = body.Eccentricity*np.sin(body.Obliquity*np.pi/180.0) *np.sin((body.ArgP+body.LongA+body.PrecA)*np.pi/180.0)  
  #   COPP = body.Eccentricity *np.sin((body.ArgP+body.LongA+body.PrecA)*np.pi/180.0)  
    datacopp = COPP - np.mean(COPP)
    orbitdata = 1
  except:
    orbitdata = 0
    e = getattr(out.log.initial,plname).Eccentricity
    obl = getattr(out.log.initial,plname).Obliquity
    COPP = e*np.sin(obl) *np.sin((body.PrecA)*np.pi/180.0)  
    datacopp = COPP - np.mean(COPP)
    
  period = 1./freqs
  
  if orbitdata:
    freqs0, powobliq = sig.periodogram(dataobliq,fs=0.001,window='bartlett')
    freqs0, poweccen = sig.periodogram(dataeccen,fs=0.001,window='bartlett')

    powobliq *= 1./np.max(powobliq)
    poweccen *= 1./np.max(poweccen)
    #-------------------------------------------
  
    period0 = 1./freqs0
    eccxs = [1,1]*period0[poweccen ==np.max(poweccen)]
    eccxs1 = [1,1]*period0[poweccen ==np.max(poweccen[period0<30000])]
    oblxs = [1,1]*period0[powobliq == np.max(powobliq)]
    
#   import pdb;pdb.set_trace()

  ys = [-0.2, 100]
  freqs0, powcopp = sig.periodogram(datacopp,fs=0.001,window='bartlett')
  period0 = 1./freqs0
  powcopp *= 1./np.max(powcopp)
  coppxs0 = [1,1]*period0[powcopp == np.max(powcopp)]
  coppxs1 = [1,1]*period0[powcopp == np.max(powcopp[powcopp<.99])]


  fig = plt.figure(figsize=(10,8))
  fig.subplots_adjust(hspace=0.0)
  ax1= plt.subplot(3,1,1)
  plt.semilogx(period,powsouth,linestyle='--',color=vplorg,marker='None',lw=2,label=r'Ice height (65$^{\circ}$ S)',zorder = 3)
  plt.semilogx(period,pownorth,linestyle='-',color=vpldbl,marker='None',lw=1.5,label='Ice height (65$^{\circ}$ N)')
  if orbitdata:
    plt.plot(eccxs,ys,linestyle='--',color=vplred,marker='None',lw=1,zorder=1)
    plt.plot(eccxs1,ys,linestyle='--',color=vplred,marker='None',lw=1,zorder=1)
    plt.plot(oblxs,ys,linestyle='--',color=vplpur,marker='None',lw=1,zorder=1)
   
  plt.plot(coppxs0,ys,linestyle='--',color=vpllbl,marker='None',lw=1,zorder=1)
  plt.plot(coppxs1,ys,linestyle='--',color=vpllbl,marker='None',lw=1,zorder=1)
#   plt.title('Power spectrum (normalized to peak)')
  plt.legend(loc='upper right')
  plt.xticks(visible = False)
  plt.xlim(10000,3e5)  
  if log:
    plt.ylim(1e-10,10)
    ax1.set_yscale('log')
    plt.yticks([1e-9,1e-8,1e-6,1e-7,1e-5,1e-4,1e-3,1e-2,1e-1,1])
  else:
    plt.ylim(-0.2,1.2)
    plt.yticks([0,.2,.4,.6,.8,1])

  ax2=plt.subplot(3,1,2)
  plt.semilogx(period,powtotal,linestyle='-',color='k',marker='None',lw=2,label='Global ice mass')
#   plt.xlabel('Period [years]')
  if orbitdata:
    plt.plot(eccxs,ys,linestyle='--',color=vplred,marker='None',lw=1,zorder=1)
    plt.plot(eccxs1,ys,linestyle='--',color=vplred,marker='None',lw=1,zorder=1)
    plt.plot(oblxs,ys,linestyle='--',color=vplpur,marker='None',lw=1,zorder=1)
  
  plt.plot(coppxs0,ys,linestyle='--',color=vpllbl,marker='None',lw=1,zorder=1)
  plt.plot(coppxs1,ys,linestyle='--',color=vpllbl,marker='None',lw=1,zorder=1)
  plt.legend(loc='upper right')
  plt.xticks(visible = False)
  plt.xlim(10000,3e5) 
  if log:
    plt.ylim(1e-10,10)
    ax2.set_yscale('log')
    plt.yticks([1e-9,1e-8,1e-6,1e-7,1e-5,1e-4,1e-3,1e-2,1e-1,1])
  else:
    plt.ylim(-0.2,1.2)
    plt.yticks([0,.2,.4,.6,.8,1])
  
  ax3=plt.subplot(3,1,3)
  if orbitdata:
    plt.semilogx(period0,powobliq,linestyle='-',color=vplpur,marker='None',lw=2,label='Obliquity')
    plt.semilogx(period0,poweccen,linestyle='-',color=vplred,marker='None',lw=2,label='Eccentricity')
  
  plt.semilogx(period0,powcopp,linestyle='-',color=vpllbl,marker='None',lw=2,label='COPP')
  plt.legend(loc='upper right')
  plt.xlabel('Period [years]')
  plt.xlim(10000,3e5)
  if log:
    plt.ylim(1e-11,10)
    ax3.set_yscale('log')
    plt.yticks([1e-10,1e-9,1e-8,1e-6,1e-7,1e-5,1e-4,1e-3,1e-2,1e-1,1])
  else:
    plt.ylim(-0.2,1.2)
    plt.yticks([0,.2,.4,.6,.8,1])

  if log:
    plt.savefig('periodogram_'+dir+'_log.pdf')
  else: 
    plt.savefig('periodogram_'+dir+'.pdf')
  plt.close()
 
def init_icelines(line):
  line.set_data([],[])
  line1.set_data([],[])
  line2.set_data([],[])
  return line, line1, line2
  
def ice_lines(i,body,line,line1,line2):
#   print(i)
  line.set_data(body.Latitude[i],(body.IceHeight+body.BedrockH)[i])
  line1.set_data(body.Latitude[i],body.BedrockH[i])
  line2.set_data(body.Latitude[i],body.IceFlow[i]*3.15576e7*1000)
  return line, line1, line2

def ice_video(plname,dir='.'):
  if not isinstance(dir,(list,tuple)):
    dir = [dir]
  
  nfiles = len(dir)

  if nfiles > 1 and orbit == True:
    raise Exception("Error: cannot plot multiple files when orbit = True")
  
  for ii in np.arange(nfiles):
    out = vplot.GetOutput(dir[ii])
    
    #pdb.set_trace()
  
    ctmp = 0
    for p in range(len(out.bodies)):
      if out.bodies[p].name == plname:
        body = out.bodies[p]
        ctmp = 1
      else:
        if p == len(out.bodies)-1 and ctmp == 0:
          raise Exception("Planet %s not found in folder %s"%(plname,dir[ii]))
  
  
  os.system('mkdir '+dir[0]+'/ice_frames')
  
  max = np.max(body.IceHeight+body.BedrockH)
  min = np.min(body.BedrockH)
  
  fig = plt.figure(figsize=(10,8))
  ax = fig.add_subplot(111)
  plt.xlim(-85,85)
  plt.ylim(-1000,2000)
  line, = ax.plot([],[],'-',color='k')
  line1, = ax.plot([],[],'-',color='r')
  line2, = ax.plot([],[],'-',color='b')
#   pdb.set_trace()
  anim = animation.FuncAnimation(fig, ice_lines, np.arange(0,len(body.Time)),interval=100,blit=True,fargs=(body,line,line1,line2))
  
  plt.show()
  
  # for i in np.arange(len(body.Time)):
#     
#     fig = plt.figure(figsize=(10,5))
#     fig.subplots_adjust(wspace=0.3)
#     plt.subplot(1,2,1)
#     plt.plot(body.Latitude[i],body.IceHeight[i]+body.BedrockH[i],color=vpllbl)
#     plt.plot(body.Latitude[i],body.BedrockH[i],color=vplred)
#     
#     plt.xlabel('Height [m]')
#     plt.ylabel('Latitude (deg)')
#     plt.xlim(-90,-40)
#     plt.ylim(min,max)
#     plt.text(-70,0.95*max,'Time=%#.2f year'%body.Time[i])
#     
#     plt.subplot(1,2,2)
#     plt.plot(body.Latitude[i],body.IceHeight[i]+body.BedrockH[i],color=vpllbl)
#     plt.plot(body.Latitude[i],body.BedrockH[i],color=vplred)
#     
#     plt.xlabel('Height [m]')
#     plt.ylabel('Latitude (deg)')
#     plt.xlim(40,90)
#     plt.ylim(min,max)
# #     plt.text(-20,0.95*max,'Time=%#.2f year'%body.Time[i])
#     
#     plt.savefig(dir[0]+'/ice_frames/frame%05d.png'%i)
#     plt.close()



def ice_evol(plname,dir='.',log = False):
  out = vplot.GetOutput(dir)

  ctmp = 0
  for p in range(len(out.bodies)):
      if out.bodies[p].name == plname:
        body = out.bodies[p]
        ctmp = 1
      else:
        if p == len(out.bodies)-1 and ctmp == 0:
          raise Exception("Planet %s not found in folder %s"%(plname,dir[ii]))

  icelatsouth, icelatnorth = ice_lats(body.Time, body.Latitude, body.IceHeight)
  n65 = np.where(np.abs(body.Latitude[0]-65)==np.min(np.abs(body.Latitude[0]-65)))[0]
  s65 = np.where(np.abs(body.Latitude[0]+65)==np.min(np.abs(body.Latitude[0]+65)))[0]
  
  n40 = np.where(np.abs(body.Latitude[0]-40)==np.min(np.abs(body.Latitude[0]-40)))[0]
  s40 = np.where(np.abs(body.Latitude[0]+40)==np.min(np.abs(body.Latitude[0]+40)))[0]
  
  n10 = np.where(np.abs(body.Latitude[0]-10)==np.min(np.abs(body.Latitude[0]-10)))[0]
  s10 = np.where(np.abs(body.Latitude[0]+10)==np.min(np.abs(body.Latitude[0]+10)))[0]
  
  norths = np.where(body.Latitude[0]>=0)[0]
  souths = np.where(body.Latitude[0]<0)[0]

  icehsouth = body.IceHeight[:,s65[0]]
  icehnorth = body.IceHeight[:,n65[0]]
  insol65n = body.AnnInsol[:,n65[0]]
  insol65s = body.AnnInsol[:,s65[0]]
  
  
  insol40n = body.AnnInsol[:,n40[0]]
  insol40s = body.AnnInsol[:,s40[0]]
  
  insol10n = body.AnnInsol[:,n10[0]]
  insol10s = body.AnnInsol[:,s10[0]]
  

  icemsouth = np.sum(body.IceMass[:,souths],1)
  icemnorth = np.sum(body.IceMass[:,norths],1)
  
  fig = plt.figure(figsize=(10,12))
  fig.suptitle('ice mass (northern hemi)')
  ax1 = plt.subplot(3,1,1)
  plt.plot(body.Time,icemnorth,linestyle='-',color=vpldbl)
  ax2 = ax1.twinx()
  plt.plot(body.Time,body.Obliquity,linestyle='--',color=vpllbl)
  
  ax1 = plt.subplot(3,1,2)
  plt.plot(body.Time,icemnorth,linestyle='-',color=vpldbl)
  ax2 = ax1.twinx()
  plt.plot(body.Time,body.Eccentricity,linestyle='--',color=vpllbl)
  
  
  COPP = body.Eccentricity*np.sin(body.Obliquity*np.pi/180.0) *np.sin((body.ArgP+body.LongA+body.PrecA)*np.pi/180.0) 
  
  ax1 = plt.subplot(3,1,3)
  plt.plot(body.Time,icemnorth,linestyle='-',color=vpldbl)
  ax2 = ax1.twinx()
  plt.plot(body.Time,COPP,linestyle='--',color=vpllbl)
  
  
  fig = plt.figure(figsize=(10,12))
  fig.suptitle('ice mass (southern hemi)')
  ax1 = plt.subplot(3,1,1)
  plt.plot(body.Time,icemsouth,linestyle='-',color=vplpur)
  ax2 = ax1.twinx()
  plt.plot(body.Time,body.Obliquity,linestyle='--',color=vpllbl)
  
  ax1 = plt.subplot(3,1,2)
  plt.plot(body.Time,icemsouth,linestyle='-',color=vpldbl)
  ax2 = ax1.twinx()
  plt.plot(body.Time,body.Eccentricity,linestyle='--',color=vplpur)
  
  
  COPP = body.Eccentricity*np.sin(body.Obliquity*np.pi/180.0) *np.sin((body.ArgP+body.LongA+body.PrecA)*np.pi/180.0) 
  
  
  ax1 = plt.subplot(3,1,3)
  plt.plot(body.Time,icemsouth,linestyle='-',color=vplpur)
  ax2 = ax1.twinx()
  plt.plot(body.Time,COPP,linestyle='--',color=vpllbl)
  
  
  fig = plt.figure(figsize=(10,12))
  fig.suptitle('ice mass (total)')
  ax1 = plt.subplot(3,1,1)
  plt.plot(body.Time,body.TotIceMass,linestyle='-',color=vplred)
  ax2 = ax1.twinx()
  plt.plot(body.Time,body.Obliquity,linestyle='--',color=vpllbl)
  
  ax1 = plt.subplot(3,1,2)
  plt.plot(body.Time,body.TotIceMass,linestyle='-',color=vplred)
  ax2 = ax1.twinx()
  plt.plot(body.Time,body.Eccentricity,linestyle='--',color=vplpur)
  
  
  COPP = body.Eccentricity*np.sin(body.Obliquity*np.pi/180.0) *np.sin((body.ArgP+body.LongA+body.PrecA)*np.pi/180.0) 
  
  ax1 = plt.subplot(3,1,3)
  plt.plot(body.Time,body.TotIceMass,linestyle='-',color=vplred)
  ax2 = ax1.twinx()
  plt.plot(body.Time,COPP,linestyle='--',color=vpllbl)
  
  
  fig = plt.figure(figsize=(10,12))
  fig.suptitle('ice diff (n-s)')
  ax1 = plt.subplot(3,1,1)
  plt.plot(body.Time,icemnorth-icemsouth,linestyle='-',color=vplorg)
  ax2 = ax1.twinx()
  plt.plot(body.Time,body.Obliquity,linestyle='--',color=vpllbl)
  
  ax1 = plt.subplot(3,1,2)
  plt.plot(body.Time,icemnorth-icemsouth,linestyle='-',color=vplorg)
  ax2 = ax1.twinx()
  plt.plot(body.Time,body.Eccentricity,linestyle='--',color=vplpur)
  
  
  COPP = body.Eccentricity*np.sin(body.Obliquity*np.pi/180.0) *np.sin((body.ArgP+body.LongA+body.PrecA)*np.pi/180.0) 
  
  ax1 = plt.subplot(3,1,3)
  plt.plot(body.Time,icemnorth-icemsouth,linestyle='-',color=vplorg)
  ax2 = ax1.twinx()
  plt.plot(body.Time,COPP,linestyle='--',color=vpllbl)
  
  
  fig = plt.figure(figsize=(10,12))
  fig.suptitle('ice mass (northern hemi)')
  ax1 = plt.subplot(3,1,1)
  plt.plot(body.Time,icemnorth,linestyle='-',color=vpldbl)
  ax2 = ax1.twinx()
  plt.plot(body.Time,-insol65n,linestyle='--',color=vpllbl)
  plt.ylabel(r'insolation (65$^{\circ}$)')
  
  ax1 = plt.subplot(3,1,2)
  plt.plot(body.Time,icemnorth,linestyle='-',color=vpldbl)
  ax2 = ax1.twinx()
  plt.plot(body.Time,-insol40n,linestyle='--',color=vpllbl)
  plt.ylabel(r'insolation (40$^{\circ}$)')

  ax1 = plt.subplot(3,1,3)
  plt.plot(body.Time,icemnorth,linestyle='-',color=vpldbl)
  ax2 = ax1.twinx()
  plt.plot(body.Time,-insol10n,linestyle='--',color=vpllbl)
  plt.ylabel(r'insolation (10$^{\circ}$)')


  plt.show()
  

def comp2huybers(plname,dir='.',xrange=False,show=True):
  """
  Creates plots of insolation, temperature, albedo, ice mass,
  and bed rock height over the length of the simulation
  
  Parameters
  ----------
  plname : string
    The name of the planet with .Climate data
    
  Keyword Arguments
  -----------------
  dir : string
    Directory of vplanet simulation (default = '.')
  xrange : float tuple, list, or numpy array
    Range of x-values (time) to restrict plot
    (default = False (no restriction))
  orbit : bool
    Plot orbital data (obliquity, eccentricity, COPP)
    (default = False)
  show : bool
    Show plot in Python (default = True)
  
  Output
  ------
  PDF format plot with name 'evol_<dir>.pdf'
  
  """
  if not isinstance(dir,(list,tuple)):
    dir = [dir]
  
  nfiles = len(dir)

  if nfiles > 1 and orbit == True:
    raise Exception("Error: cannot plot multiple files when orbit = True")


  fig = plt.figure(figsize=(8,12))

  fig.subplots_adjust(wspace=0.3,top=0.9,hspace=0.3)

  for ii in np.arange(nfiles):
    out = vplot.GetOutput(dir[ii])
    
    #pdb.set_trace()
  
    ctmp = 0
    for p in range(len(out.bodies)):
      if out.bodies[p].name == plname:
        body = out.bodies[p]
        ctmp = 1
      else:
        if p == len(out.bodies)-1 and ctmp == 0:
          raise Exception("Planet %s not found in folder %s"%(plname,dir[ii]))
  
    try:
      ecc = body.Eccentricity
    except:
      ecc = np.zeros_like(body.Time)+getattr(out.log.initial,plname).Eccentricity
    
    try:
      inc = body.Inc
    except:
      inc = np.zeros_like(body.Time)
    
    try:
      obl = body.Obliquity
    except:
      obltmp = getattr(out.log.initial,plname).Obliquity
      if obltmp.unit == 'rad':
        obltmp *= 180/np.pi
      obl = np.zeros_like(body.Time)+obltmp

    f = open(dir[ii]+'/'+plname+'.in','r')
    lines = f.readlines()
    f.close()
    pco2 = 0
    #pdb.set_trace()
    for i in range(len(lines)):
      if lines[i].split() != []:
        if lines[i].split()[0] == 'dRotPeriod':
          P = -1*np.float(lines[i].split()[1]) 
        if lines[i].split()[0] == 'dSemi':
          semi = np.float(lines[i].split()[1]) 
          if semi < 0:
            semi *= -1
        if lines[i].split()[0] == 'dpCO2':
          pco2 = np.float(lines[i].split()[1])

    try:
      longp = (body.ArgP + body.LongA + body.PrecA)*np.pi/180.0
    except:
      longp = body.PrecA*np.pi/180.0
    
    esinv = ecc*np.sin(longp)#*np.sin(obl*np.pi/180.)
    
    # titlestr = []
#     titlestr.append(r'$a = %f, pCO_2 = %f$'%(semi,pco2))
#     titlestr.append(r'$e_0 = %f, i_0 = %f^{\circ}, \psi_0 = %f^{\circ}, P_{rot} = %f$ d'%(ecc[0],inc[0],obl[0],P))
    fig.subplots_adjust(wspace=0.3,hspace=0.4)

    lats = np.unique(body.Latitude)
    nlats = len(lats)
    ntimes = len(body.Time)
    
    # plot temperature
    temp = np.reshape(body.TempLandLat,(ntimes,nlats))

    ax1 = plt.subplot(6,1,4)
    pos = ax1.figbox.get_points()
    c = plt.contourf(body.Time,lats[lats>60],temp.T[lats>60],20,cmap='jet')
    plt.ylabel('Latitude')
    plt.title(r'Surface Temp ($^{\circ}$C)',fontsize=12)
    plt.ylim(60,85)
    plt.yticks([60,70,80])
    if xrange == False:
      left = 0
    else:
      left = xrange[0]
#     plt.text(left,140,'\n'.join(titlestr),fontsize=20) 
    if xrange:
      plt.xlim(xrange)
    plt.contour(body.Time,lats[lats>60],temp.T[lats>60],levels=[0],colors='w')

    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
# plot ice accumulation
#     alb = np.reshape(body.AlbedoLat,(ntimes,nlats))
#     if orbit == True:
#       ax2 = plt.subplot(4,2,3)
#     else:
#       ax2 = plt.subplot(5,nfiles,ii+2*nfiles+1)
#     pos = ax2.figbox.get_points()
#     c = plt.contourf(body.Time,lats,alb.T,cmap='Blues_r')
#     plt.ylabel('Latitude')
#     plt.title('Albedo (TOA)')
#     plt.ylim(-90,90)
#     plt.yticks([-60,-30,0,30,60])
#     if xrange:
#       plt.xlim(xrange)
#     plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  

    # plot ice height
    ice = np.reshape(body.IceHeight+body.BedrockH,(ntimes,nlats))

    ax3 = plt.subplot(6,1,3)
    pos = ax3.figbox.get_points()
#     pdb.set_trace()
    c = plt.contourf(body.Time,lats[lats>60],ice.T[lats>60,:],20,cmap='jet')
    plt.ylabel('Latitude')
    plt.title('Ice sheet height (m)',fontsize=12)
    plt.ylim(60,85)
  #   plt.xlim(0,2e6)
    plt.yticks([60,70,80])
    if xrange:
      plt.xlim(xrange)
    plt.contour(body.Time,lats[lats>60],ice.T[lats>60],levels=[0],colors='w')

    clb=plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
    #clb.set_label('Ice height (m)')
    # ax3p = ax3.twinx()
  #   plt.plot(body.Time,esinv,linestyle = 'solid',marker='None',color='salmon',linewidth=2)
  

    ax4 = plt.subplot(6,1,5)
    pos = ax4.figbox.get_points()
    acc = body.IceAccum
    c = plt.contourf(body.Time,lats[lats>60],acc.T[lats>60],20,cmap='jet')
    plt.ylabel('Latitude')
    plt.title(r'Ice Accumulation (m year$^{-1}$)',fontsize=12)
    plt.ylim(60,85)
    plt.yticks([60,70,80])
    if xrange == False:
      left = 0
    else:
      left = xrange[0]
#     plt.text(left,140,'\n'.join(titlestr),fontsize=20) 
    if xrange:
      plt.xlim(xrange)
    clb=plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
    tloc = ticker.MaxNLocator(nbins=5)
    clb.locator=tloc
    clb.update_ticks()
    
    ax5 = plt.subplot(6,1,6)
    pos = ax5.figbox.get_points()
    abl = body.IceAblate
    c = plt.contourf(body.Time,lats[lats>60],-abl.T[lats>60],20,cmap='jet')
    plt.ylabel('Latitude')
    plt.title(r'Ice Ablation (m year$^{-1}$)',fontsize=12)
    plt.ylim(60,85)
    plt.yticks([60,70,80])
    plt.xlabel('Time [years]')
    if xrange == False:
      left = 0
    else:
      left = xrange[0]
#     plt.text(left,140,'\n'.join(titlestr),fontsize=20) 
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
# plot bedrock
#     brock = np.reshape(body.BedrockH,(ntimes,nlats))
#     if orbit == True:
#       ax4 = plt.subplot(4,2,7)
#     else:
#       ax4 = plt.subplot(5,nfiles,ii+4*nfiles+1)
#     pos = ax4.figbox.get_points()
#     c = plt.contourf(body.Time,lats,brock.T,cmap='Reds_r')
#     plt.ylabel('Latitude')
#     plt.title('Bedrock height [m]')
#     plt.ylim(-90,90)
#     plt.yticks([-60,-30,0,30,60])
#     plt.xlabel('Time [years]')
#     if xrange:
#       plt.xlim(xrange)
#     plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  

    # plot insolation
#     insol = np.reshape(body.AnnInsol,(ntimes,nlats))
#     if orbit == True:
#       ax5 = plt.subplot(4,2,2)
#     else:
#       ax5 = plt.subplot(5,nfiles,ii+nfiles+1)
#     pos = ax5.figbox.get_points()
#     c = plt.contourf(body.Time,lats,insol.T,cmap='plasma')
#     plt.ylabel('Latitude')
#     plt.title(r'Annual average insolation [w/m$^2$]')
#     plt.ylim(-90,90)
#     plt.yticks([-60,-30,0,30,60])
#     if xrange:
#       plt.xlim(xrange)
#     plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))

#     if orbit == True:
      #obliquity
    plt.subplot(6,1,2)
    plt.plot(body.Time,obl,linestyle = 'solid',marker='None',color='darkblue',linewidth =2)
    plt.ylabel('Obliquity')
    if xrange:
      plt.xlim(xrange)

    #eccentricity
  #   plt.subplot(4,2,6)
#       plt.plot(body.Time,ecc,linestyle = 'solid',marker='None',color='darkorchid',linewidth =2)
#       plt.ylabel('Eccentricity')
#       if xrange:
#         plt.xlim(xrange)

    #e sin(obl) sin varpi
    plt.subplot(6,1,1)
    plt.plot(body.Time,esinv,linestyle = 'solid',marker='None',color='salmon',linewidth=2)
    plt.ylabel('CPP')

    if xrange:
      plt.xlim(xrange)

    if dir[ii] == '.':
      dir[ii] = 'cwd'
  
  #fig.suptitle('\n'.join(titlestr),fontsize=20) 
  
  if xrange:
    sfile = 'comp2huy_'+'_'.join(dir)+'_%d_%d.pdf'%(xrange[0],xrange[1])
  else:
    sfile = 'comp2huy_'+'_'.join(dir)+'.pdf'
  plt.savefig(sfile)

  if show:
    plt.show()
  else:
    plt.close()

def globavg(plname,dir='.',xrange=False,orbit=False,show=True):

  if not isinstance(dir,(list,tuple)):
    dir = [dir]
  
  nfiles = len(dir)

  if nfiles > 1 and orbit == True:
    raise Exception("Error: cannot plot multiple files when orbit = True")
  
  

#   if orbit == True:
#     fig = plt.figure(figsize=(16,13))
#   else:
#     fig = plt.figure(figsize=(10*nfiles,15))
# 
#   fig.subplots_adjust(wspace=0.3,top=0.9)

  for ii in np.arange(nfiles):
    out = vplot.GetOutput(dir[ii])
    
    #pdb.set_trace()
  
    ctmp = 0
    for p in range(len(out.bodies)):
      if out.bodies[p].name == plname:
        body = out.bodies[p]
        ctmp = 1
      else:
        if p == len(out.bodies)-1 and ctmp == 0:
          raise Exception("Planet %s not found in folder %s"%(plname,dir[ii]))
  
    try:
      ecc = body.Eccentricity
    except:
      ecc = np.zeros_like(body.Time)+getattr(out.log.initial,plname).Eccentricity
    
    try:
      inc = body.Inc
    except:
      inc = np.zeros_like(body.Time)
    
    try:
      obl = body.Obliquity
    except:
      obltmp = getattr(out.log.initial,plname).Obliquity
      if obltmp.unit == 'rad':
        obltmp *= 180/np.pi
      obl = np.zeros_like(body.Time)+obltmp

    f = open(dir[ii]+'/'+plname+'.in','r')
    lines = f.readlines()
    f.close()
    pco2 = 0
    #pdb.set_trace()
    for i in range(len(lines)):
      if lines[i].split() != []:
        if lines[i].split()[0] == 'dRotPeriod':
          P = -1*np.float(lines[i].split()[1]) 
        if lines[i].split()[0] == 'dSemi':
          semi = np.float(lines[i].split()[1]) 
          if semi < 0:
            semi *= -1
        if lines[i].split()[0] == 'dpCO2':
          pco2 = np.float(lines[i].split()[1])

    try:
      longp = (body.ArgP + body.LongA + body.PrecA)*np.pi/180.0
    except:
      longp = body.PrecA*np.pi/180.0
    
    esinv = ecc*np.sin(longp)*np.sin(obl*np.pi/180.)
    
    plt.figure(figsize=(8,8))
    plt.subplot(2,1,1)
    plt.plot(body.Time,body.TGlobal)
    plt.ylabel('Global average temperature')
    
    plt.subplot(2,1,2)
    plt.plot(body.Time,body.AlbedoGlobal)
    plt.ylabel('Global average albedo')
    plt.xlabel('Time (year)')
    
    plt.savefig('global_'+dir[ii]+'.pdf')
    if show == True:
      plt.show()
    else:
      plt.close()
 
# def init_icelines(line):
#   line.set_data([],[])
#   line1.set_data([],[])
#   line2.set_data([],[])
#   return line, line1, line2
#   
def ice_curves(i,Time,Obliquity,Eccentricity,Latitude,IceMass,TempMaxWater,phi,delta,alpha,qref,line,line0,line1,line2,line3,line4,time,obl,ecc,line5,line6,line7,line8,line9,line10,dir):
  beta = Obliquity[i]
  s20  = -5./16*(2-3*(np.sin(beta*np.pi/180))**2)

  qs = ebm.q(np.sin(phi*np.pi/180),delta,s20,alpha)
  #pdb.set_trace()
  #find ice line
  icenorth = IceMass[i,Latitude>=0]
  phinorth = Latitude[Latitude>=0]
  tempwnorth = TempMaxWater[i,Latitude>=0]
  if (icenorth==0).all():
    phi_ice = 90.0
  else:
    phi_ice = (phinorth[np.where(icenorth>0)])
  
  if (tempwnorth>-2).all():
    phi_sea = 90.0
  else:
    phi_sea = (phinorth[np.where(tempwnorth<=-2)])
    
  icesouth = IceMass[i,Latitude<=0]
  phisouth = Latitude[Latitude<=0]
  tempwsouth = TempMaxWater[i,Latitude<=0]
  if (icesouth==0).all():
    phi_ice2 = 90.0
  else:
   #  if icesouth[-1] > 0 and icesouth[0] == 0:
#       #ice cap
    phi_ice2 = np.abs((phisouth[np.where(icesouth>0)]))
    # elif icesouth[-1] == 0 and icesouth[0] > 0:
#       #ice belt
#     
#     elif   
    
  
  if (tempwsouth>-2).all():
    phi_sea2 = 90.0
  else:
    phi_sea2 = np.abs((phisouth[np.where(tempwsouth<=-2)]))
  
#   pdb.set_trace()
  qcurr = qref/np.sqrt(1-Eccentricity[i]**2)
  plt.subplot(121)
  if s20 < 0:
    line.set_data(qs,phi)
    line0.set_data([],[])
  else: 
    line.set_data([],[])
    line0.set_data(qs,phi)
#   line1.set_data(qcurr,phi_ice)
#   line2.set_data(qcurr,phi_sea)
#   line3.set_data(qcurr,phi_ice2)
#   line4.set_data(qcurr,phi_sea2)
  time.set_text('%#.0f years'%(Time[i]))
  obl.set_text(r'$\psi=%#.2f^{\circ}$'%(Obliquity[i]))
  #ecc.set_text(r'$e=%#.5f$'%(Eccentricity[i]))
  plt.subplot(122)
  if s20 > 0:
    line5.set_data(qs,phi)
    line10.set_data([],[])
  else:
    line5.set_data([],[])
    line10.set_data(qs,phi)
#   line6.set_data(qcurr,phi_ice)
#   line7.set_data(qcurr,phi_sea)
#   line8.set_data(qcurr,phi_ice2)
#   line9.set_data(qcurr,phi_sea2)

  plt.savefig('ice_images_'+dir+'/%05d.png'%i)
  return line,line0, line1, line2, line3, line4, time,obl,ecc,line5,line6,line7,line8,line9,line10
      
def ice_stab(plname,dir='.'):
  if not isinstance(dir,(list,tuple)):
    dir = [dir]

  nfiles = len(dir)

  # if nfiles > 1 and orbit == True:
#     raise Exception("Error: cannot plot multiple files when orbit = True")
  
  hdffile = 0
  if os.path.exists('climatedata.hdf5'):
    hdff = h5py.File('climatedata.hdf5','r')
    hdffile = 1
  
  for ii in np.arange(nfiles):
    if hdffile:
      Time = hdff[dir[ii]]['Time'][...]
      Latitude = np.unique(hdff[dir[ii]]['Latitude'][...])
      Obliquity = hdff[dir[ii]]['Obliquity'][...]
      Eccentricity = hdff[dir[ii]]['Eccentricity'][...]
      #pdb.set_trace()
      IceMass = np.reshape(hdff[dir[ii]]['IceMass'][...],(len(Time),len(Latitude)))
      TempMaxWater = np.reshape(hdff[dir[ii]]['TempMaxWater'][...],(len(Time),len(Latitude)))
    else:
      out = vplot.GetOutput(dir[ii])
    
    #pdb.set_trace()
  
      ctmp = 0
      for p in range(len(out.bodies)):
        if out.bodies[p].name == plname:
          body = out.bodies[p]
          ctmp = 1
        else:
          if p == len(out.bodies)-1 and ctmp == 0:
            raise Exception("Planet %s not found in folder %s"%(plname,dir[ii]))
      
      Time = body.Time
      Obliquity = body.Obliquity
      Eccentricity = body.Eccentricity
      Latitude = body.Latitude[0]
      IceMass = body.IceMass
      TempMaxWater = body.TempMaxWater
    
  f = open(dir[ii]+'/star.in','r')
  fl = f.readlines()
  f.close()
  for jj in np.arange(len(fl)):
    if fl[jj].split() != []:
      if fl[jj].split()[0] == 'dLuminosity':  
        lum = np.float(fl[jj].split()[1])
    
  #pdb.set_trace()  
  #calc some shiz
  phi = np.linspace(0.1,89.9,100)
  xs = np.sin(phi*np.pi/180)
  #coalbedos
  al = (1-0.363)
  aw = (1-0.263)
  ai = (1-0.6)
  a0 = 0.34*al+0.66*aw
  alpha = (a0-ai)/a0

  #diffusion parameter
  D = 0.58
  B = 2.09
  delta = D/B
  
  Tref = -2.0
  Q = lum/(4*np.pi*(1.0031*1.49598e11)**2)/4
  A = 203.3
  olr = A+B*Tref
  qref = a0*Q/olr  #need to divide by sqrt(1-e**2)
  
  #set up plot
  fig = plt.figure(figsize=(14,6))
  ax = fig.add_subplot(121)
  plt.xlim(1.1,1.6)
  plt.ylim(0,90)
  line, = ax.plot([],[],'-',color='k')
  line0, = ax.plot([],[],'-',color='r')
  line1, = ax.plot([],[],'o',linestyle='None',mfc='None',mec='k',ms=5,mew=2)
  line2, = ax.plot([],[],'o',linestyle='None',mfc='None',mec='b',ms=5,mew=2)
  line3, = ax.plot([],[],'^',linestyle='None',mfc='None',mec='k',ms=5,mew=2)
  line4, = ax.plot([],[],'^',linestyle='None',mfc='None',mec='b',ms=5,mew=2)
  time = ax.text(1.4,70,'',fontsize=20)
  obl = ax.text(1.4,60,'',fontsize=20)
  ecc = ax.text(1.4,50,'',fontsize=20)
  plt.xlabel('$q$',fontsize=20)
  plt.ylabel('Ice edge latitude (degrees)',fontsize=20)
  
  ax = fig.add_subplot(122)
  plt.xlim(1.1,1.6)
  plt.ylim(90,0,-1)
  line5, = ax.plot([],[],'-',color='k')
  line10, = ax.plot([],[],'-',color='r')
  line6, = ax.plot([],[],'o',linestyle='None',mfc='None',mec='k',ms=5,mew=2)
  line7, = ax.plot([],[],'o',linestyle='None',mfc='None',mec='b',ms=5,mew=2)
  line8, = ax.plot([],[],'^',linestyle='None',mfc='None',mec='k',ms=5,mew=2)
  line9, = ax.plot([],[],'^',linestyle='None',mfc='None',mec='b',ms=5,mew=2)
  plt.xlabel('$q$',fontsize=20)

  if not os.path.exists('ice_images_'+dir[0]):
    os.mkdir('ice_images_'+dir[0])
  anim = animation.FuncAnimation(fig, ice_curves, np.arange(0,len(Time)),interval=50,blit=True,fargs=(Time,Obliquity,Eccentricity,Latitude,IceMass,TempMaxWater,phi,delta,alpha,qref,line,line0,line1,line2,line3,line4,time,obl,ecc,line5,line6,line7,line8,line9,line10,dir[0]))
  
  if os.path.exists('climatedata.hdf5'):
    hdff.close()
  plt.show()
  
def icestabregion(plname,dir='.',stepsinst=2):
  if not isinstance(dir,(list,tuple)):
    dir = [dir]
  
  nfiles = len(dir)

  # if nfiles > 1 and orbit == True:
#     raise Exception("Error: cannot plot multiple files when orbit = True")
  
  hdffile = 0
  if os.path.exists('climatedata.hdf5'):
    hdff = h5py.File('climatedata.hdf5','r')
    hdffile = 1
  
  for ii in np.arange(nfiles):
    if hdffile:
      Time = hdff[dir[ii]]['Time'][...]
      Latitude = np.unique(hdff[dir[ii]]['Latitude'][...])
      Obliquity = hdff[dir[ii]]['Obliquity'][...]
      Eccentricity = hdff[dir[ii]]['Eccentricity'][...]
      #pdb.set_trace()
      IceMass = np.reshape(hdff[dir[ii]]['IceMass'][...],(len(Time),len(Latitude)))
      TempMaxWater = np.reshape(hdff[dir[ii]]['TempMaxWater'][...],(len(Time),len(Latitude)))
      Snow = hdff[dir[ii]]['Snowball'][...]
    else:
      out = vplot.GetOutput(dir[ii])
    
    #pdb.set_trace()
  
      ctmp = 0
      for p in range(len(out.bodies)):
        if out.bodies[p].name == plname:
          body = out.bodies[p]
          ctmp = 1
        else:
          if p == len(out.bodies)-1 and ctmp == 0:
            raise Exception("Planet %s not found in folder %s"%(plname,dir[ii]))
      
      Time = body.Time
      Obliquity = body.Obliquity
      Eccentricity = body.Eccentricity
      Latitude = body.Latitude[0]
      IceMass = body.IceMass
      TempMaxWater = body.TempMaxWater
      Snow = body.Snowball
    
  #pdb.set_trace()  
  #calc some shiz
  phi = np.linspace(0.1,89.9,100)
  xs = np.sin(phi*np.pi/180)
  #coalbedos
  al = (1-0.363)
  aw = (1-0.263)
  ai = (1-0.6)
  a0 = 0.34*al+0.66*aw
  alpha = (a0-ai)/a0

  #diffusion parameter
  D = 0.58
  B = 2.09
  delta = D/B
  
  Tref = -2.0
  Q = 3.77e26/(4*np.pi*(1.0031*1.49598e11)**2)/4
  A = 203.3
  olr = A+B*Tref
  qref = a0*Q/olr  #need to divide by sqrt(1-e**2)
  
  beta = Obliquity  
  s20_max  = -5./16*(2-3*(np.sin(np.max(beta)*np.pi/180))**2)
  s20_min  = -5./16*(2-3*(np.sin(np.min(beta)*np.pi/180))**2)
  e2 = Eccentricity[beta == np.max(beta)]
  q2 = qref/np.sqrt(1.0-e2**2)

  e1 = Eccentricity[beta == np.min(beta)] 
  q1 = qref/np.sqrt(1.0-e1**2)
  
  fig = plt.figure(figsize=(7,6))
  plt.subplot(1,1,1)
  
  plt.ylim(-1,91)
  
  if s20_max > 0:
    s20_max = -1e-3
    
  qs_max = ebm.q(np.sin(phi*np.pi/180),delta,s20_max,alpha)
  qs_min = ebm.q(np.sin(phi*np.pi/180),delta,s20_min,alpha)
  
  plt.xlim(0.99*np.min([np.min(qs_max),np.min(qs_min)]),np.max([np.max(qs_max),np.max(qs_min)]))

#     qcurr = qref/np.sqrt(1-Eccentricity[jj]**2)
  
  plt.plot(qs_max,phi,c=vplred,alpha=0.9)
  plt.plot(qs_min,phi,c=vpllbl,alpha=0.9)
  plt.hlines(0,0.99*np.min([np.min(qs_max),np.min(qs_min)]), qs_max[0],color=vplred,alpha=0.9,lw=4)
  plt.hlines(0,0.99*np.min([np.min(qs_max),np.min(qs_min)]), qs_min[0],color=vpllbl,alpha=0.9,zorder=100)
  
  plt.hlines(90,qs_max[-1],np.max([np.max(qs_max),np.max(qs_min)]), color=vplred,alpha=0.9,lw=4)
  plt.hlines(90,qs_min[-1],np.max([np.max(qs_max),np.max(qs_min)]), color=vpllbl,alpha=0.9,zorder=100)
  
  plt.vlines(q1,0,90,color=vpllbl,alpha=0.9,linestyle='--')
  plt.vlines(q2,0,90,color=vplred,alpha=0.9,linestyle='--')
  plt.fill_betweenx(phi,qs_min,qs_max,facecolor='0.5',alpha=0.5)
  
  if (Snow==1).any():
    snowtime = np.max(np.where(Snow==0)[0])-stepsinst
    colorp = vpldbl    
    s20_inst = -5./16*(2-3*(np.sin(np.max(beta[snowtime])*np.pi/180))**2)
    qs_inst = ebm.q(np.sin(phi*np.pi/180),delta,s20_inst,alpha)
    plt.plot(qs_inst,phi,c=vpldbl,alpha=1,lw=3)
    plt.hlines(0,0.99*np.min([np.min(qs_max),np.min(qs_min)]), qs_inst[0],color=vpldbl,alpha=0.9,lw=3)
    plt.hlines(90,qs_inst[-1],np.max([np.max(qs_max),np.max(qs_min)]), color=vpldbl,alpha=0.9,lw=3)

  
  else:
    snowtime = np.where(beta==np.max(beta))[0]
    colorp = vplred

  q_inst = qref/(np.sqrt(1.0-Eccentricity[snowtime]**2))

  icenorth = IceMass[snowtime,Latitude>=0]
  phinorth = Latitude[Latitude>=0]
  tempwnorth = TempMaxWater[snowtime,Latitude>=0]
  if (icenorth==0).all():
    phi_ice = 90.0
  else:
    phi_ice = np.min(phinorth[np.where(icenorth>0)])

  if (tempwnorth>-2).all():
    phi_sea = 90.0
  else:
    phi_sea = np.min(phinorth[np.where(tempwnorth<=-2)])

  icesouth = IceMass[snowtime,Latitude<=0]
  phisouth = Latitude[Latitude<=0]
  tempwsouth = TempMaxWater[snowtime,Latitude<=0]
  if (icesouth==0).all():
    phi_ice2 = 90.0
  else:
    phi_ice2 = np.abs(np.max(phisouth[np.where(icesouth>0)]))

  if (tempwsouth>-2).all():
    phi_sea2 = 90.0
  else:
    phi_sea2 = np.abs(np.max(phisouth[np.where(tempwsouth<=-2)]))

  plt.plot(q_inst,phi_sea,'o',ms=10,mfc=colorp,mec=colorp,mew=3,label='Ocean (N)')
  plt.plot(q_inst,phi_sea2,'o',ms=10,mfc='None',mec=colorp,mew=3,label='Ocean (S)')
  plt.plot(q_inst,phi_ice,'^',ms=10,mfc=colorp,mec=colorp,mew=3,label='Land (N)')
  plt.plot(q_inst,phi_ice2,'^',ms=10,mfc='None',mec=colorp,mew=3,label='Land (S)')
  plt.legend(loc='upper right',numpoints=1)
  plt.ylabel('Ice edge latitude')
  plt.xlabel(r'$q$')
  
  plt.savefig('ice_stab_'+dir[ii]+'.pdf')
#   plt.show()
  
  #-----time evolution---------------------
  plt.figure(figsize=(8,10))
  plt.subplot(3,1,1)
  plt.plot(Time,Obliquity,'k-')
  plt.ylabel('Obliquity (deg)')
  
  plt.subplot(3,1,2)
  plt.plot(Time,qref/(np.sqrt(1-Eccentricity**2)),'k-')
  plt.ylabel('q')
  
  icenorth = IceMass[:,Latitude>=0]
  tempwnorth = TempMaxWater[:,Latitude>=0]
  phi_ice = np.zeros_like(Time)
  phi_sea = np.zeros_like(Time)
  
  icesouth = IceMass[:,Latitude<=0]
  tempwsouth = TempMaxWater[:,Latitude<=0]
  phi_ice2 = np.zeros_like(Time)
  phi_sea2 = np.zeros_like(Time)
  
  for jj in np.arange(len(Time)):
    if (icenorth[jj]==0).all():
      phi_ice[jj] = 90.0
    else:
      phi_ice[jj] = np.min(phinorth[np.where(icenorth[jj]>0)])

    if (tempwnorth[jj]>-2).all():
      phi_sea[jj] = 90.0
    else:
      phi_sea[jj] = np.min(phinorth[np.where(tempwnorth[jj]<=-2)])

  
    if (icesouth[jj]==0).all():
      phi_ice2[jj] = 90.0
    else:
      phi_ice2[jj] = np.abs(np.max(phisouth[np.where(icesouth[jj]>0)]))

    if (tempwsouth[jj]>-2).all():
      phi_sea2[jj] = 90.0
    else:
      phi_sea2[jj] = np.abs(np.max(phisouth[np.where(tempwsouth[jj]<=-2)]))  
  
  plt.subplot(3,1,3)
  plt.plot(Time,phi_ice,'-',color=vplred,label='Land (N)')
  plt.plot(Time,phi_ice2,'-',color=vplorg,label='Land (S)') 
  
  plt.plot(Time,phi_sea,'-',color=vpldbl,label='Ocean (N)')
  plt.plot(Time,phi_sea2,'-',color=vpllbl,label='Ocean (S)')   
  plt.ylabel('Ice edge (deg)')
  plt.legend(loc='upper right',fontsize=10)
  plt.savefig('ice_edge_q_'+dir[ii]+'.pdf')    
  
  #-------parametric plot dq/dx vs dq----------
  
  s20t = -5./16*(2-3*(np.sin((beta)*np.pi/180))**2)
  phi_ice[phi_ice==0] = 0.1
  phi_ice2[phi_ice2==0] = 0.1
  qNL = ebm.q(np.sin(phi_ice*np.pi/180),delta,s20t,alpha)
  qSL = ebm.q(np.sin(phi_ice2*np.pi/180),delta,s20t,alpha)
  
  phi_sea[phi_sea==0] = 0.1
  phi_sea2[phi_sea2==0] = 0.1
  qNW = ebm.q(np.sin(phi_sea*np.pi/180),delta,s20t,alpha)
  qSW = ebm.q(np.sin(phi_sea2*np.pi/180),delta,s20t,alpha) 
  
  qtime = qref/(np.sqrt(1-Eccentricity**2))
  
  dqNL = qtime-qNL
  dqSL = qtime-qSL
  dqNW = qtime-qNW
  dqSW = qtime-qSW
     
  dqdxNL = -qNL**2*ebm.dqinvdx(np.sin(phi_ice*np.pi/180),delta,s20t,alpha)
  dqdxSL = -qSL**2*ebm.dqinvdx(np.sin(phi_ice2*np.pi/180),delta,s20t,alpha)
    
  dqdxNW = -qNW**2*ebm.dqinvdx(np.sin(phi_sea*np.pi/180),delta,s20t,alpha)
  dqdxSW = -qSW**2*ebm.dqinvdx(np.sin(phi_sea2*np.pi/180),delta,s20t,alpha)
  
  if (s20t>=0).any():
    dqdxNL[s20t>=0] *= -1
    dqdxSL[s20t>=0] *= -1
    dqdxNW[s20t>=0] *= -1
    dqdxSW[s20t>=0] *= -1
  
#   pdb.set_trace()    
  
  if (Snow==1).any():
    snowt = np.min(np.where(Snow==1)[0])
  else: 
    snowt = len(Time)
  
  plt.figure(figsize=(8,8))
  plt.plot(dqNL[:snowt],dqdxNL[:snowt],'^',color=vplred,alpha=0.5)
  plt.plot(dqNL[snowt-1],dqdxNL[snowt-1],'o',color=vplred,ms=10)
  plt.plot([0,0],[-1,1],'--',color='0.5')
  plt.plot([-1,1],[0,0],'--',color='0.5')
  plt.xlim(-0.1,0.1)
  plt.ylim(-.5,.5)
  plt.xlabel(r'$\Delta q$',fontsize=18)
  plt.ylabel(r'$dq/dx_s$',fontsize=18)
  
#   plt.figure(figsize=(8,8))
  plt.plot(dqSL[:snowt],dqdxSL[:snowt],'^',color=vplorg,alpha=0.5)
  #plt.plot(dqSL[snowt-1],dqdxSL[snowt-1],'o',color=vplorg,ms=10)
#   plt.plot([0,0],[-1,1],'--',color='0.5')
#   plt.plot([-1,1],[0,0],'--',color='0.5')
#   plt.xlim(-0.1,0.1)
#   plt.ylim(-.5,.5)
  
#   plt.figure(figsize=(8,8))
  plt.plot(dqNW[:snowt],dqdxNW[:snowt],'o',color=vpldbl,alpha=0.5)
  #plt.plot(dqNW[snowt-1],dqdxNW[snowt-1],'o',color=vpldbl,ms=10)
#   plt.plot([0,0],[-1,1],'--',color='0.5')
#   plt.plot([-1,1],[0,0],'--',color='0.5')
#   plt.xlim(-0.1,0.1)
#   plt.ylim(-.5,.5)
  
#   plt.figure(figsize=(8,8))
  plt.plot(dqSW[:snowt],dqdxSW[:snowt],'o',color=vpllbl,alpha=0.5)
  #plt.plot(dqSW[snowt-1],dqdxSW[snowt-1],'o',color=vpllbl,ms=10)
#   plt.plot([0,0],[-1,1],'--',color='0.5')
#   plt.plot([-1,1],[0,0],'--',color='0.5')
#   plt.xlim(-0.1,0.1)
#   plt.ylim(-.5,.5)
  plt.savefig('parametric_ice_'+dir[ii]+'.pdf')

  if Time[:snowt][-1] >= 1e6:
      Time /= 1e6
      tunit = ' (Myr)'
      # if xrange:
#         xrange /= 1e6
  else:
      tunit = ' (years)'

  fig = plt.figure(figsize=(8,8))
  fig.subplots_adjust(hspace=0.15)
  plt.subplot(2,1,1)
  plt.plot(Time[:snowt],dqdxNL[:snowt],'-',color=vplred)
  plt.plot(Time[:snowt],dqdxSL[:snowt],'-',color=vplorg)
  plt.plot(Time[:snowt],dqdxNW[:snowt],'-',color=vpldbl)
  plt.plot(Time[:snowt],dqdxSW[:snowt],'-',color=vpllbl)
  plt.plot(Time[:snowt],np.zeros_like(Time[:snowt]),'--',color='0.5')
  plt.ylabel(r'$dq/dx_s$',fontsize=18)
  plt.ylim(-.5,.5)
  
  plt.subplot(2,1,2)
  plt.plot(Time[:snowt],dqNL[:snowt],'-',color=vplred)
  plt.plot(Time[:snowt],dqSL[:snowt],'-',color=vplorg)
  plt.plot(Time[:snowt],dqNW[:snowt],'-',color=vpldbl)
  plt.plot(Time[:snowt],dqSW[:snowt],'-',color=vpllbl)
  plt.plot(Time[:snowt],np.zeros_like(Time[:snowt]),'--',color='0.5')
  plt.ylabel(r'$\Delta q$',fontsize=18)
  plt.xlabel('Time'+tunit,fontsize=18)
  plt.ylim(-.1,.1)
  plt.savefig('stability_time_'+dir[ii]+'.pdf')
  plt.show()
  
def peak_insol(plname,dir='.',xrange=False,orbit=False,show=True):
  """
  Creates plots of insolation, temperature, albedo, ice mass,
  and bed rock height over the length of the simulation
  
  Parameters
  ----------
  plname : string
    The name of the planet with .Climate data
    
  Keyword Arguments
  -----------------
  dir : string
    Directory of vplanet simulation (default = '.')
  xrange : float tuple, list, or numpy array
    Range of x-values (time) to restrict plot
    (default = False (no restriction))
  orbit : bool
    Plot orbital data (obliquity, eccentricity, COPP)
    (default = False)
  show : bool
    Show plot in Python (default = True)
  
  Output
  ------
  PDF format plot with name 'evol_<dir>.pdf'
  
  """
  if not isinstance(dir,(list,tuple)):
    dir = [dir]
  
  nfiles = len(dir)

  if nfiles > 1 and orbit == True:
    raise Exception("Error: cannot plot multiple files when orbit = True")
  
  if orbit == True:
    fig = plt.figure(figsize=(10*nfiles,10))
  else:
    fig = plt.figure(figsize=(8*nfiles,7))

  fig.subplots_adjust(wspace=0.3,top=0.9)

  hdffile = 0
  if os.path.exists('climatedata.hdf5'):
    hdff = h5py.File('climatedata.hdf5')
    hdffile = 1

  for ii in np.arange(nfiles):
    if hdffile:
      time = hdff[dir[ii]]['Time'][...]
      try:
        ecc = hdff[dir[ii]]['Eccentricity'][...]
      except:
        ecc = np.zeros_like(time)+hdff[dir[ii]]['EccInit'][...]
        
      try:
        inc = (hdff[dir[ii]]['Inc'][...])
      except:
        inc = np.zeros_like(time)
        
      try:  
        obl = (hdff[dir[ii]]['Obliquity'][...])
      except:
        obl = np.zeros_like(time)+hdff[dir[ii]]['OblInit'][...]
        
      preca = hdff[dir[ii]]['PrecA'][...]
      try:
        argp = (hdff[dir[ii]]['ArgP'][...])
        longa = hdff[dir[ii]]['LongA'][...]
        longp = (argp+longa+preca)*np.pi/180.0
      except:
        longp = preca*np.pi/180.0
        
      lats = np.unique(hdff[dir[ii]]['Latitude'][...])
      TempLat = hdff[dir[ii]]['TempLat'][...]
      AlbedoLat = hdff[dir[ii]]['AlbedoLat'][...]
      IceHeight = hdff[dir[ii]]['IceHeight'][...]
      BedrockH = hdff[dir[ii]]['BedrockH'][...]
      AnnInsol = hdff[dir[ii]]['AnnInsol'][...]
      
    else:
      out = vplot.GetOutput(dir[ii])
    
    #pdb.set_trace()
  
      ctmp = 0
      for p in range(len(out.bodies)):
        if out.bodies[p].name == plname:
          body = out.bodies[p]
          ctmp = 1
        else:
          if p == len(out.bodies)-1 and ctmp == 0:
            raise Exception("Planet %s not found in folder %s"%(plname,dir[ii]))
  
      try:
        ecc = body.Eccentricity
      except:
        ecc = np.zeros_like(body.Time)+getattr(out.log.initial,plname).Eccentricity
    
      try:
        inc = body.Inc
      except:
        inc = np.zeros_like(body.Time)
    
      try:
        obl = body.Obliquity
      except:
        obltmp = getattr(out.log.initial,plname).Obliquity
        if obltmp.unit == 'rad':
          obltmp *= 180/np.pi
        obl = np.zeros_like(body.Time)+obltmp

      try:
        longp = (body.ArgP + body.LongA + body.PrecA)*np.pi/180.0
      except:
        longp = body.PrecA*np.pi/180.0
    
      lats = np.unique(body.Latitude)
      time = body.Time
      TempLat = body.TempLat
      AlbedoLat = body.AlbedoLat
      IceHeight = body.IceHeight
      BedrockH = body.BedrockH
      PeakInsol = body.PeakInsol
    
    f = open(dir[ii]+'/'+plname+'.in','r')
    lines = f.readlines()
    f.close()
    pco2 = 0
    #pdb.set_trace()
    for i in range(len(lines)):
      if lines[i].split() != []:
        if lines[i].split()[0] == 'dRotPeriod':
          P = -1*np.float(lines[i].split()[1]) 
        if lines[i].split()[0] == 'dSemi':
          semi = np.float(lines[i].split()[1]) 
          if semi < 0:
            semi *= -1
        if lines[i].split()[0] == 'dpCO2':
          pco2 = np.float(lines[i].split()[1])
  
    esinv = ecc*np.sin(longp)*np.sin(obl*np.pi/180.)
    
    titlestr = []
    titlestr.append(r'$a = %f, pCO_2 = %f$'%(semi,pco2))
    titlestr.append(r'$e_0 = %f, i_0 = %f^{\circ}, \psi_0 = %f^{\circ}, P_{rot} = %f$ d'%(ecc[0],inc[0],obl[0],P))
    fig.subplots_adjust(wspace=0.3)

    nlats = len(lats)
    ntimes = len(time)

    # plot insolation
    insol = np.reshape(PeakInsol,(ntimes,nlats))
    if orbit == True:
      ax5 = plt.subplot(4,1,1)
    else:
      ax5 = plt.subplot(1,nfiles,ii+nfiles)
    pos = ax5.figbox.get_points()
    c = plt.contourf(time/1e6,lats,insol.T,cmap='plasma',rasterized=True)
    plt.ylabel('Latitude ($^{\circ}$)',fontsize=22,fontweight='bold')
    plt.xlabel('Time (Myr)',fontsize=22,fontweight='bold')
    plt.title(r'Peak insolation [W m$^{-2}$]',fontsize=22,fontweight='bold')
    plt.ylim(np.min(lats),np.max(lats))
    plt.yticks([-60,-30,0,30,60],fontsize=22,fontweight='bold')
    plt.xticks(fontsize=20,fontweight='bold')

    if xrange:
      plt.xlim(xrange)
    clb=plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
    clb.ax.set_yticklabels(['%d'%val for val in np.linspace(200,520,9)],fontsize=20,weight='bold')

    if orbit == True:
      #obliquity
      plt.subplot(4,1,2)
      plt.plot(time,obl,linestyle = 'solid',marker='None',color='darkblue',linewidth =2,rasterized=True)
      plt.ylabel('Obliquity')
      if xrange:
        plt.xlim(xrange)

      #eccentricity
      plt.subplot(4,1,1)
      plt.plot(time,ecc,linestyle = 'solid',marker='None',color='darkorchid',linewidth =2,rasterized=True)
      plt.ylabel('Eccentricity')
      if xrange:
        plt.xlim(xrange)

      #e sin(obl) sin varpi
      plt.subplot(4,1,1)
      plt.plot(time,esinv,linestyle = 'solid',marker='None',color='salmon',linewidth=2,rasterized=True)
      plt.ylabel('COPP')
      plt.xlabel('Time [years]')
      if xrange:
        plt.xlim(xrange)

    if dir[ii] == '.':
      dir[ii] = 'cwd'
  
  #fig.suptitle('\n'.join(titlestr),fontsize=20) 
  
  if xrange:
    sfile = 'peak_'+'_'.join(dir)+'_%d_%d.pdf'%(xrange[0],xrange[1])
  else:
    sfile = 'peak_'+'_'.join(dir)+'.pdf'
  plt.savefig(sfile,dpi =300)

  if show:
    plt.show()
  else:
    plt.close()
