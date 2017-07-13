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
    
    titlestr = []
    titlestr.append(r'$a = %f, pCO_2 = %f$'%(semi,pco2))
    titlestr.append(r'$e_0 = %f, i_0 = %f^{\circ}, \psi_0 = %f^{\circ}, P_{rot} = %f$ d'%(ecc[0],inc[0],obl[0],P))
    fig.subplots_adjust(wspace=0.3)

    lats = np.unique(body.Latitude)
    nlats = len(lats)
    ntimes = len(body.Time)
    
    # plot temperature
    temp = np.reshape(body.TempLat,(ntimes,nlats))
    if orbit == True:
      ax1 = plt.subplot(4,2,1)
    else:
      ax1 = plt.subplot(5,nfiles,ii+1)
    pos = ax1.figbox.get_points()
    c = plt.contourf(body.Time,lats,temp.T,cmap='plasma')
    plt.ylabel('Latitude')
    plt.title(r'Surface Temp [$^{\circ}$C]')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange == False:
      left = 0
    else:
      left = xrange[0]
    plt.text(left,140,'\n'.join(titlestr),fontsize=20) 
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    # plot albedo
    alb = np.reshape(body.AlbedoLat,(ntimes,nlats))
    if orbit == True:
      ax2 = plt.subplot(4,2,3)
    else:
      ax2 = plt.subplot(5,nfiles,ii+2*nfiles+1)
    pos = ax2.figbox.get_points()
    c = plt.contourf(body.Time,lats,alb.T,cmap='Blues_r')
    plt.ylabel('Latitude')
    plt.title('Albedo (TOA)')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  

    # plot ice height
    ice = np.reshape(body.IceHeight,(ntimes,nlats))
    if orbit == True:
      ax3 = plt.subplot(4,2,5)
    else:
      ax3 = plt.subplot(5,nfiles,ii+3*nfiles+1)
    pos = ax3.figbox.get_points()
    c = plt.contourf(body.Time,lats,ice.T,cmap='Blues_r')
    plt.ylabel('Latitude')
    plt.title('Ice sheet height [m]')
    plt.ylim(-90,90)
  #   plt.xlim(0,2e6)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
    # ax3p = ax3.twinx()
  #   plt.plot(body.Time,esinv,linestyle = 'solid',marker='None',color='salmon',linewidth=2)
  

    # plot bedrock
    brock = np.reshape(body.BedrockH,(ntimes,nlats))
    if orbit == True:
      ax4 = plt.subplot(4,2,7)
    else:
      ax4 = plt.subplot(5,nfiles,ii+4*nfiles+1)
    pos = ax4.figbox.get_points()
    c = plt.contourf(body.Time,lats,brock.T,cmap='Reds_r')
    plt.ylabel('Latitude')
    plt.title('Bedrock height [m]')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    plt.xlabel('Time [years]')
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  

    # plot insolation
    insol = np.reshape(body.AnnInsol,(ntimes,nlats))
    if orbit == True:
      ax5 = plt.subplot(4,2,2)
    else:
      ax5 = plt.subplot(5,nfiles,ii+nfiles+1)
    pos = ax5.figbox.get_points()
    c = plt.contourf(body.Time,lats,insol.T,cmap='plasma')
    plt.ylabel('Latitude')
    plt.title(r'Annual average insolation [w/m$^2$]')
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
  
  #fig.suptitle('\n'.join(titlestr),fontsize=20) 
  
  if xrange:
    sfile = 'evol_'+'_'.join(dir)+'_%d_%d.pdf'%(xrange[0],xrange[1])
  else:
    sfile = 'evol_'+'_'.join(dir)+'.pdf'
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
    maxscl = np.max(body.TempMaxLat)
    minscl = np.min(body.TempMinLat)
    norm = mplcol.Normalize(vmin=minscl,vmax=maxscl)
  
  
    temp = np.reshape(body.TempLat,(ntimes,nlats))
    if orbit == True:
      ax1 = plt.subplot(3,2,1)
    else:
      ax1 = plt.subplot(3,nfiles,ii+1)
    pos = ax1.figbox.get_points()
    c = plt.contourf(body.Time,lats,temp.T,cmap='plasma',norm = norm)
    plt.ylabel('Latitude')
    plt.title('Surface Temp [$^{\circ}$C]')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    tempmin = np.reshape(body.TempMinLat,(ntimes,nlats))
    if orbit == True:
      ax3 = plt.subplot(3,2,3)
    else:
      ax3 = plt.subplot(3,nfiles,ii+nfiles+1)
    pos = ax3.figbox.get_points()
    c = plt.contourf(body.Time,lats,tempmin.T,cmap='plasma',norm = norm)
    plt.ylabel('Latitude')
    plt.title('Min annual surface Temp [$^{\circ}$C]')
    plt.ylim(-90,90)
    plt.yticks([-60,-30,0,30,60])
    if xrange:
      plt.xlim(xrange)
    plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
    tempmax = np.reshape(body.TempMaxLat,(ntimes,nlats))
    if orbit == True:
      ax5 = plt.subplot(3,2,5)
    else:
      ax5 = plt.subplot(3,nfiles,ii+2*nfiles+1)
    pos = ax5.figbox.get_points()
    c = plt.contourf(body.Time,lats,tempmax.T,cmap='plasma',norm = norm)
    plt.ylabel('Latitude')
    plt.title('Max annual surface Temp [$^{\circ}$C]')
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
