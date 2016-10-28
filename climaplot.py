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

vplred = '#C91111'
vplorg = '#E09401'
vpllbl = '#13AED5'
vpldbl = '#1321D8'
vplpur = '#642197'

def seasonal_maps(time, dir = '.'):
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
        check = 1
        
    if check == 0:
      raise StandardError('Climate data not found for time %f'%time)
    
    insol = np.loadtxt(insolf,unpack=True)
    temp = np.loadtxt(tempf,unpack=True)
    ice = np.loadtxt(icef,unpack=True)
    output = vplot.GetOutput(dir)

    ctmp = 0
    for p in range(len(output.bodies)):
      if output.bodies[p].name == plname:
        body = output.bodies[p]
        ctmp = 1
      else:
        if p == len(output.bodies)-1 and ctmp == 0:
          raise Exception("Planet %s not found in folder %s"%(plname,folders[j]))
        
    lats = body.Latitude[0]
    try:
      obl = planet.Obliquity[np.where(body.Time==time)[0]]
    except:
      obl = getattr(output.log.initial,plname).Obliquity
      if obl.unit == 'rad':
        obl *= 180/np.pi

    try:
      ecc = planet.Eccentricity[np.where(planet.Time==time)[0]]
    except:
      ecc = getattr(output.log.initial,plname).Eccentricity
      
    try:
      longp = (planet.LongP+planet.PrecA)[np.where(planet.Time==time)[0]]
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

    fig = plt.figure(figsize=(16,6))
    fig.suptitle('Time = %f, Obl = %f, Ecc = %f, LongP = %f'%(time,obl,ecc,longp),fontsize=20)
    plt.subplot(1,3,1)
    plt.title('insolation [W/m^2] (1 orbit)',fontsize=12)
    c1=plt.contourf(range(np.shape(insol)[1]),lats,insol,cmap='plasma')
    plt.colorbar(c1)
    plt.ylim(lats[0],lats[-1])
    
    plt.subplot(1,3,2)
    c2=plt.contourf(range(np.shape(temp)[1]),lats,temp,cmap='plasma')
    plt.title('surface temp [C] (4 orbits)',fontsize=12)
    plt.colorbar(c2)
    plt.ylim(lats[0],lats[-1])
    
    plt.subplot(1,3,3)
    c3=plt.contourf(range(np.shape(ice)[1]),lats,ice,cmap='Blues_r')
    plt.title('ice balance [kg/m^2/s] (1 orbit)',fontsize=12)
    plt.colorbar(c3)
    plt.ylim(lats[0],lats[-1])
    
    plt.savefig('surf_seas_%.0f.pdf'%time)
    plt.close()
    
def clim_evol(plname,dir='.',xrange=False):
  out = vplot.GetOutput(dir)
  
  ctmp = 0
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
    ecc = np.zeros_like(body.Time)
    
  try:
    inc = body.Inc
  except:
    inc = np.zeros_like(body.Time)
    
  try:
    obl = body.Obliquity
  except:
    obl = np.zeros_like(body.Time)

  f = open(dir+'/'+plname+'.in','r')
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

  fig = plt.figure(figsize=(16,12))
  fig.suptitle('$e_0 = %f, i_0 = %f, \psi_0 = %f, P_{rot} = %f$'%(ecc[0],inc[0],obl[0],P)) 
  fig.subplots_adjust(wspace=0.3)

  lats = np.unique(body.Latitude)
  nlats = len(lats)
  ntimes = len(body.Time)

  # plot temperature
  temp = np.reshape(body.TempLat,(ntimes,nlats))
  ax1 = plt.subplot(4,2,1)
  pos = ax1.figbox.get_points()
  c = plt.contourf(body.Time,lats,temp.T,cmap='plasma')
  plt.ylabel('Latitude')
  plt.title(r'Surface Temp [$^{\circ}$C]')
  plt.ylim(-90,90)
  plt.yticks([-60,-30,0,30,60])
  if xrange:
    plt.xlim(xrange)
  plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  

  # plot albedo
  alb = np.reshape(body.AlbedoLat,(ntimes,nlats))
  ax2 = plt.subplot(4,2,3)
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
  ax3 = plt.subplot(4,2,5)
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
  ax4 = plt.subplot(4,2,7)
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
  ax5 = plt.subplot(4,2,2)
  pos = ax5.figbox.get_points()
  c = plt.contourf(body.Time,lats,insol.T,cmap='plasma')
  plt.ylabel('Latitude')
  plt.title(r'Annual average insolation [w/m$^2$]')
  plt.ylim(-90,90)
  plt.yticks([-60,-30,0,30,60])
  if xrange:
    plt.xlim(xrange)
  plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))


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

  if dir == '.':
    dir = 'cwd'
  
  if xrange:
    sfile = 'evol_'+dir+'_%d_%d.pdf'%(xrange[0],xrange[1])
  else:
    sfile = 'evol_'+dir+'.pdf'
  plt.savefig(sfile)
  plt.close()



def tempminmax(plname,dir='.',xrange=False):
  out = vplot.GetOutput(dir)
  
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
    ecc = np.zeros_like(body.Time)
    
  try:
    inc = body.Inc
  except:
    inc = np.zeros_like(body.Time)
    
  try:
    obl = body.Obliquity
  except:
    obl = np.zeros_like(body.Time)

  f = open(dir+'/'+plname+'.in','r')
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

  fig = plt.figure(figsize=(16,12))
  fig.suptitle('$e_0 = %f, i_0 = %f, \psi_0 = %f, P_{rot} = %f$'%(ecc[0],inc[0],obl[0],P)) 
  fig.subplots_adjust(wspace=0.3)

  lats = np.unique(body.Latitude)
  nlats = len(lats)
  ntimes = len(body.Time)

  # plot temperature
  maxscl = np.max(body.TempMaxLat)
  minscl = np.min(body.TempMinLat)
  norm = mplcol.Normalize(vmin=minscl,vmax=maxscl)
  
  
  temp = np.reshape(body.TempLat,(ntimes,nlats))
  ax1 = plt.subplot(3,2,1)
  pos = ax1.figbox.get_points()
  c = plt.contourf(body.Time,lats,temp.T,cmap='plasma',norm = norm)
  plt.ylabel('Latitude')
  plt.title('Surface Temp [K]')
  plt.ylim(-90,90)
  plt.yticks([-60,-30,0,30,60])
  if xrange:
    plt.xlim(xrange)
  plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
  tempmin = np.reshape(body.TempMinLat,(ntimes,nlats))
  ax3 = plt.subplot(3,2,3)
  pos = ax3.figbox.get_points()
  c = plt.contourf(body.Time,lats,tempmin.T,cmap='plasma',norm = norm)
  plt.ylabel('Latitude')
  plt.title('Min annual surface Temp [K]')
  plt.ylim(-90,90)
  plt.yticks([-60,-30,0,30,60])
  if xrange:
    plt.xlim(xrange)
  plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
  tempmax = np.reshape(body.TempMaxLat,(ntimes,nlats))
  ax5 = plt.subplot(3,2,5)
  pos = ax5.figbox.get_points()
  c = plt.contourf(body.Time,lats,tempmax.T,cmap='plasma',norm = norm)
  plt.ylabel('Latitude')
  plt.title('Max annual surface Temp [K]')
  plt.ylim(-90,90)
  plt.yticks([-60,-30,0,30,60])
  if xrange:
    plt.xlim(xrange)
  plt.colorbar(c,cax=plt.axes([pos[1,0]+0.01,pos[0,1],0.01,pos[1,1]-pos[0,1]]))
  
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

  if xrange:
    sfile = 'temp_'+dir+'_%d_%d.pdf'%(xrange[0],xrange[1])
  else:
    sfile = 'temp_'+dir+'.pdf'
  plt.savefig(sfile)
  plt.close()