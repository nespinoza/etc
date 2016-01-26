from math import modf
import ephem as E
import os
import sys
import numpy as np
import jdcal

"""
       Interactive Transit Ephemerides Calculator
     Author: Nestor Espinoza (nespino@astro.puc.cl)

 Usage: 

   From console, do:  python intercative_calc.py. 
   Insert the year, month and day in which you want a transit 
   and that's it!

 Input data:

   Data from exoplanets.org, with the following format:

   NAME,PER (days),T14 (days),DEPTH,TT,RA (hr),DEC (deg),V,VURL,VREF,TRANSIT,TRANSITURL,TRANSITREF

   This can be easily obtained by asking for Name, Orbital Period, Duration of Transit, Transit Depth, 
   Epoch of Transit Center, RA, DEC, V mag and Transit in the Tables of the webpage.

 Parameters:

   Note that you can change the options below, and choose a maximum
   airmass and magnitude, as well as the observatory data you are 
   proposing time to (current data is from Las Campanas Observatory).
"""
################## Constants #############################
G = 6.67384e-11      # SI
Mjup = 1.89813e27    # kg
Rjup = 69911e3       # m
Rsun = 6.955e8       # m
AU = 149597871e3     # m
k_B = 1.3806488e-23  # J/K
amu = 1.66053892e-27 # kg

################## OPTIONS ###############################
max_airmass = 3.0
max_mag = 20.0
obs = { 'longitude' : '70:42.1', # Data for LCO
        'latitude'  : '-29:0.2',
        'elevation' : 2282,
        'temperature': 7.0,
        'pressure': 760.0,
      }
##########################################################


def getAirmass(ra,dec,day,observatory):
    star = E.FixedBody()
    star._ra = ra
    star._dec = dec
 
    observer = E.Observer()
    observer.date = day
    observer.long = -1*E.degrees(observatory['longitude'])
    observer.lat = E.degrees(observatory['latitude'])
    observer.elevation = observatory['elevation']
    observer.temp = observatory['temperature']
    observer.pressure = observatory['pressure']
    observer.epoch = 2000.0
    observer.horizon = -1*np.sqrt(2*observer.elevation/E.earth_radius)

    star.compute(observer)
    airmass = secz(star.alt)
    return airmass

def secz(alt):
        """Compute airmass"""
        if alt < E.degrees('03:00:00'):
            alt = E.degrees('03:00:00')
        sz = 1.0/np.sin(alt) - 1.0
        xp = 1.0 + sz*(0.9981833 - sz*(0.002875 + 0.0008083*sz))
        return xp

def getPrintDate(year,month,day,hh,mm,ss):
          if(month<10):
             s_month = '0'+str(month)
          else:
             s_month = str(month)
          if(day<10):
             s_day = '0'+str(day)
          else:
             s_day = str(day)
          return str(year)+'-'+s_month+'-'+s_day+' '+hh+':'+mm+':'+ss

def read_data(fname):
        name = []
        period = []    # days
        tduration = [] # days
        tdepth = []    # fractional flux 
        tcenter = []   # days
        ra = []        # hours
        dec = []       # deg
        Vmag = []      # vmag
        a = []         # semi-major axis (au)
        teff = []      # teff (K)
        rstar = []     # star radius (rsun)
        rplanet = []   # planet radius (rjup)
        inc = []       # inclination of orbit (deg)
        msini = []     # mass*sin(i) (mjup)
        f = open(fname)
        reading = True
        counter = 0
        while(reading):
             line = f.readline()
             if(line==''):
                break
             if(line[0]=='#' and counter == 0):
                counter = 1
                all_names = line.split(',')
                for i in range(len(all_names)):
                    if('NAME' in all_names[i]):
                        idx_name = i
                    elif('PER' in all_names[i]):
                        idx_period = i
                    elif('T14' in all_names[i]):
                        idx_duration = i
                    elif('DEPTH' in all_names[i]):
                        idx_depth = i
                    elif('TT' in all_names[i]):
                        idx_center = i
                    elif('RA' in all_names[i]):
                        idx_ra = i
                    elif('DEC' in all_names[i]):
                        idx_dec = i
                    elif(all_names[i]=='V'):
                        idx_vmag = i
                    elif(all_names[i]=='A'):
                        idx_a = i
                    elif('TEFF' in all_names[i]):
                        idx_teff = i
                    elif('RSTAR' in all_names[i]):
                        idx_rstar = i
                    elif(all_names[i]=='I'):
                        idx_inc = i
                    elif('MSINI' in all_names[i]):
                        idx_msini = i
                    elif('R' == all_names[i]):
                        idx_rplanet = i

                # Search every name for each data value 
             elif(line[0]!='#'):
                splitted = line.split(',')
                planet_ok = True
                # If basic info used in this code is missing from the planet...
                if(splitted[idx_a]=='' or splitted[idx_teff]=='' or splitted[idx_rstar]=='' or splitted[idx_inc]==''\
                   or splitted[idx_msini]=='' or splitted[idx_vmag]=='' or splitted[idx_duration]==''\
                   or splitted[idx_center]=='' or splitted[idx_depth]==''):
                     planet_ok = False
                if(planet_ok): # It transits and all info is known:
                   name.append(splitted[idx_name])
                   period.append(np.float(splitted[idx_period]))
                   tduration.append(np.float(splitted[idx_duration]))
                   tdepth.append(np.float(splitted[idx_depth]))
                   tcenter.append(np.float(splitted[idx_center]))
                   ra.append(np.float(splitted[idx_ra]))
                   dec.append(np.float(splitted[idx_dec]))
                   Vmag.append(np.float(splitted[idx_vmag]))
                   a.append(AU*np.float(splitted[idx_a]))
                   teff.append(np.float(splitted[idx_teff]))
                   rstar.append(Rsun*np.float(splitted[idx_rstar]))
                   rplanet.append(Rjup*np.float(splitted[idx_rplanet]))
                   inc.append(np.float(splitted[idx_inc]))
                   msini.append(Mjup*np.float(splitted[idx_msini]))
        return name,period,tduration,tdepth,tcenter,ra,dec,Vmag,a,teff,rstar,rplanet,inc,msini

def deg2deg(deg):
    degrees = int(modf(deg)[1])
    minutes = np.abs(modf(deg)[0]*60.0)
    seconds = np.abs(modf(minutes)[0]*60)
    mm = int(modf(minutes)[1])
    ss = np.round(seconds,2) #int(modf(seconds)[1])
    if(np.abs(degrees)<10):
       if(degrees>0):
          s_deg = '0'+str(degrees)
       else:
          s_deg = '-0'+str(-degrees)
    else:
       s_deg = str(degrees)
    if(mm<10):
       s_min = '0'+str(mm)
    else:
       s_min = str(mm)
    if(ss<10):
       s_sec = '0'+str(ss)
    else:
       s_sec = str(ss)
    return s_deg+':'+s_min+':'+s_sec

def getCalDay(JD):
	year, month, day, hour= jdcal.jd2gcal(JD,0.0)
	hour = hour*24
	minutes = modf(hour)[0]*60.0
	seconds = modf(minutes)[0]*60.0
	hh = int(modf(hour)[1])
	mm = int(modf(minutes)[1])
	ss = seconds
	if(hh<10):
	   hh = '0'+str(hh)
	else:
	   hh = str(hh)
	if(mm<10):
	   mm = '0'+str(mm)
	else:
	   mm = str(mm)
	if(ss<10):
	   ss = '0'+str(np.round(ss,1))
	else:
	   ss = str(np.round(ss,1))
        return year,month,day,hh,mm,ss

def getTime(year,month,day,hh,mm,ss):
    return str(year)+'/'+str(month)+'/'+str(day)+' '+hh+':'+mm+':'+ss

chosen_planet = raw_input('Which planet?')
short_pname = "".join((chosen_planet.lower()).split())
name,period,tduration,tdepth,tcenter,ra,dec,Vmag,a,teff,rstar,rplanet,inc,msini = read_data('data.dat')
for i in range(len(name)):
        n = name[i]
        if "".join((n.lower()).split()) == short_pname:
                name,period,tduration,tdepth,tcenter,ra,dec,Vmag,a,teff,rstar,rplanet,inc,msini = \
                [name[i]],[period[i]],[tduration[i]],[tdepth[i]],[tcenter[i]],[ra[i]],[dec[i]],[Vmag[i]],\
                [a[i]],[teff[i]],[rstar[i]],[rplanet[i]],[inc[i]],[msini[i]]
                break
if len(name)>1:
        print 'Planet '+chosen_planet+' not found'
##########################################################
u_year = int(raw_input('Insert year on which you want to search for transits: '))
u_month = int(raw_input('Month number? (1: January, 12:December): '))
#u_day = int(raw_input('Day?: '))
##########################################################
names = np.array([])
tdepths = np.array([])
ras = np.array([])
decs = np.array([])
airmasses = np.array([])
airmasses_e = np.array([])
airmasses_mt = np.array([])
i_dates = np.array([])
m_dates = np.array([])
e_dates = np.array([])
vmags = np.array([])
t = 0
all_days = [31,28,31,30,31,30,31,31,30,31,30,31]
for u_day in range(1,all_days[u_month-1]+1):
   year,month,day,hh,mm,ss = getCalDay(tcenter[t])
   imax = int((u_year+1-year)*366.0/period[t])
   for i in range(imax):
       # Obtain date data for ingress:
       year,month,day,hh,mm,ss = getCalDay(tcenter[t]+i*period[t]-(tduration[t]/2.0))
       # Obtain date data for mid-transit:
       year_mt,month_mt,day_mt,hh_mt,mm_mt,ss_mt = getCalDay(tcenter[t]+i*period[t])
       # Obtain date data for egress:
       year_e,month_e,day_e,hh_e,mm_e,ss_e = getCalDay(tcenter[t]+i*period[t]+(tduration[t]/2.0))  
       time = getTime(year,month,day,hh,mm,ss)
       time_e = getTime(year_e,month_e,day_e,hh_e,mm_e,ss_e)
       time_mt = getTime(year_mt,month_mt,day_mt,hh_mt,mm_mt,ss_mt)
       airmass = getAirmass(deg2deg(ra[t]),deg2deg(dec[t]),time,obs)
       airmass_e = getAirmass(deg2deg(ra[t]),deg2deg(dec[t]),time_e,obs)
       airmass_mt = getAirmass(deg2deg(ra[t]),deg2deg(dec[t]),time_mt,obs)
       if(( (year==u_year and month==u_month and day==u_day and np.float(hh)>=0) or (year==u_year and month==u_month and day==u_day-1 and np.float(hh)>=23))):
          names = np.append(names,name[t])
          tdepths = np.append(tdepths,tdepth[t])
          ras = np.append(ras,deg2deg(ra[t]))
          decs = np.append(decs,deg2deg(dec[t]))
          airmasses = np.append(airmasses,airmass)
          airmasses_e = np.append(airmasses_e,airmass_e)
          airmasses_mt = np.append(airmasses_mt,airmass_mt)
          vmags = np.append(vmags,Vmag[t])
          i_dates = np.append(i_dates,getPrintDate(year,month,day,hh,mm,ss))
          m_dates = np.append(m_dates,getPrintDate(year_mt,month_mt,day_mt,hh_mt,mm_mt,ss_mt))
          e_dates = np.append(e_dates,getPrintDate(year_e,month_e,day_e,hh_e,mm_e,ss_e))

idx = tdepths.argsort()
idx = idx[::-1]
names = names[idx]
ras = ras[idx]
decs = decs[idx]
i_dates = i_dates[idx]
m_dates = m_dates[idx]
e_dates = e_dates[idx]
tdepths = tdepths[idx]
airmasses = airmasses[idx]
airmasses_e = airmasses_e[idx]
airmasses_mt = airmasses_mt[idx]
vmags = vmags[idx]
for i in range(len(tdepths)):
    print 'Planet name         '+names[i]
    print 'Vmag              : '+str(vmags[i])
    print 'Coordinates       : '+ras[i]+' '+decs[i]
    print 'Transit starts at : '+i_dates[i]+' (Airmass:'+str(np.round(airmasses[i],2))+')'
    print 'Mid-transit at    : '+m_dates[i]+' (Airmass:'+str(np.round(airmasses_mt[i],2))+')'
    print 'Transit ends at   : '+e_dates[i]+' (Airmass:'+str(np.round(airmasses_e[i],2))+')'
    print 'Transit depth     : '+str(tdepths[i])
    raw_input('Next? (press enter, ctrl+z to exit)')
    os.system('clear')
