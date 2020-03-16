import matplotlib
matplotlib.use("TkAgg")
from pylab import *
from scipy.interpolate import interp1d
import sys,os,numpy,scipy,subprocess
from math import ceil
import tables # for reading matlab file format (HDF5)
from filter import lowpass,bandpass
from scipy.signal import decimate
import pickle
import h5py
import matplotlib.patches as mpatches
import sys

#rcParams['lines.markersize'] = 15
#rcParams['lines.linewidth'] = 1
#rcParams['font.size'] = 25

#
def downsample (olddata,oldrate,newrate):  
  ratio=oldrate/float(newrate) # Calculate ratio of sampling rates
  newdata = scipy.signal.decimate(olddata, int(ratio), ftype='fir',zero_phase=True)
  return newdata  

def index2ms (idx, sampr): return 1e3*idx/sampr
def ms2index (ms, sampr): return int(sampr*ms/1e3)

# read the matlab .mat file and return the sampling rate and electrophys data
def rdmat (fn,samprds=0):  
  fp = h5py.File(fn,'r') # open the .mat / HDF5 formatted data
  sampr = fp['craw']['adrate'][0][0] # sampling rate
  dt = 1.0 / sampr # time-step in seconds
  dat = fp['craw']['cnt'] # cnt record stores the electrophys data
  npdat = zeros(dat.shape)
  tmax = ( len(npdat) - 1.0 ) * dt # use original sampling rate for tmax - otherwise shifts phase
  dat.read_direct(npdat) # read it into memory
  fp.close()
  if samprds > 0.0: # resample the LFPs
    dsfctr = sampr/samprds
    dt = 1.0 / samprds
    siglen = max((npdat.shape[0],npdat.shape[1]))
    nchan = min((npdat.shape[0],npdat.shape[1]))
    npds = [] # zeros((int(siglen/float(dsfctr)),nchan))
    # print dsfctr, dt, siglen, nchan, samprds, ceil(int(siglen / float(dsfctr))), npds.shape
    for i in range(nchan): 
      print('resampling channel', i)
      npds.append(downsample(npdat[:,i], sampr, samprds))
    npdat = np.array(npds)
    npdat = npdat.T
    sampr = samprds
  tt = numpy.linspace(0,tmax,len(npdat)) # time in seconds
  return sampr,npdat,dt,tt # sampling rate, LFP data, time-step, time (in seconds)

#
def getHDF5values (fn,key):
  fp = h5py.File(fn,'r')
  hdf5obj = fp[key]
  x = np.array(fp[hdf5obj.name])
  val = [y[0] for y in fp[x[0,0]].value]
  fp.close()
  return val

#
def getTriggerTimes (fn):
  fp = h5py.File(fn,'r')
  hdf5obj = fp['trig/anatrig']
  x = np.array(fp[hdf5obj.name])
  val = [y[0] for y in fp[x[0,0]].value]
  fp.close()
  return val  

#
def getTriggerIDs (fn):
  fp = h5py.File(fn,'r')
  hdf5obj = fp['trig/ttype']
  x = np.array(fp[hdf5obj.name])
  val = fp[x[0,0]].value[0] 
  val = [int(x) for x in val]
  fp.close()
  return val

# bandpass filter the items in lfps. lfps is a list or numpy array of LFPs arranged spatially by column
def getbandpass (lfps,sampr,minf=0.05,maxf=300):
  datband = []
  for i in range(len(lfps[0])): datband.append(bandpass(lfps[:,i],minf,maxf,df=sampr,zerophase=True))
  datband = numpy.array(datband)
  return datband

# Vaknin correction for CSD analysis
# Allows CSD to be performed on all N contacts instead of N-2 contacts
# See Vaknin et al (1989) for more details
def Vaknin(x):
    # Preallocate array with 2 more rows than input array
    x_new = np.zeros((x.shape[0]+2, x.shape[1]))

    # Duplicate first and last row of x into first and last row of x_new
    x_new[0, :] = x[0, :]
    x_new[-1, :] = x[-1, :]

    # Duplicate all of x into middle rows of x_new
    x_new[1:-1, :] = x

    return x_new

# get CSD - first do a lowpass filter. lfps is a list or numpy array of LFPs arranged spatially by column
def getCSD (lfps,sampr,spacing_um,minf=0.05,maxf=300):

  # convert from uV to mV
  lfps = lfps/1000

  datband = getbandpass(lfps,sampr,minf,maxf)
  if datband.shape[0] > datband.shape[1]: # take CSD along smaller dimension
    ax = 1
  else:
    ax = 0

  spacing_mm = spacing_um/1000
  # when drawing CSD make sure that negative values (depolarizing intracellular current) drawn in red,
  # and positive values (hyperpolarizing intracellular current) drawn in blue
  CSD = -numpy.diff(datband,n=2,axis=ax)/spacing_mm**2 # now each column (or row) is an electrode -- CSD along electrodes

  CSD = Vaknin(CSD)

  return CSD

# get MUA - first do a bandpass filter then rectify. 
#  lfps is a list or numpy array of LFPs arranged spatially by column
def getMUA (lfps,sampr,minf=300,maxf=5000):
  datband = getbandpass(lfps,sampr,minf,maxf)
  return numpy.abs(datband)

#
def loadfile (fn,samprds,spacing_um=100):
  # load a .mat data file (fn) using sampling rate samprds (should be integer factor of original sampling rate 44000),
  # returns: sampr is sampling rate after downsampling
  #          LFP is laminar local field potential data
  #          dt is time-step (redundant with sampr)
  #          tt is time array (in seconds)
  #          CSD is laminar current source density signal
  #          trigtimes is array of stimulus trigger indices (indices into arrays)
  sampr,LFP,dt,tt=rdmat(fn,samprds=samprds) #
  sampr,dt,tt[0],tt[-1] # (2000.0, 0.001, 0.0, 1789.1610000000001)
  CSD = getCSD(LFP,sampr,spacing_um)
  divby = 44e3 / samprds
  trigtimes = None
  try: # not all files have stimuli
    trigtimes = [int(round(x)) for x in np.array(getTriggerTimes(fn)) / divby] # divby since downsampled signals by factor of divby
  except:
    pass
  #trigIDs = getTriggerIDs(fn)
  LFP = LFP.T # make sure each row is a channel
  return sampr,LFP,dt,tt,CSD,trigtimes

samprds = 11000.0 # downsampling to this frequency
divby = 44e3 / samprds

#first
#fn = '1-rb067068027@os.mat'
#second
#fn = 'data/1-rb067068029@os.mat'
#third
#fn = 'data/1-rb067068030@os.mat'
#fourth
#fn = 'data/1-rb067068031@os.mat'

spacing_um = 100
#sampr,LFP,dt,tt,CSD,trigtimes = loadfile(fn,samprds, spacing_um)

def calPosThresh(dat, sigmathresh):
  #return dat.mean() + sigmathresh * dat.std()
  return 100

def calNegThresh(dat, sigmathresh):
  #return dat.mean() - sigmathresh * dat.std()
  return -100

# remove noise, where noise < negthres < dat < posthres < noise
def badEpoch (dat, sigmathresh):
  badValues = len(np.where(dat <= calNegThresh(dat, sigmathresh))[0]) + \
              len(np.where(dat >= calPosThresh(dat, sigmathresh))[0])
  if badValues > 0:
    return True
  else:
    return False

def removeBadEpochs (dat, sampr, trigtimes, swindowms, ewindowms, sigmathresh):
  nrow = dat.shape[0]
  swindowidx = ms2index(swindowms,sampr) # could be negative
  ewindowidx = ms2index(ewindowms,sampr)

  # trigByChannel could be returned for removing different epochs on each channel
  trigByChannel = [x for x in range(nrow)]
  badEpochs = []
  for chan in range(nrow): # go through channels
    trigByChannel[chan] = []
    for trigidx in trigtimes: # go through stimuli
      sidx = max(0,trigidx+swindowidx)
      eidx = min(dat.shape[1],trigidx+ewindowidx)
      if not badEpoch(dat[chan, sidx:eidx], sigmathresh):
        trigByChannel[chan].append(trigidx)
      else:
        badEpochs.append(trigidx)
    print('Found %d bad epochs in channel %d. Range: [%.2f, %.2f]'%
          (len(trigtimes) - len(trigByChannel[chan]), chan,
           calNegThresh(dat[chan, sidx:eidx], sigmathresh),
           calPosThresh(dat[chan, sidx:eidx], sigmathresh)))

  # combine bad epochs into a single sorted list (without duplicates)
  badEpochs = sort(list(set(badEpochs)))
  print('%d bad epochs:'%len(badEpochs),[x for x in badEpochs])

  # remove the associated trigger times before returning
  trigtimes = np.delete(trigtimes,[trigtimes.index(x) for x in badEpochs])

  return trigtimes

# get the average ERP (dat should be either LFP or CSD)
def getAvgERP (dat, sampr, trigtimes, swindowms, ewindowms):
  nrow = dat.shape[0]
  tt = np.linspace(swindowms, ewindowms,ms2index(ewindowms - swindowms,sampr))
  swindowidx = ms2index(swindowms,sampr) # could be negative
  ewindowidx = ms2index(ewindowms,sampr)
  avgERP = np.zeros((nrow,len(tt)))
  for chan in range(nrow): # go through channels
    for trigidx in trigtimes: # go through stimuli
      sidx = max(0,trigidx+swindowidx)
      eidx = min(dat.shape[1],trigidx+ewindowidx)
      avgERP[chan,:] += dat[chan, sidx:eidx]
    avgERP[chan,:] /= float(len(trigtimes))
  return tt,avgERP 


  # get the average ERP (dat should be either LFP or CSD)
def getERPstats (dat, sampr, trigtimes, swindowms, ewindowms):
  nrow = dat.shape[0]
  tt = np.linspace(swindowms, ewindowms,load.ms2index(ewindowms - swindowms,sampr))
  swindowidx = load.ms2index(swindowms,sampr) # could be negative
  ewindowidx = load.ms2index(ewindowms,sampr)
  avgERP = np.zeros((nrow,len(tt)))
  semERP = np.zeros((nrow,len(tt)))  

  for chan in range(nrow): # go through channels
    for trigidx in trigtimes: # go through stimuli
      sidx = max(0,trigidx+swindowidx)
      eidx = min(dat.shape[1],trigidx+ewindowidx)
      avgERP[chan,:] += dat[chan, sidx:eidx]

    #Divide to get avg      
    avgERP[chan,:] /= float(len(trigtimes))
     
       
    for trigidx in trigtimes:
      sidx = max(0,trigidx+swindowidx)
      eidx = min(dat.shape[1],trigidx+ewindowidx)
      semERP[chan,:] += (dat[chan, sidx:eidx] - avgERP[chan,:]) ** 2 #Sum of squares

    #Divide by sample size
    semERP[chan,:] /= float(len(trigtimes))
    semERP[chan,:] /= np.sqrt(float(len(trigtimes)))



  return tt,avgERP, semERP

# draw the average ERP (dat should be either LFP or CSD)
def drawAvgERP (dat, sampr, trigtimes, windowms, whichchan=None, yl=None, clr=None,lw=1):
  ttavg,avgERP = getAvgERP(dat,sampr,trigtimes,windowms)
  nrow = avgERP.shape[0]
  for chan in range(nrow): # go through channels
    if whichchan is None:
      subplot(nrow,1,chan+1)
      plot(ttavg,avgERP[chan,:],color=clr,linewidth=lw)
    elif chan==whichchan:
      plot(ttavg,avgERP[chan,:],color=clr,linewidth=lw)
    xlim((-windowms,windowms))
    if yl is not None: ylim(yl)
  
# draw the event related potential (or associated CSD signal), centered around stimulus start (aligned to t=0)
def drawERP (dat, sampr, trigtimes, windowms, whichchan=None, yl=None,clr=None,lw=1):
  if clr is None: clr = 'gray'
  nrow = dat.shape[0]
  tt = np.linspace(-windowms,windowms,ms2index(windowms*2,sampr))
  windowidx = ms2index(windowms,sampr)
  for trigidx in trigtimes: # go through stimuli
    for chan in range(nrow): # go through channels
      sidx = max(0,trigidx-windowidx)
      eidx = min(dat.shape[1],trigidx+windowidx)
      if whichchan is None:
        subplot(nrow,1,chan+1)
        plot(tt,dat[chan, sidx:eidx],color=clr,linewidth=lw)
      elif chan==whichchan:
        plot(tt,dat[chan, sidx:eidx],color=clr,linewidth=lw)
      xlim((-windowms,windowms))
      if yl is not None: ylim(yl)
      #xlabel('Time (ms)')

