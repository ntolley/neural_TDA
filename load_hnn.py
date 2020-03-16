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

debug = 1

def downsample (olddata,oldrate,newrate):  
  ratio=oldrate/float(newrate) # Calculate ratio of sampling rates
  newdata = scipy.signal.decimate(olddata, int(ratio), ftype='fir',zero_phase=True)
  return newdata  

def index2ms (idx, sampr): return 1e3*idx/sampr
def ms2index (ms, sampr): return int(sampr*ms/1e3)

# bandpass filter the items in lfps. lfps is a list or numpy array of LFPs arranged spatially by column
def getbandpass (lfps,sampr,minf=0.05,maxf=300):
  datband = []
  for i in range(len(lfps[0])): datband.append(bandpass(lfps[:,i],minf,maxf,df=sampr,zerophase=True))
  datband = numpy.array(datband)
  return datband

# lowpass filter the items in lfps. lfps is a list or numpy array of LFPs arranged spatially by row
def getlowpass (lfps,sampr,maxf):
  datband = []
  for i in range(len(lfps[0])): datband.append(lowpass(lfps[:,i],maxf,df=sampr,zerophase=True))
  datband = numpy.array(datband)
  return datband

# gets 2nd spatial derivative of voltage as approximation of CSD.
# performs lowpass filter on voltages before taking spatial derivative
# input dlfp is dictionary of LFP voltage time-series keyed by (trial, electrode)
# output dCSD is keyed by trial
def getCSD (lfps,sampr,maxlfp,ntrial,spacing_um,minf=0.1,maxf=300.0):
  spacing_mm = spacing_um/1000

  CSD = {}
  for trial in range(ntrial): # go through trials
    temp = lfps[trial].T

    # convert from uV to mV
    temp = temp/1000

    datband = getlowpass(temp,sampr,maxf)
    if datband.shape[0] > datband.shape[1]: # take CSD along smaller dimension
      ax = 1
    else:
      ax = 0


    # when drawing CSD make sure that negative values (depolarizing intracellular current) drawn in red,
    # and positive values (hyperpolarizing intracellular current) drawn in blue
    CSD[trial] = -numpy.diff(datband,n=2,axis=ax)/spacing_mm**2 # now each column (or row) is an electrode -- CSD along electrodes

    CSD[trial] = Vaknin(CSD[trial])

  return CSD


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

# get MUA - first do a bandpass filter then rectify. 
#  lfps is a list or numpy array of LFPs arranged spatially by column
def getMUA (lfps,sampr,minf=300,maxf=5000):
  datband = getbandpass(lfps,sampr,minf,maxf)
  return numpy.abs(datband)

#

def readLFPs (basedir):
  lfps = {}
  LFP = {}
  lfile = os.listdir(basedir)
  maxlfp = 0; ntrial_indx = 0; tvec = None
  for f in lfile:
    if f.count('lfp_') > 0 and f.endswith('.txt'):
      lf = f.split('.txt')[0].split('_')
      try:
        trial = int(lf[1])
        nlfp = int(lf[2])
      except IndexError:
        trial = 0
        nlfp = int(lf[1])
      maxlfp = max(nlfp,maxlfp)
      ntrial_indx = max(trial,ntrial_indx)
      fullpath = os.path.join(basedir,f)
      try:
        temp = np.loadtxt(fullpath)
        lfps[(trial,nlfp)]=np.array(temp[:,1])
        if tvec is None: # only need to do this once
          tvec = temp[:,0]
      except:
        print('exception!')

  for trial in range(ntrial_indx+1): # go through trials
    LFP[trial]=np.zeros((maxlfp+1,len(tvec)))
    for nlfp in range(maxlfp+1):
      LFP[trial][nlfp,:] = lfps[(trial,nlfp)]

  sampr = 1e3 / (tvec[1]-tvec[0])
  dt = 1.0 / sampr # time-step in seconds
  return sampr, LFP, dt, tvec, maxlfp, ntrial_indx+1

def loadHNNdir (basedir,spacing_um):
  sampr, LFP, dt, tt, maxlfp, ntrial = readLFPs(basedir)
  CSD = getCSD(LFP,sampr,maxlfp,ntrial,spacing_um)
  return sampr,LFP,dt,tt,CSD,maxlfp,ntrial

  # get the average ERP (dat should be either LFP or CSD)
def getAvgERP (dat, sampr, tvec, maxlfp, ntrial):
  avgERP = np.zeros((maxlfp+1,len(tvec)))
  for chan in range(maxlfp+1): # go through channels
    for trial in range(ntrial): # go through trials
      avgERP[chan,:] += dat[trial][chan,:]
    avgERP[chan,:] /= ntrial
  return avgERP
