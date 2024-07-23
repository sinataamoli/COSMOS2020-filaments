import numpy as np
import scipy as sp
import pandas as pd
import math
import glob

from scipy.stats import norm
from scipy.stats import binned_statistic
from scipy import special
from scipy.integrate import quad
from scipy import integrate
from scipy import stats

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import matplotlib as mpl
from matplotlib import ticker, cm
from matplotlib.colors import LinearSegmentedColormap
# mpl.matplotlib_fname()

from matplotlib import rc
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rc('font',**{'family':'sans-serif','sans-serif':'Helvetica'})
mpl.rc('text', usetex = True)
# mpl.rcParams['font.family'] = ['sans-serif','sans-serif':'Helvetica']
mpl.rc('xtick', labelsize=18)
mpl.rc('ytick', labelsize=18)
mpl.rc('axes', labelsize=22)

title = {'color':  'black', 'weight': 'normal',
         'size': 16}
label = {'color':  'black', 'weight': 'normal',
         'size': 23}

# import matplotlib.animation as animation
# # matplotlib.rcParams.update(matplotlib.rcParamsDefault)
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman"]

from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import Planck15
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astroML.stats import binned_statistic_2d

# from tqdm import tqdm_notebook
from tqdm.auto import tqdm

import os
import time

# import seaborn as sns
import concurrent.futures