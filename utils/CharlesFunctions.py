import matplotlib.pyplot as plt
import numpy as np
import scipy
import os
import re
import math
import warnings
from scipy import optimize
from scipy import interpolate
from scipy import integrate
from scipy import special



def readfile(filename):
    doc = list(open(filename, "r"))
    holder = [[],[]]
    for i in range(len(doc)):
        if(doc[i]):
            txt = doc[i].split()
            holder[0].append(float(txt[0]))
            holder[1].append(float(txt[1]))
    return np.array(holder)

def AvgFuncs(arrs):
    temp = [[],[]]
    for i in range(len(arrs[0][0])):
        holder = 0
        temp[0].append(arrs[0][0][i])
        for arr in arrs:
            holder = holder + arr[1][i]
        temp[1].append(holder)
    return np.array(temp)


def find_closest_index(array, value):
    closest_index = None
    min_difference = float('inf')

    if len(array) == 0:
        raise ValueError("The array is empty.")
    if len(set(array)) != len(array):
        raise ValueError("The array has duplicate values.")

    for i, num in enumerate(array):
        difference = abs(num - value)
        if difference < min_difference:
            min_difference = difference
            closest_index = i
        elif difference == min_difference:
            raise ValueError(f"Multiple occurrences of the closest value in the array. Check indices {closest_index} and {i}")

    return closest_index

############################################################################
def BroadeningLorentzian(arr, broad):
    # arr: [x,y] array to be broadened
    # broad: FWHM
    temp = np.zeros(len(arr[0]))  
    stepsize = np.abs(arr[0][1]-arr[0][0])
    for i in range(len(arr[0])):
        kernel = broad/(2*np.pi) * 1/((arr[0]-arr[0][i])**2 + (broad/2)**2)
        temp[i] = np.dot(kernel,arr[1]) * stepsize
    arr[1] = np.array(temp)
    return arr

def BroadeningGaussian(arr, broad):
    # arr: [x,y] array to be broadened
    # broad: FWHM
    temp = np.zeros(len(arr[0]))
    stepsize = np.abs(arr[0][1]-arr[0][0])
    sigma = broad/(2*np.sqrt(2*np.log(2)))
    for i in range(len(arr[0])):
        kernel = 1/(sigma*np.sqrt(2*np.pi)) * np.exp(-1/2 * (arr[0]-arr[0][i])**2 / sigma**2)
        temp[i] = np.dot(kernel,arr[1]) * stepsize 
    arr[1] = np.array(temp)
    return arr

def EdepBroadeningLorentzian(arr, broadarr):
    # arr: [x,y] array to be broadened
    # broadarr: [FWHM] array with same length as arr[0]
    temp = np.zeros(len(arr[0]))
    stepsize = np.round(np.abs(arr[0][1]-arr[0][0]),6)
    for i in range(len(arr[0])):
        kernel = 0
        for j in range(-3,4):
            kernel += 1/7 * broadarr/(2*np.pi) * 1/(((arr[0]-arr[0][i]+j/7*stepsize))**2 + (broadarr/2)**2)
        temp[i] = (np.sum(kernel*arr[1]) * stepsize)
    arr[1] = np.array(temp)
    return arr
############################################################################

def DS_Function(x, alpha, E, F, Amplitude):
    #print(alpha)
    DS = Amplitude*1/np.pi*(np.cos(np.pi*alpha/2 + (1-alpha)*np.arctan((x-E)/F)))/(F**2 + (x-E)**2)**((1-alpha)/2)
    #DS = Normalize(np.asarray([x,DS]))[1]
    return DS

def DS_FitParameters(arr):
    x = []
    y = []
    for i in range(len(arr[1])):
        x.append(arr[0][i])
        y.append(arr[1][i])
    params, params_covariance = optimize.curve_fit(DS_Function, x, y, p0=[0.02,0.0,0.32,1],
            bounds=([0,-np.inf,0,0],[np.inf,np.inf,np.inf,np.inf]))
    #print("Covariance of parameters: " + str(params_covariance))
    print("Error: " + str(np.sqrt(np.diag(params_covariance))))
    return params

def Lorentzian(x, N, W, center):
    Lorentz = N*W/(W**2 + (x-center)**2)
    return Lorentz

def Gaussian(x, N, W, center):
    #Gauss = N*np.exp(-W(x-center)**2)
    Gauss = N * np.exp(-1*W*(x-center)**2)
    return Gauss

def SkewGaussian(x, SN, SW, Scenter, alpha):
    phi = Gaussian(x, SN, SW, Scenter)
    Phi = 0.5*(1+special.erf(alpha*x/np.sqrt(2)))
    return 2*phi*Phi
  
def Voigt(x,A,cen,sigma,gamma):
    return A*scipy.special.voigt_profile(x-cen,sigma,gamma)

def SkewVoigt(x,A,cen,sigma,gamma,skew):
    v = Voigt(x,A,cen,sigma,gamma)
    s = (1+scipy.special.erf(skew*(x-cen)/(sigma*np.sqrt(2))))
    temp = []
    for i in range(len(v)):
        temp.append(v[i]*s[i])
    return np.asarray(temp)

def SGDouble_Function(x, SN, SW, Scenter, alpha, N, W, center):
    return SkewGaussian(x, SN, SW, Scenter, alpha) + Lorentzian(x, N, W, center)
    #return SkewGaussian(x, SN, SW, Scenter, alpha) + Gaussian(x, N, W, center)

def SG_FitParameters(arr):
    x = []
    y = []
    for i in range(len(arr[1])):
        x.append(arr[0][i])
        y.append(arr[1][i])
    params, params_covariance = optimize.curve_fit(SkewGaussian, x, y, 
            p0=[1.3, 4, 0, 0], bounds=([-np.inf,-np.inf,-np.inf,-np.inf],
                [np.inf,np.inf,np.inf,0]))
    return params

def SGDouble_FitParameters(arr, Nguess, Wguess, centerguess):
    x = []
    y = []
    for i in range(len(arr[1])):
        x.append(arr[0][i])
        y.append(arr[1][i])
    params, params_covariance = optimize.curve_fit(SGDouble_Function, x, y, 
            p0=[1,0.1,0,-15,nguess,wguess,centerguess], bounds=([0,0,-np.inf,-np.inf,0,0,-np.inf],
                [np.inf,np.inf,np.inf,0,np.inf,np.inf,0]))
    return params

def binning(arr,left,right,binnum):
    step = np.round(((right-left)/binnum),6)
    x = np.arange((left+step/2),right,step)
    temp = [[],[]]
    for i in range(len(x)):
        ytemp = 0
        lefttemp = np.round(x[i]-step,8)
        righttemp = np.round(x[i]+step,8)
        for j in range(len(arr[0])):
            if(lefttemp <= arr[0][j] < righttemp):
                ytemp = np.round(ytemp + arr[1][j],8)
        temp[0].append(x[i])
        temp[1].append(ytemp)
    return temp

def matchbinning(main,arr):
    holder = [[],[]]
    temp = main[0]
    for i in range(len(temp)):
        j = 0
        while j < len(arr[0])-1:
            if arr[0][j] < temp[i] < arr[0][j+1]:
                b = arr[1][j]
                m = (arr[1][j+1]-arr[1][j])/(arr[0][j+1]-arr[0][j])
                y = m*(temp[i]-arr[0][j])+b
                holder[0].append(temp[i])
                holder[1].append(y)
                break
            if arr[0][j] == temp[i]:
                holder[0].append(temp[i])
                holder[1].append(arr[1][j])
            j = j+1
    return np.asarray(holder)

def UnitStep(condition, arr):
    for i in range(len(arr[0])):
        if(arr[0][i]>condition):
            arr[1][i] = 1.0
    return arr

def Flipping(arr):
    temp = [[],[]]
    for i in range(len(arr[0])):
        temp[0].append(arr[0][i]*-1)
        temp[1].append(arr[1][i])
    arr = np.asarray(temp)
    arr = arr.T
    arr = arr[np.argsort(arr[:,0])]
    arr = arr.T
    return arr

def xshift(arr, shift):
    for i in range(len(arr[0])):
        arr[0][i] = arr[0][i] + shift
    return arr

def yshift(arr, shift):
    for i in range(len(arr[0])):
        arr[1][i] = arr[1][i] + shift
    return arr

def PlotTrim(arr, left, right):
    temp = [[],[]]
    for i in range(len(arr[0])):
        if left <= arr[0][i] <= right:
            temp[0].append(arr[0][i])
            temp[1].append(arr[1][i])
    return np.asarray(temp)

def alignpeak(arr):
    maxpoint = arr[1][0]
    maxpointpos = arr[0][0]
    for i in range(len(arr[1])):
        if arr[1][i] > maxpoint:
            maxpointpos = arr[0][i]
            maxpoint = arr[1][i]
    temp = [[],[]]
    for i in range(len(arr[0])):
        temp[0].append(np.round((arr[0][i]-maxpointpos),6))
        temp[1].append(arr[1][i])
    return np.asarray(temp)

def PeakMatch(main,arr):
    arr[1] = arr[1]*np.amax(main[1])/np.amax(arr[1])
    return arr

def PeakPosition(arr):
    peakval = arr[1][0]
    peakpos = arr[0][0]
    for i in range(len(arr[1])):
        if arr[1][i] > peakval:
            peakval = arr[1][i]
            peakpos = arr[0][i]
    return [peakpos,peakval]

def FillInXPS(arr):
    i=0
    arr[0] = np.round(arr[0], 6)
    arr[1] = np.round(arr[1], 6)
    stepsize= np.round(np.abs(arr[0][1]-arr[0][0]), 6)
    while i<30:
        if(np.round(np.abs(arr[0][i+1]-arr[0][i]),6) < stepsize):
            stepsize = np.round(np.abs(arr[0][i+1]-arr[0][i]), 6)
        i = i+1
    temp=[[],[]]
    for i in range(len(arr[0])-1):
        temp[0].append(arr[0][i])
        temp[1].append(arr[1][i])
        if(np.round(np.abs(arr[0][i+1]-arr[0][i]), 6) == stepsize):
            skip = True
        else:
            avg = np.round(((arr[1][i]+arr[1][i+1])/2), 6)
            temp[0].append(arr[0][i]+stepsize)
            temp[1].append(avg)
    temp[0].append(arr[0][-1])
    temp[1].append(arr[1][-1])
    return np.asarray(temp)
            
def Normalize(arr):
    arr[1] = arr[1]/integrate.trapz(arr[1],arr[0])
    return np.asarray(np.round(arr,8))

def PeakNormalize(arr):
    arr[1] = arr[1]/np.max(arr[1])
    return arr

def FWHM(arr,peakindex):
        ''' 
        Newer and better way of calculating the FWHM, using linear interpolation. 
            
        This works better for situations where you don't have ton of data points.
        '''
        #Finds half the height of the peak
        half = arr[1][peakindex]/2
        lindex_outer=peakindex
        rindex_outer=peakindex
        #Loops until it hits the right edge of the peak
        while half < arr[1][rindex_outer]:
            rindex_outer = rindex_outer+1
            if half > arr[1][rindex_outer]: 
                rindex_inner = rindex_outer - 1
            elif half == arr[1][rindex_outer]:
                print("In FWHM function there happened to be a discrete point, exactly at half the maximum height")
                rindex_inner = rindex_outer 
            
        #Loops until it hits the left edge of the peak
        while half < arr[1][lindex_outer]:
            lindex_outer = lindex_outer - 1
            if half > arr[1][lindex_outer]: 
                lindex_inner = lindex_outer + 1
            elif half == arr[1][lindex_outer]:
                print("In FWHM function there happened to be a discrete point, exactly at half the maximum height")
                lindex_inner = lindex_outer
        # Find left and right x values by linear interpolation between inner and outer data points
        if rindex_inner == rindex_outer:
            rightval = arr[0][rindex_inner]
        else:
            # Define a line to find approximate x position of halfmax crossing point, y = mx + b
            m = (arr[1][rindex_outer] - arr[1][rindex_inner])/(arr[0][rindex_outer] - arr[0][rindex_inner])
            b = arr[1][rindex_outer] - m*arr[0][rindex_outer]
            rightval = (half - b)/m
            
        if lindex_inner == lindex_outer:
            leftval = arr[0][lindex_inner]
        else:
            # Define a line to find approximate x position of halfmax crossing point, y = mx + b
            m = (arr[1][lindex_inner] - arr[1][lindex_outer])/(arr[0][lindex_inner] - arr[0][lindex_outer])
            b = arr[1][lindex_outer] - m*arr[0][lindex_outer]
            leftval = (half - b)/m
        return np.abs(rightval-leftval)

        ''' Quick and dirty way of finding the FWHM '''
        
        ##Finds half the height of the peak
        #half = arr[1][peakindex]/2
        #lindex=peakindex
        #rindex=peakindex
        ##Loops until it hits the right edge of the peak
        #while half < arr[1][rindex]:
        #    rindex = rindex+1
        ##Loops until it hits the left edge of the peak
        #while half < arr[1][lindex]:
        #    lindex = lindex-1
        #rightval = (arr[0][rindex+1]+arr[0][rindex])/2
        #leftval = (arr[0][lindex-1]+arr[0][lindex])/2
        #return np.abs(rightval-leftval)

def badGREP(file_contents,start_str,start,end,which="first"):
    '''
    Search for the beginning of a particular string in a file
    and return a set of lines in the vicinity of that string. This is
    my shitty approximation of the grep command in linux.

        Paramters
        ---------
        file_contents: list
            list of every line in the file
        start_str : str
            The first few letters (or whole string) that you are looking for in the file
            -> start_index
        start : int
            the first line number, relative to the start_str line number, that you want to keep
            -> [start_index + start:...]
        end : int
            one more than last line number, relative to the start_str line number, that you want to keep
            -> [start_index + start:start_index + end]
        which : str or int
            'fist', 'last', or integer clump that you want to keep

    Examples:
        File Contents:
        1 Hi, my name
        2 is Bob
        3 I love pie
        4
        5 I want to go home
        6 I don't want to
        7 be an example anymore
        8
        9 Please, someone help me

        with open('filename','r') as f:
            file_contents = f.readlines()

        badGREP(file_contents, start_str="is", start=0, end=3)
            -> ["is Bob", "I love pie", ""]

        badGREP(file_contents, start_str="I", start=1, end=3, which="last")
            -> ["be an example anymore", ""]
    '''

    arr = []
    hit_start = False
    for key,line in enumerate(file_contents):
        if line.strip().startswith(start_str):
            arr.append(file_contents[key+start:key+end])
    arr = [[line.replace("\n","") for line in clump] for clump in arr]
    if which == "first":
        return arr[0]
    elif which == "last":
        return arr[-1]
    elif isinstance(which,int):
        return arr[which-1]
    else:
        raise ValueError("which must be defined as either \"first\", \"last\", or an integer.")

def sort_along_x(arr):
    """
    Sorts a numpy array along the x-axis and returns it

    Parameters
    ---------
        arr: numpy array
        Has the format [x_values, y_values]
    """

    if not isinstance(arr, np.ndarray):
        raise ValueError("array for sorting must be a numpy array")
    arr = arr.T
    arr = arr[np.argsort(arr[:,0])]
    arr = arr.T
   
    return arr

def LightenColor(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    mc.cnames.update(mc.BASE_COLORS)
    mc.cnames.update(mc.TABLEAU_COLORS)
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

def Convert_CompName_ToLaTeX(name):
    num_indexes = [m.start(0) for m in re.finditer('(?<=[a-zA-Z])[0-9]+', name)]
    offset = 0
    for ind in num_indexes:
        name = name[0:ind+offset] + "$_\mathdefault{" + name[ind+offset] + "}$" + name[ind+offset+1:]
        name = r"{}".format(name)
        added_character_length = 17 # from '$_\mathdefault{}'
        offset += added_character_length
    return name
