# -*- coding: utf-8 -*-
"""
Created on Sun Nov  9 15:36:25 2025

@author: abram
"""
#########importing########

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

########End import########

#############functions############

#Making a model of function of gaussian function + linear line function
#which can be used in curve_fit function to find the best fit curve
#this function used core worksheet3 solution as reference
def fit(x, a, mu, s, m, c):
      gau = a * np.exp(-(x-mu)**2/(2*s**2)) #gaussian function
      y = m * x + c #linear function
      z=y+gau #function of gaussian distribution+linear
      return z
  
#function used to calculated the redshift velocity with formula of redshift
#frequency observed is the input value
def velocity(freq):
    c_light=299792458 #speed of light
    lam=c_light/freq #wavelength = speed of light/frequency
    ratio_lam=(lam/(486.1*10**-9)) #ratio between emitted and observed freqency
    v=((ratio_lam**2-1)/(ratio_lam**2+1))*c_light #redshift velocity
    return v

###########End functions##########

#######################Main codes(A,B)###########################

###A.plotting best fit curve for noisy data and estimating the velocities

df=pd.read_csv(r"C:/Users/abram/OneDrive - Imperial College London/yr1/Yr1 computing/project1/Data/SpectralData_Hbeta.csv", skiprows=4) 
#data of observed wavelengths and intensities are read from excel 
#and header is skipped 
v_list=[]
err_v=[]
c_light=299792458
#empty list for velocity observed and error of velocity

for i in range(0, df.shape[1], 2): 
#storing coloumns of data every 2 coloumn(each observation)  
#looping through every coloumn, there are 4 part inside the loop

#1.frequency and intensity for each observation are stored into lists
      f_list=df.iloc[:,i].values #frequency 
      int_list=df.iloc[:,i+1].values #intensity
      
#2.initial parameters for curve_fit function are estimated 
      m, c = np.polyfit(f_list, int_list, 1) 
      #linear line of best fit for the noisy data
      max_int = max(int_list) 
      #value of max intensity
      x=np.where(int_list == max_int)[0][0] 
      #index of where max intensity is in the list of intensity
      y=m*f_list[x]+c 
      #value of max intensity frequency on the linear line 
      #(for finding the amplitude)     
      plt.figure() #seperating figures for each observation
      amp=max_int-y
      #estimation of amplitude of gaussian distribution
      #(max intensity-intensity of same frequency on the linear best fit line)
      plt.plot(f_list,int_list, label='data') 
      #plot of the noisy data
      sd_guess=(max(f_list)-min(f_list))/10 
      #estimation of standard deviation(suggested by chatgpt)
      guess=[amp,f_list[x], sd_guess,m,np.min(f_list)]
      #initial guessed parametre for the curve_fit function

#3.Curve fitting with scipy curve_fit,and plot the fitting curve  
   #this function used core worksheet3 solution as reference   
      po, po_cov = curve_fit(fit, f_list, int_list, guess) 
      #curve fitting
      plt.plot(f_list, fit(f_list, po[0], po[1], po[2], po[3], po[4]),
               label='fitting curve')
      #plot of the fitting curve to check the validity
      plt.xlabel("Frequency(hz)")#labeling x-axis
      plt.ylabel("intensity(A.U.)")#labeling y-axis
      plt.legend()#legend
      #plt.savefig(f"C:/Users/abram/OneDrive - Imperial College London/yr1/Yr1 computing/project1/spectral_fit_{i//2 + 1:02d}.png",
            #dpi=300, bbox_inches="tight")
      plt.close()
#4.estimating the mean 
      v = velocity(po[1])
      #estimateing velocity 
      #by inputting mean observed frequency of gaussian distribution graph
      #to the redshift velocity function at the top
      v_list.append(v/1000)
      #storing velocity into a list(1/1000 changes the unit from m/s to km/s)
      #estimating error of velocity:
      sigma_freq = np.sqrt(po_cov[1,1])
      wavelength_0 = 486.1e-9
      r = c_light / (po[1] * wavelength_0)
      sigma_v = abs(-4*c_light*r**2 /(po[1]*(r**2+1)**2)) * sigma_freq
      err_v.append(sigma_v/1000)

###B.eliminating unwanted data,
###  matching distance to velocity with the same observation number
dis_df=pd.read_csv(r'C:/Users/abram/OneDrive - Imperial College London/yr1/Yr1 computing/project1/Data/dist_data.txt',sep='\t',skiprows=5)
#data of distance for each observation number are read and header is skipped
ob = df.columns[::2].astype(int).tolist()
#observation numbers from the excel file are stored into a list
#the order of observation number is important in matching the data later
#every 2 number are stored as each number repeat twice
dist=[]
err_vplot=np.array([])
v_plot=[]
#empty list of distance and velocity for plotting
#empty array for error of velocity
ins = dis_df['Instrument_response']
distance=dis_df['Distance']
obs=dis_df['Observation_number']
#storing row of instrument response, distance and observation number 

for i in range(0,dis_df.shape[0],1): 
#looping through every row
    if ins[i]=="E":
        continue
    #ignore the data if E appears in the instrument response
    else:
       dist.append(distance[i])  
       #storing the value of distance into list for plotting
       idx = ob.index(int(obs.iloc[i]))
       #index for where the observation number in the excel file matches 
       #the one in the distance dataframe
       v_plot.append(v_list[idx])
       #storing the matched velocity from excel file 
       #to list of velocity for plotting
       err_vplot=np.append(err_vplot,float(err_v[idx]))
       #storing the matched error of velocity into list
plt.figure()
plt.scatter(dist,v_plot, color='black')
#distance against velocity
plt.xlabel('distance(Mpc)')
plt.ylabel('Redshift(km/s)')

dist = np.array(dist, dtype=float)
v_plot = np.array(v_plot, dtype=float)
err_vplot = np.array(err_vplot, dtype=float)

h, e = np.polyfit(dist,v_plot,1,w=1/err_vplot,cov=True)
err_h = np.sqrt(e[0,0])
#find the error of hubble constant and the gradient of best fit line
plt.plot(dist,h[0]*dist + h[1],color='red',linewidth=2,label='Weighted fitting line')
#plotting best fit line
plt.legend()#legend
plt.show()
print(f'Estimate of Hubble"s constant: {h[0]:.2f} +/- {err_h:.2f} km/s/Mpc')  
#stating the final result 

#########################End Main##########################
