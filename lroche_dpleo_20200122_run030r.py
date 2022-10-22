#PhD reseach 2020
#Kittipong Wangnok, D6010218
#School of Physics, Institute of Science, Suranaree University of Technology

#Import all module
import sys
import os
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
from statistics import stdev
from statistics import mean
np.seterr(divide='ignore', invalid='ignore')

#Latex font
import matplotlib as mpl
from matplotlib import rc
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=16)

#################################################################################
'''
1. Input file: lcurve_dpleo_data
'''
#################################################################################
#Please change the input file
lcurve_dpleo_data = open("dpleo_20200122_run030r.dat",'r').readlines()
N_lcurve_dpleo_data = len(lcurve_dpleo_data)

dat_MJD_time = []
dat_MJD_time_err = []
dat_Flux = []
dat_Flux_err = []

for line in open("dpleo_20200122_run030r.dat"):
    li=line.strip()
    if not li.startswith("#"):
        dat_MJD_time.append(float(li.split(" ")[0]))
        dat_MJD_time_err.append(float(li.split(" ")[1]))
        dat_Flux.append(float(li.split(" ")[2]))
        dat_Flux_err.append(float(li.split(" ")[3]))

data_result = []
for i in range (len(lcurve_dpleo_data)):
#    print ('%0.6f\t%0.6f\t%0.6f\t%0.6f' %(dat_MJD_time[i],dat_MJD_time_err[i],dat_Flux[i],dat_Flux_err[i]))
    data_result.append('%0.6f\t%0.6f\t%0.6f\t%0.6f' %(dat_MJD_time[i],dat_MJD_time_err[i],dat_Flux[i],dat_Flux_err[i]))
    
dat = data_result
f = open('data_dpleo_20200122_run030r.txt', 'w')
#for upper_result in upper_result:
for i in range(len(dat)):
    f.write(str(dat[i])+ '\n')
f.close()

#################################################################################
'''
2. Input file: lcurve_dpleo_output
'''
#################################################################################
#Please change the input file
lcurve_dpleo_output = open("dpleo_20200122_run030r.out",'r').readlines()
N_lcurve_dpleo_output = len(lcurve_dpleo_output)

out_MJD_time = []
out_MJD_time_err = []
out_Flux = []
out_Flux_err = []

for line in open("dpleo_20200122_run030r.out"):
    li=line.strip()
    if not li.startswith("#"):
        out_MJD_time.append(float(li.split(" ")[0]))
        out_MJD_time_err.append(float(li.split(" ")[1]))
        out_Flux.append(float(li.split(" ")[2]))
        out_Flux_err.append(float(li.split(" ")[3]))
        
output_result = []
Res = []
residual_result = []
Chi_sqr_a = [i for i in range(len(lcurve_dpleo_data))]

for i in range (len(lcurve_dpleo_data)):
#    print (out_MJD_time[0])
#    print ('%0.6f\t%0.6f\t%0.6f\t%0.6f' %(out_MJD_time[i],out_MJD_time_err[i],out_Flux[i],out_Flux_err[i]))
    output_result.append('%0.6f\t%0.6f\t%0.6f\t%0.6f' %(out_MJD_time[i],out_MJD_time_err[i],out_Flux[i],out_Flux_err[i]))
    Res = dat_Flux[i] - out_Flux[i]
    residual_result.append('%0.6f' %(Res))
    Chi_sqr = (dat_Flux[i] - out_Flux[i])**2/out_Flux[i]
    Chi_sqr_a[i] = Chi_sqr
    chisq = sum(Chi_sqr_a)
#    print ('%0.6f' %(Res))
#    print ('%0.6f' %(Chi_sqr))
#print ('%0.2f' %(chisq))
'''
3. Save the output
'''

out = output_result
f = open('output_dpleo_20200122_run030r.txt', 'w')
#for upper_result in upper_result:
for i in range(len(out)):
    f.write(str(out[i])+ '\n')
f.close()

res = residual_result
f = open('residual_dpleo_20200122_run030r.txt', 'w')
#for upper_result in upper_result:
for i in range(len(res)):
    f.write(str(res[i])+ '\n')
f.close()

#################################################################################
'''
4. Plot the result
'''
#################################################################################
InputFile_1 = "data_dpleo_20200122_run030r.txt"
Data_1   = np.genfromtxt(InputFile_1)

InputFile_2 = "output_dpleo_20200122_run030r.txt"
Data_2   = np.genfromtxt(InputFile_2)

InputFile_3 = "residual_dpleo_20200122_run030r.txt"
Data_3   = np.genfromtxt(InputFile_3)

#Read data
E = out_MJD_time[0]

MJD_time_1 = Data_1[:,0] - E
Flux_1 = Data_1[:,2]
Flux_err_1 = Data_1[:,3]

MJD_time_2 = Data_2[:,0] - E
Flux_2 = Data_2[:,2]
Flux_err_2 = Data_2[:,3]

Res = Flux_1 - Flux_2
#print (Res)

fig, (ax0, ax1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1]}, sharex=True, sharey=False, figsize=(8, 7), tight_layout=True)
plt.xlim(MJD_time_1[0], MJD_time_1[-1])
#plt.xlim(0.6405, 0.646)
plt.xlabel('MJD - '+str(E))

ax0.tick_params(direction='in', which='both', bottom='on',top='on', right = 'on')
ax1.tick_params(direction='in', which='both', bottom='on',top='on', right = 'on')


ax0.errorbar(MJD_time_1, Flux_1, yerr=Flux_err_1, fmt='o', markersize = 4, alpha = 0.25, color='red', label = 'data\_20200122\_run030r')
ax0.plot(MJD_time_2, Flux_2, 'b-', label='model\_fitting, $\chi^2$ = '+str('%0.2f' %(chisq)))
ax0.legend(loc="best")
ax0.set_ylabel('Relative flux')

ax1.errorbar(MJD_time_1, Res, yerr=Flux_err_1, fmt='o',markersize = 4, color='red', alpha = 0.25, label = 'data\_star24')
ax1.hlines(y=0, xmin=MJD_time_1[0], xmax=MJD_time_1[-1], colors='black', linestyles='-')
ax1.set_ylabel('Residual')

fig.align_ylabels()
output_filename = os.path.splitext(__file__)[0] + '.png'
plt.savefig(output_filename, dpi=1000)
plt.show()

