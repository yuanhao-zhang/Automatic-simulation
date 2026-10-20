#!/usr/bin/python
# Script:  domainsize.py
# -------------------------------------------------------

import sys
import numpy as np
from math import *
import fileinput

# INPUT PARAMETERs
nconf =    10
nskip =    80 # 跳过前50帧
M = np.array([0, 1.0, 1.0, 1.0, 1.0])
M_ps = 10000

nbin =     26   
nlayers =  2
mg = nbin * nlayers
half =     True
infiles = ['3_ht.lammpstrj']
file = 'klong.data' if half else 'denprof.csv'

# 数据存储容器
all_configs_densities = np.zeros((nconf, 4, nbin + 1))
all_configs_slopes = np.zeros(nconf)
all_configs_regs = np.zeros(nconf) # 存储每一帧的 0.5 * (reg1 + reg2)

def dennum():
    global reg1, reg2
    for ii in range(1, natoms + 1):
        xr = xc[ii]
        if xr >= xcm1:
            ig = int((xr - xcm1) / drg) + 1
        else:
            ig = int((xr + xbox - xcm1) / drg) + 1
        
        if 1 <= ig <= mg:
            if typea[ii] == 1: dnum1[ig] += 1
            elif typea[ii] == 2: dnum2[ig] += 1
            elif typea[ii] == 3: dnum3[ig] += 1
            elif typea[ii] == 4: dnum4[ig] += 1

        if typea[ii] == 1:
            if 10 <= xr < 20 or 30 <= xr < 40: reg1 += 1
        if typea[ii] == 2:
            if 0 <= xr < 10 or 20 <= xr < 30: reg2 += 1

IN = fileinput.input(infiles)

# 跳过前 nskip 帧
for loopnum in range(0, nskip):
    for i in range(3): IN.readline()
    line = IN.readline()
    natoms = int(line.split()[0])
    for i in range(5 + natoms): IN.readline()

avgLz = 0 

# 处理后续 nconf 帧
for kconf in range(0, nconf):
    IN.readline() 
    line = IN.readline()
    IN.readline() 
    line = IN.readline()
    natoms = int(line.split()[0])
    dim = natoms + 1
    
    if kconf == 0:
        zc, zz, yc, xc = [0]*dim, [0]*dim, [0]*dim, [0]*dim
        typea = [0]*dim

    IN.readline() 
    [xm, xp] = map(float, IN.readline().split())
    [ym, yp] = map(float, IN.readline().split())
    [zm, zp] = map(float, IN.readline().split())
    IN.readline() 
    
    xbox, ybox, zbox = xp - xm, yp - ym, zp - zm
    drg = xbox / mg
    vol = xbox * ybox
    xboxi = xbox / nlayers
    
    cirzcm = cirycm = 0
    
    for j in range(1, dim):
        line = IN.readline()
        fields = line.split()
        k = int(fields[0])
        typea[k] = int(fields[2])
        x1, x2, x3 = float(fields[4]), float(fields[5]), float(fields[6])
        n1, n2, n3 = float(fields[7]), float(fields[8]), float(fields[9])
        
        xc[k] = xbox * x1 + xbox * n1
        yc[k] = ybox * x2 + ybox * n2
        zz[k] = zbox * x3 + zbox * n3

        while xc[k] < 0: xc[k] += xbox
        while xc[k] > xbox: xc[k] -= xbox
        
        xc_relative = xc[k] % xboxi
        if typea[k] == 2:
            theta = (2 * pi / xboxi) * xc_relative
            cirz = (xboxi / (2 * pi)) * sin(theta)
            ciry = (xboxi / (2 * pi)) * cos(theta)
            cirzcm += cirz / M_ps
            cirycm += ciry / M_ps

    thetacm = atan2(-cirzcm, -cirycm) + pi
    xcm1 = (xboxi / (2 * pi)) * thetacm
    avgLz += xboxi / nconf
    
    dnum1, dnum2, dnum3, dnum4 = [0]*(mg+1), [0]*(mg+1), [0]*(mg+1), [0]*(mg+1)
    reg1 = reg2 = 0
    dennum()
    
    zbreal = max(zz) - min(zz)
    binsize = (vol * zbreal) / mg

    # 记录当前帧密度和 Regs
    current_den2 = np.zeros(nbin + 1)
    for ig in range(1, mg + 1):
        bin_id = ((ig - 1) % nbin) + 1
        val2 = float(dnum2[ig]) / (binsize * nlayers)
        all_configs_densities[kconf, 1, bin_id] += val2
        all_configs_densities[kconf, 0, bin_id] += float(dnum1[ig]) / (binsize * nlayers)
        current_den2[bin_id] += val2

    # 每一帧的梯度计算 (基于对称化)
    if half:
        halfbin = int(nbin/2) + 1
        den2_h = (current_den2[halfbin:] + current_den2[1:halfbin][::-1]) / 2
        z_h = (np.arange(len(den2_h)) + 0.5) / (2 * len(den2_h))
        cur_slopes = np.gradient(den2_h, z_h)
        all_configs_slopes[kconf] = cur_slopes[int(len(cur_slopes)/2)]

    # 每一帧的 0.5 * (reg1 + reg2)
    all_configs_regs[kconf] = 0.5 * (reg1 + reg2) / M_ps

# --- 最终统计计算 ---
ddnum1_avg = np.mean(all_configs_densities[:, 0, :], axis=0)
ddnum2_avg = np.mean(all_configs_densities[:, 1, :], axis=0)
stderr1 = np.std(all_configs_densities[:, 0, :], axis=0) / np.sqrt(nconf)
stderr2 = np.std(all_configs_densities[:, 1, :], axis=0) / np.sqrt(nconf)

slope_avg = np.mean(all_configs_slopes)
slope_err = np.std(all_configs_slopes) / np.sqrt(nconf)

reg_avg = np.mean(all_configs_regs)
reg_err = np.std(all_configs_regs) / np.sqrt(nconf)

# --- 文件输出 ---
OUT = open(file, 'w')
if half:
    halfbin = int(nbin/2) + 1
    ddnum1 = (ddnum1_avg[halfbin:] + ddnum1_avg[1:halfbin][::-1]) / 2
    ddnum2 = (ddnum2_avg[halfbin:] + ddnum2_avg[1:halfbin][::-1]) / 2
    err1 = np.sqrt(stderr1[halfbin:]**2 + stderr1[1:halfbin][::-1]**2) / np.sqrt(2)
    err2 = np.sqrt(stderr2[halfbin:]**2 + stderr2[1:halfbin][::-1]**2) / np.sqrt(2)
    nbin_final = len(ddnum1)

    OUT.write("# z, den1, err1, den2, err2\n")
    for i in range(nbin_final):
        z_coord = (i + 0.5) / (2 * nbin_final)
        OUT.write("%1.4f %8.9f %8.9f %8.9f %8.9f\n" % (z_coord, ddnum1[i], err1[i], ddnum2[i], err2[i]))
    
    # 最后一行依次输出要求的四个统计值
    OUT.write("%8.9f %8.9f %8.9f %8.9f\n" % (slope_avg, slope_err, reg_avg, reg_err))

OUT.close()
print('Processing complete. Output written to %s' % file)
