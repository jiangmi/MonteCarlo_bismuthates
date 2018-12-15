import commands
import shutil
import os
import sys
import time
import linecache
import math
import numpy as np

do_submit = True

Nxs = [4]  # e.g. 5 supercell = 10x10 lattice
Nys = Nxs
mus = np.arange(-0.14, -0.131, 0.1)#, 1.2, 1.1, 1]#, 2.4, 2.6, 2.8, 3, 3.2, 3.4]#, 3.6, 3.8, 4, 4.2, 4.4]

es = 6.42
ep = 2.42
tsp = 2.08 
tpp = 0.056
dXamp = 0.2

beta_max = 10.0
beta_min = 10.0
num_beta_steps = 1

nwarm = 500
ninv = 10    # print warmup progress per ninv steps
nmeas = 2000
nbin = 10
alpha = 300
spring_const = 10

def prepare_file(Nx, mu, fname, dir):

    file = open("parameters.f90", "r")
    text = file.read()
    file.close()

    text = text.replace("Nxval",  str(Nx))
    text = text.replace("Nyval",  str(Nx))
    text = text.replace("muval",  str(mu))
    text = text.replace("esval",  str(es))
    text = text.replace("epval",  str(ep))
    text = text.replace("tspval",  str(tsp))
    text = text.replace("tppval",  str(tpp))
    text = text.replace("dXampval",  str(dXamp))
    text = text.replace("bemaxval",  str(beta_max))
    text = text.replace("beminval",  str(beta_min))
    text = text.replace("nbetaval",  str(num_beta_steps))
    text = text.replace("nwarmval",  str(nwarm))
    text = text.replace("nmeasval",  str(nmeas))
    text = text.replace("ninvval",  str(ninv))
    text = text.replace("nbinval",  str(nbin))
    text = text.replace("alphaval" ,  str(alpha))
    text = text.replace("springval",  str(spring_const))
    text = text.replace("fnameval",  str(fname))
    text = text.replace("dirval",  str(dir))

    file = open("parameters.f90", "w")
    file.write(text)
    file.close()

for Nx in Nxs:
    for mu in mus:
        print "Nx = ", Nx, "mu = ", mu
        fname = 'Nx'+str(Nx)+'_mu'+str(mu)+'_es'+str(es)+'_ep'+str(ep)+'_tsp'+str(tsp)+'_tpp'+str(tpp)+'_alpha'+str(alpha)+'_c'+str(spring_const)+'_dX'+str(dXamp)
        dir = "/home/mijiang/bismuth_MC/Nx"+str(Nx)+"_mu" + str(mu)

        os.chdir('/home/mijiang/bismuth_MC')
        cmd = "cp parameters_tmp.f90 parameters.f90"
        os.system(cmd)
        prepare_file(Nx,mu,fname,dir)

        cmd = "ifort -O3 main.f90 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
        os.system(cmd)
        cmd = "cp batch_script.slm job.slm"                                              
        os.system(cmd)                 

        file = open("job.slm", "r")
        text = file.read()
        file.close()                
        text = text.replace("fnameval",  str(fname))
        text = text.replace("dirval",  str(dir))
        file = open("job.slm", "w")
        file.write(text)
        file.close()
    
        if not os.path.exists(dir):
          cmd = "mkdir " + dir
          os.system(cmd)
          cmd = "cp a.out job.slm " + dir
          os.system(cmd) 

        if(do_submit):
          os.chdir(dir)
          cmd = "sbatch job.slm"
          os.system(cmd)

