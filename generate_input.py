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
mus = np.arange(0.0, 1.001, 0.05)#, 1.2, 1.1, 1]#, 2.4, 2.6, 2.8, 3, 3.2, 3.4]#, 3.6, 3.8, 4, 4.2, 4.4]

es = 6.42
ep = 2.42
tsp = 2.08 
tpp = 0.056
dXamp = 0.15

nbeta = 1
betas = '10.0'

nwarm = 100000
ninv = 10    # print warmup progress per ninv steps
nmeas = 500000
nbin = 10
ks = [35]#, 50, 60, 70, 80]#1, 10, 30, 50, 70, 90]
alphas = [400, 500, 600, 700]

def prepare_file(Nx, mu, k, al, fname, dir):

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
    text = text.replace("nbetaval",  str(nbeta))
    text = text.replace("betasval",  '(/ '+str(betas)+' /)')
    text = text.replace("nwarmval",  str(nwarm))
    text = text.replace("nmeasval",  str(nmeas))
    text = text.replace("ninvval",  str(ninv))
    text = text.replace("nbinval",  str(nbin))
    text = text.replace("alphaval" ,  str(al))
    text = text.replace("springval",  str(k))
    text = text.replace("fnameval",  str(fname))
    text = text.replace("dirval",  str(dir))

    file = open("parameters.f90", "w")
    file.write(text)
    file.close()

for Nx in Nxs:
    for mu in mus:
        for k in ks:
            for al in alphas:
                print "Nx = ", Nx, "mu = ", mu, 'k = ', k, 'alpha = ', al
                fname = 'Nx'+str(Nx)+'_mu'+str(mu)+'_es'+str(es)+'_ep'+str(ep)+'_tsp'+str(tsp)+'_tpp'+str(tpp)+'_k'+str(k)+'_alpha'+str(al)+'_dX'+str(dXamp)
                dir = "/home/mijiang/bismuth_MC/run/"+fname

                os.chdir('/home/mijiang/bismuth_MC/run')
                cmd = "cp parameters_tmp.f90 parameters.f90"
                os.system(cmd)
                prepare_file(Nx,mu,k,al,fname,dir)

                cmd = "ifort -O3 main.f90 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
                os.system(cmd)
                cmd = "cp batch_script.slm job.slm"                                              
                os.system(cmd)                 

                file = open("job.slm", "r")
                text = file.read()
                file.close()                
                text = text.replace("fnameval",  str(fname))
                text = text.replace("dirval",  str(dir))
		text = text.replace("muval",  str(mu))
		text = text.replace("kval",  str(k))
		text = text.replace("alval",  str(al))

                file = open("job.slm", "w")
                file.write(text)
                file.close()

                if not os.path.exists(dir):
                  cmd = "mkdir " + dir
                  os.system(cmd)
                  cmd = "cp a.out job.slm parameters.f90 " + dir
                  os.system(cmd) 

                if(do_submit):
                  os.chdir(dir)
                  cmd = "sbatch job.slm"
                  os.system(cmd)

