import commands
import shutil
import os
import sys
import time
import linecache
import math
import numpy as np

do_submit = True

Nxs = [10]  # e.g. 5 supercell = 10x10 lattice
Nys = Nxs

# travelling cluster size:
Ncx = 2
Ncy = 2

#mus = np.arange(0.1, -0.151, -0.01)#, 1.2, 1.1, 1]#, 2.4, 2.6, 2.8, 3, 3.2, 3.4]#, 3.6, 3.8, 4, 4.2, 4.4]
mus = [-0.22, -0.24, -0.26, -0.28, -0.3, -0.2, -0.19, -0.18, -0.17, -0.16, -0.15, -0.14, -0.13, -0.12, -0.11, -0.1, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0.0, \
       0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2]
#mus = [0.1]#, 0.0, -0.22]#, -0.24, -0.26, -0.28, -0.3]
mus = [-0.2, -0.18, 0.18, 0.0, 0.1]#, -0.18, -0.16, -0.14, -0.12, -0.1, -0.08, -0.06, -0.04, -0.02, 0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12]

es = 6.42
ep = 2.42
tsp = 2.08 
tpp = 0.056
dXamp = 0.05

nbeta = 3
betas = '1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 200.0, 300.0, 500.0, 1000.0'
betas = '3.0, 10.0, 100.0'

nwarm = 1000
ninv = 10    # print warmup progress per ninv steps
nmeas = 20000
measinv = 20
nbin = 10
ks = [65]#, 50, 60, 70, 80]#1, 10, 30, 50, 70, 90]
alphas = [800]#, 1500]#, 500, 600, 700]

def prepare_file(Nx, mu, k, al, fname, dir):

    file = open("parameters.f90", "r")
    text = file.read()
    file.close()

    text = text.replace("Nxval",  str(Nx))
    text = text.replace("Nyval",  str(Nx))
    text = text.replace("Ncxval",  str(Ncx))
    text = text.replace("Ncyval",  str(Ncy))
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
    text = text.replace("measinvval",  str(measinv))
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
                fname = 'Nx'+str(Nx)+'_mu'+str(mu)+'_es'+str(es)+'_ep'+str(ep)+'_tsp'+str(tsp)+'_tpp'+str(tpp)+'_k'+str(k)+'_al'+str(al)+'_dX'+str(dXamp)
                dir = "/home/mijiang/bismuth_MC/test/"+fname

                os.chdir('/home/mijiang/bismuth_MC/test')
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
                text = text.replace("Nxval",  str(Nx))

                file = open("job.slm", "w")
                file.write(text)
                file.close()

                if os.path.exists(dir):
                  cmd = "rm -r " + dir

                cmd = "mkdir " + dir
                os.system(cmd)
                cmd = "cp a.out job.slm parameters.f90 " + dir
                os.system(cmd) 
                cmd = "cp ./Xreadin/X_readin_mu"+str(mu)+"* "+dir
                os.system(cmd)

                if(do_submit):
                  os.chdir(dir)
                  cmd = "sbatch job.slm"
                  os.system(cmd)

