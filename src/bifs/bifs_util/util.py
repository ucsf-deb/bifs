###                  Utilities for BIFS                     ###
###                                                         ###
###  A set of utility functions for reading a jsonpickled   ###
###  BIFS object, examing the components of a BIFS object   ###
###  (e.g. image, transform representation, likelihood and  ###
###  prior at a voxel) during various stages                ###
###  of a BIFS analysis.                                    ###

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.stats as stats
from scipy.optimize import minimize_scalar as ms
import jsonpickle

def read_bifs(filename):
    bifs_handle = open(filename,"r")
    bifsj = bifs_handle.read()
    bifs_handle.close()
    return(jsonpickle.decode(bifsj))
 
def voxel_dist(bifs_obj,vox_coord,do_plots=False):
    """
    voxel_dist() is a function for examining particular 
    choices of prior, likelihood, and posterior distributions at 
    a particular transform space point. It uses a Gaussian prior
    and focuses on illustrating the differences between the results
    of using the analytic form for a conjugate Gaussian likelihood and
    prior, or just directly calculating, via optimaization, use of a
    Gaussian or Rician likelihood. It outputs a number of parameters
    (e.g. distribution means and MAP estimates) and optionally plots
    the distributions, and assumes that an image has been loaded.
    
    Inputs:

    bifs_obj - a bifs object containing the Fourier representation
               of an image (i.e. an image has to have been loaded)

    vox_coord - an array containing the 1D, 2D, or 3D coordinates
                of a transform space voxel at which you want to
                examine the the prior, likelihood, and posterior
                distributions and some of their properties (e.g. mean)

    do_plots = a flag determining whehter or not to display the
               distributions (True,False); if False will just 
               print mean and scale for Gaussian prior and , Gaussian 
               and Rician likelihoods (only distributions currently used 
               for this test)

    Output:

    Prints  -  Prior mean
               Prior scale
               K-space data value
               Likelihood scale,
               Maximum of Rice Likelihood
               Maximum of Gaussian Likelihood, 
               Conjugate Gaussian value, 
               Gaussian Likelihood MAP estimate,
               Rician Likelihood MAP estimate
               k-space co-ordinate

    Optionally plots - Prior Distribution,
                       Gaussian Likelihood Distribution,
                       Rician Likelihood Distribution,
                       Posterior with Gaussian Likelihood,
                       Posterior with Rician Likelihood
        
    """
    # Normalization options:
    # 1. Match power in prior to power in image
    # 2. No normalization
    # norm = (np.sum(mybifs.mod_image))/(np.sum(mybifs.prior_mean))
    norm = 1.0
        
    prior_mn = norm*bifs_obj.prior_mean[tuple(vox_coord)]
    prior_sc = norm*bifs_obj.prior_std[tuple(vox_coord)]
    data_mn = bifs_obj.mod_image[tuple(vox_coord)]
    data_sc = bifs_obj.likelihood_scale*data_mn
    
    # Shape parmeter for Rician
    b = 1.0/bifs_obj.likelihood_scale
    
    # Guassian conjugate prior case
    ksd2 = data_sc**2
    priorsd2 = prior_sc**2
    gauss_guass = ((prior_mn/priorsd2+data_mn/ksd2)/(1./priorsd2 + 1./ksd2))
        
    err_code = 0
    if np.size(vox_coord) != bifs_obj.imdim:
        print("Co-ordinate dimension:",np.size(vox_coord),"is not the same as image dimension:",bifs_obj.imdim)
        err_code = 1
    else:
        print("Prior mean: ",prior_mn)
        print("Prior scale: ",prior_sc)
        print("K-space data value:",data_mn)
        print("Likelihood scale:",data_sc)

        # Distributions
        glike = lambda x: stats.norm.pdf(x,loc=data_mn,scale=data_sc)
        rlike = lambda x: stats.rice.pdf(x,b,scale=data_sc)
        prior = lambda x: stats.norm.pdf(x,loc=prior_mn,scale=prior_sc)
            
        findmin = lambda x: -rlike(x)
        res = ms(findmin, method='brent')
        print("Maximum of Rice Likelihood at: ",res.x)

        findmin = lambda x: -glike(x)
        res = ms(findmin, method='brent')
        print("Maximum of Gaussian Likelihood at: ",res.x)

        gposterior = lambda x: prior(x)*glike(x)
        rposterior = lambda x: prior(x)*rlike(x)

        print("Conjugate Gaussian value: ",gauss_guass)

        findming = lambda x: -gposterior(x)
        resg = ms(findming, method='brent')
        print("Gaussian Likelihood MAP estimate:",resg.x,"at k-space co-ordinates:",tuple(vox_coord))

        findminr = lambda x: -rposterior(x)
        resr = ms(findminr, method='brent')
        print("Rician Likelihood MAP estimate:",resr.x,"at k-space co-ordinates:",tuple(vox_coord))
        
        if do_plots:
            n_points = 100 # points to plot
            
            plot_start = prior_mn - 4.*prior_sc
            plot_end = prior_mn + 4.*prior_sc
            params_prior = np.linspace(plot_start,plot_end, n_points)
            plot_prior = np.array([prior(p) for p in params_prior])
            
            plot_start = data_mn - 4.*data_sc
            plot_end = data_mn + 4.*data_sc
            params_glike = np.linspace(plot_start,plot_end, n_points)
            plot_glike = np.array([glike(p) for p in params_glike])
            
            params_rlike = np.linspace(stats.rice.ppf(0.01,b,scale=data_sc),stats.rice.ppf(0.99,b,scale=data_sc), 100)
            plot_rlike = np.array([rlike(p) for p in params_rlike])
                
            plot_start = data_mn - 4.*data_sc
            plot_end = data_mn + 4.*data_sc
            params_glike_post = np.linspace(plot_start,plot_end, n_points)
            plot_glike_post = np.array([gposterior(p) for p in params_glike_post])
    
            plot_start = data_mn - 4.*data_sc
            plot_end = data_mn + 4.*data_sc
            params_rlike_post = np.linspace(plot_start,plot_end, n_points)
            plot_rlike_post = np.array([rposterior(p) for p in params_rlike_post])

            fig, axes = plt.subplots(5, 1, sharex=True, figsize=(8,8))
            axes[0].plot(params_prior, plot_prior)
            axes[0].set_title("Prior Distribution")
            axes[1].plot(params_glike, plot_glike)
            axes[1].set_title("Gaussian Likelihood Distribution")
            axes[2].plot(params_rlike, plot_rlike)
            axes[2].set_title("Rician Likelihood Distribution")
            axes[3].plot(params_glike_post, plot_glike_post)
            axes[3].set_title("Posterior with Gaussian Likelihood")
            axes[4].plot(params_rlike_post, plot_rlike_post)
            axes[4].set_title("Posterior with Rician Likelihood")
            plt.tight_layout()
            plt.show()
            plt.close("all")
    return err_code

def plot_param_func(bifs_obj):
    """
    plot_param_func() plots the current transform space parameter 
    function given a bifs object.

    Inputs:

    bifs_obj - a bifs object containing the Fourier representation
               of an image (i.e. an image has to have been loaded)

    Outputs:

    Plots current paramter function

    """

    if bifs_obj.imdim == 1 or bifs_obj.imdim == 3:
        plt.title("K Space Parameter Function") 
        plt.xlabel("kx")
        KX = np.arange(bifs_obj.imdim1)
        # For now have to do the following by hand
        # Should figure some way to automate
        # when new parameter functions are added
        if bifs_obj.param_func_type == "Inverse Power Decay":
            Z = bifs_obj.bas.ixsc(bifs_obj.bvec,KX,bifs_obj.decay)
        elif bifs_obj.param_func_type == "Banded Inverse Power Decay":
            Z = bifs_obj.bas.ixscbanded(bifs_obj.bvec,KX,bifs_obj.decay,bifs_obj.banded_cutoff)
        elif bifs_obj.param_func_type == "Linear Decay":
            Z = bifs_obj.bas.linsc(bifs_obj.bvec,KX) 
        else:
            pass
        if bifs_obj.bumps:
            Z += bifs_obj.bas.add_bumps(bifs_obj.bvec,bifs_obj.kdist,bifs_obj.bumps,np.max(bifs_obj.kdist))
        plt.plot(KX,Z)
    else: # Sort of silly since the parameter function is rotationally
          # symmetric (unless there's a trivial deviation due to kx != ky)
          # but it's pretty and might help check numerical glitches in filters
        fig = plt.figure()    
        ax = fig.gca(projection='3d')
        ax.set_title('K Space Parameter Function')
        ax.set_xlabel('kx')
        ax.set_ylabel('ky')
        KX = np.arange(bifs_obj.imdim1)
        KY = np.arange(bifs_obj.imdim2)
        KX, KY = np.meshgrid(KX, KY)
        Kdist = np.sqrt(KX**2 + KY**2)
        if bifs_obj.param_func_type == "Inverse Power Decay":
            Z = bifs_obj.bas.ixsc(bifs_obj.bvec,Kdist,bifs_obj.decay)
        elif bifs_obj.param_func_type == "Banded Inverse Power Decay":
            Z = bifs_obj.bas.ixscbanded(bifs_obj.bvec,Kdist,bifs_obj.decay,bifs_obj.banded_cutoff)
        elif bifs_obj.param_func_type == "Linear Decay":
            Z = bifs_obj.bas.linsc(bifs_obj.bvec,Kdist) 
        else:
            do_nothing = None
        if bifs_obj.bumps:
            Z += bifs_obj.bas.add_bumps_to_pf(bifs_obj.bvec,bifs_obj.kdist,bifs_obj.bumps,np.max(bifs_obj.kdist))
        ZRoll = np.roll(np.roll(Z,(Z.shape[0]//2),0),(Z.shape[1]//2),1)
        surf = ax.plot_surface(KX, KY, ZRoll, cmap=cm.coolwarm,linewidth=0, antialiased=False)

    plt.show()
    plt.close('all')
    return
