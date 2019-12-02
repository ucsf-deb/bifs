
###                        BIFS class                   ### 
###                                                     ###
###         Class for performing Bayesian Image         ###
###         Restoration in Fourier Space (BIFS)         ###

# for debugging
import traceback


import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from pylab import *
import multiprocessing as mp
from scipy import stats
import imageio
from scipy.optimize import minimize_scalar as ms
from multiprocessing import Pool, TimeoutError
from datetime import datetime
import jsonpickle as jsp

class bifs:
    """
    Class bifs for performing Bayesian image restoration
    in k-space

     Class Variables
    ----------------

    Class variables required by constructor:

    Class variables available to constructor:

    init_image - initial loaded image
    k_image - initial k-space image
    mod_image - initial modulus image in k-space
    phase_image - initial phase image in k-space
    bifsk_image - final BIFS modulus image in k-space
    final_image - final reconstructed image in image space

    NOTE: We're keeping all these images in the BIFS object for
    now re. testing and experimentation - might be useful to 
    eventually have production run options that are more
    parsimonius re. storage.

    image_file_loaded - whether an image is loaded (True,False)
    initial_image_file_name - file name of initial image
    imdim - int image dimension (1,2 or 3)

    kdist = distance funcion on the shifted k-space lattice

    view3Dslice - for 3D data this is a 2D array [a,b] where:
                  a = axis perpindicular to slice
                  b = fraction of maximum along that direction 
                      for slice location
    
    prior - string specifying prior distribution to use
            current choices:
            'Gaussian'

    prior_choices - list of current prior choices (see above)
    prior_mean_init - prior mean before paramter space function
                      is set up (used for tests)
    prior_mean - the prior mean defined at each k-space point 
                 by the k-space parameter function
    prior_std - the prior std defined at each k-space point
    prior_scale - the overall scale of the prior variance
    prior_scale_orig - prior scale at the origin - generally set huge
                           to allow data to determine overall scale

    likelihood - string specifying likelihood distribution to use
                 current choices:
                 'Gaussian'
                 'Rician'

    likelihood_choices - list of current choices (see above)
    likelihood_scale - the assumed (const) noise level in k-space

    bessel_approx_lims - limits for bessel approximtion for rice
                         distribution - see paper referenced in code

    bessel_approx_array - array for bessel approximtion for rice
                         distribution - see paper referenced in code
    
    rice_denom_cutoff - cutoff for the demoninator of the closed form 
                        of the posterior with a Gaussian prior and
                        Rician likelihood derived from bessel approximation
                        see paper referenced in code

    param_func_type - string specifying the k-space BIFS paramter
                      function to use
                      current choices:
                      "Inverse Power Decay"
                      "Banded Inverse Power Decay"
                      "Linear Decay"    
    param_func_choices - list of current choices (see above)
    decay - float decay exponent for the inverse power paramter function
    bvec - 2D float array specifying intercept and amplitude for parameter
           space functions 
    banded_cutoff - cutoff for banded inverse power k-space paramter function  

    basis - string specifying the basis to use - currently ony choice
            is "Fourier"
    basis_choices - list of current choices (see above)

    bumps - dictionary containing set of "bump" filters to implement
    bump_types - set of choices for "bump" filter types to add to k-space
                 paramter function; uses scipy.signal window types 
                 so consult that documentation for available types - 
                 currently only types that only require window type name
                 and size are used - current choices are: 
                 "boxcar"
                 "blackman"
                 "hann"
                 "bartlett"
                 "flattop"
                 "parzen"
                 "bohman"
                 "blackmanharris"
                 "nuttall"
                 "barthann"
    bump_default_type - the default window type used (currently "blackman")

    """
    
    def __init__(self,prior = "Gaussian",likelihood = "Gaussian",likelihood_scale = 0.05,prior_scale=0.0005,basis="Fourier"):
        """ 
       
        initializes class

        Inputs:
          prior - function type for prior
          likelihood - function type for likelihood
          likelihood_scale_orig - likelihood scale at k-space origin
          prior_scale - prior scale at k-space origin
          basis - basis type for transfrom space 

        Outputs:
          initializes class

        """

        self.init_image = 0
        self.final_image = 0
        self.bifsk_image = 0
        self.image_file_loaded = False
        self.initial_image_file_name = ''
        self.imdim = None
        self.view3Dslice = [0,0.5]
        self.prior = prior
        self.prior_choices = ["Gaussian"]
        self.prior_scale = prior_scale
        self.prior_scale_orig = 10.**7
        self.prior_mean_init = 0.0
        self.likelihood = likelihood
        self.likelihood_choices = ["Gaussian","Rician"]
        self.likelihood_scale = likelihood_scale
        self.rice_denom_cutoff = 0.0002
        
        # The following 3 parameters are used in conjuction
        # with a closed form for the posterior MAP estimate
        # with a Gaussian prior and Rician likelihood obtained
        # using sympy and based on an approximation to the first
        # order modified Bessel function discussed in:
        # A Simple and Efficient Approximation to the Modified Bessel
        # Functions and Its Applications to Rician Fading
        #     Ramy Salahat, Ehab Salahat, Ali Hakam and Tasneem Assaf
        # 2013 IEEE GCC Conference and exhibition, November 17-20, Doha, Qatar
        #
        self.bessel_approx_lims = np.array([0.0,11.5,20.0,37.5])
        self.bessel_approx_array  = np.array([[0.7536,0.4710,0.9807,1.144],
                                              [0.9739,-163.40,0.8672,0.995],
                                              [-0.715,0.9852,1.0795,0.5686],
                                              [0.2343,0.8554,1.0385,0.946]])
        self.basis_choices = ["Fourier"]
        self.basis = basis
        self.bas = None
        # Expand the following to if/else when more basis choices added
        # For now set default parameters here; make them editable via
        # pop widget
        if self.basis == "Fourier":
            from .bases import fourier as fb
            self.bas = fb
            self.decay = 0.5
            self.param_func_type = self.bas.param_func_type
            self.param_func_choices = self.bas.param_func_choices
            if self.param_func_type == "Inverse Power Decay":
                self.bvec = self.bas.bvec_ixsc
            elif self.param_func_type == "Banded Inverse Power Decay":
                self.bvec = self.bas.bvec_ixscbanded
            elif self.param_func_type == "Linear Decay":
                self.bvec = self.bas.bvec_linsc
            else:
                pass
            self.bvec[1] = 500.
            # self.bvec[1] = 546.
            self.banded_cutoff = 50.
            # Create an empty bump dictionary
            # keys are text strings for scipy.signal.window
            # available bump types from scipy.signal.window are 
            # (others like gaussian require more paramters):
            #
            # Due to our use in k-space picking the preferred versions
            # of these, which are the versions that end at 0
            # boxcar - kspace ring filter 
            # triang - pretty self explanatory - XXX - use bartlett instead 
            # blackman - nice decaying tails <- use as default
            # hamming - nice decaying tails - XXX
            # hann - nice decaying tails <- use this instead of hamming
            # bartlett - same as triang but ending at 0
            # flattop - fast decaying to below 0 and back <- should be 2nd def.
            # parzen - nice decaying tails
            # bohman - nice decaying tails
            # blackmanharris - good narrow, decaying tails 
            # nuttall - good narrow, decaying tail
            # barthann - kind of fat tails with pointy top
            #
            # Read about the properties of these at:
            # https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.get_window.html#scipy.signal.windows.get_window
            # functions and values are 3 element arrays containing:
            # 1) the position of the k-space filter in terms
            #    of the fraction of Kx max
            # 2) the amplitude of the k-space filter in terms
            #    of the fraction of parameter function max - bvec[1]
            # 3) the width of the filter to send to the window function
            #    specified in terms of a fraction of Kx max
            self.bumps = {}
            self.bump_types = ["boxcar","blackman","hann","bartlett","flattop","parzen","bohman","blackmanharris","nuttall","barthann"]
            self.bump_default_type = "blackman"
        ##### Test 1 ##############
        ##### - These are the 'ideal' settings for a Rician
        ##### likelihood analysis
        # self.rice_denom_cutoff = 0.00025
        # self.likelihood = "Rician"
        # self.likelihood_scale = 0.0005
        # self.prior_scale = 0.00009
        ##### Test 1 ###############
        
    def save_results(self):
        """

        pickles current bifs object and saves to a time stamped file
        using jsonpickle; in addition saves processed final image and 
        final transform space image to files with the same time stamp 
        in case user just wants copies of final images without having 
        to reload bifs object.

        Inputs:

          current bifs object
         
        Outputs:
          saves time stamped files containing images and bifs parameters

        """
        if np.isscalar(self.final_image):
            print("Error: There is no final image to output: probably need to run BIFS_MAP()")
            return
        else:
            # Date stamp for parameter output file
            date_stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
            # Output file names
            if self.image_file_loaded:
                # Assume here that input file is in standard 2 section,
                # period seperated form with file suffix as image type
                # First strip off any possible preceeding directory names
                in_file_name = self.initial_image_file_name
                in_file = in_file_name.split('/')[-1]
                file_prefix = in_file.split('.')[0]
                image_type = in_file.split('.')[1]
            else:
                # file_prefix = 'bifs_np_array_input'
                in_file_name = "Numpy Array"
                image_type = 'bmp'
            
            # Parameter file:
            out_p = 'bifs_params_'+date_stamp+'.bifspout'
            out_im = 'bifs_output_image_'+date_stamp+'.'+image_type
            out_k = 'bifs_output_kspace_'+date_stamp+'.'+image_type
            # Center final K-space image
            out_k_im = np.roll(np.roll(self.bifsk_image,self.bifsk_image.shape[0]//2+1,0),self.bifsk_image.shape[1]//2,1)
            # Output images
            plt.imsave(out_im,self.final_image,cmap = cm.Greys_r)
            plt.imsave(out_k,np.log(out_k_im),cmap = cm.Greys_r)
            # Output paramter file
            param_file = open(out_p, "w")
            param_file.write(jsp.encode(self))
            param_file.close()
            return

    def copy_params(self):
        """

        Return a new bifs object that has all the basic parameter values found in self.

        Does not copy the image or filename
        """
        newbifs = bifs()
        newbifs.param_func_type = self.param_func_type
        newbifs.decay = self.decay
        newbifs.prior = self.prior
        newbifs.prior_scale = self.prior_scale
        newbifs.likelihood = self.likelihood
        newbifs.likelihood_scale = self.likelihood_scale
        newbifs.bumps = self.bumps
        newbifs.view3Dslice = self.view3Dslice
        if hasattr(self, "prior_mean"):
            newbifs.prior_mean = self.prior_mean
            newbifs.prior_std = self.prior_std
        return newbifs


    def load_image_file(self,fileName):
        """
        
        load an image file using name fileName; can use to
        examine image without invoking full bifs initialization.
        
        Inputs:
          fileName - name of file to load

        Outputs:
          loads file into bifs object

        RB makes read_imfile persistent in some cases so bulk scans can get the headers.

        """
        self.initial_image_file_name = fileName
        try:
            # For now the convention is that 1D data sets
            # will have filennames ending in '.txt' and
            # will be in numpy text file form, e.g. as created
            # with numpy.savetxt with the conventions that
            # header comments start with # and there is no
            # sperator between lines
            if self.initial_image_file_name.lower().endswith('.txt'):
                # numpy.loadtext automatacilly reads to numpy array
                read_image = np.loadtxt(self.initial_image_file_name)
            else:
                try:
                    import nibabel
                    # if nibabel is not installed, move on
                    # if nibabel can't read file, move on
                    # nibabel.streamlines.is_supported(fname) would seem to be a good test, but files
                    # that fail it can still be read, e.g., nii.gz files.
                    # So we use the Python idiom "Better to ask forgiveness than permission"
                    self.read_imfile = nibabel.load(self.initial_image_file_name)
                    read_image = self.read_imfile.get_fdata()
                except:
                    try:
                        self.read_imfile = imageio.volread(self.initial_image_file_name)
                        assert len(self.read_imfile) > 2
                        read_image = np.asarray(self.read_imfile)
                    except:
                        # we have a 2D, or maybe 1D, image
                        self.read_imfile = imageio.imread(self.initial_image_file_name)
                        read_image = np.asarray(self.read_imfile)
            self.image_file_loaded = True
            self.load_image(read_image)
        except:
            print("BIFS Couldn't read image file: ",fileName)
            return
        return

    def load_image(self,init_image):
        """

        intializes the bifs object so as to be ready for
        analysis steps.

        load_image_file calls this

        Inputs:

           init_image - array generated by loading image file using
                       load_image_file()

        Outputs:
           sets up transformed (e.g. k-spce) data and paramter
           function based on dimension and size of init_image
           and default or specified paramters 


        """
        self.init_image = init_image
        self.init_image[np.isnan(init_image)] = 0.0
        self.final_image = 0
        self.bifsk_image = 0
        self.imdim = len(init_image.shape)
        
        myShape = init_image.shape
        if self.imdim == 1:
            self.kdist = self.bas.kdist1D(*myShape)
            self.k_image = self.bas.tx1(self.init_image) # Get k-space image
        elif self.imdim == 2:
            self.kdist = self.bas.kdist2D(*myShape)
            self.k_image = self.bas.tx2(self.init_image) # Get k-space image
        elif self.imdim == 3:
            self.kdist = self.bas.kdist3D(*myShape)
            self.k_image = self.bas.txn(self.init_image) # Get k-space image

        if self.basis == "Fourier": # Add other basis functions as else...
            self.mod_image = abs(self.k_image) # Get modulus image in k-space
            self.phase_image = sp.angle(self.k_image) # Get phase image in k-space
            # self.data_std =  self.likelihood_scale*self.mod_image
            self.image_exists = True

        # Set prior via orthogonal-space parameter function
        self.set_prior_from_param_func(self.param_func_type)
        # Set Rice parameter in case you use it
        # self.rice_arg = self.mod_image/self.likelihood_scale

        #### - The following doesn't seem to work, i.e. de-emphasizes -  ####
        #### - prior to too large an extent - need to investigate this - ####
        # Do normilization
        # # Since modulus and parameter functions are positive definite
        # don't need to square them (i.e. k-space integral) to get power
        # self.norm = (np.sum(self.mod_image))/(np.sum(self.prior_mean))
        # self.prior_mean = self.norm*self.prior_mean
        # self.prior_std = self.norm*self.prior_std
        # data std is regarded as constant, related to SNR
        # need to figure out best way to estimate and normalize
        ####                                                             ####
        ####                                                             ####
        
        return

    def add_bump(self,my_key,position,amplitude,width):
        """

        adds a bump filter to the self.bumps dictionalry

        Inputs:
           my_key - the name of an appropriate scipy.signal shape
           position -  fraction of kmax for location of bump filter
           amplitude - fraction of maximum of parameter function 
                       for amplitudeof of bump filter 
           width -     fraction of kmax for width of bump filter        

        Outputs:
          adds bump filter to self.bumps dictionary

        """
        key_orig = my_key.split(".")[0] 
        if my_key in self.bumps:
            for i in range(10): # No more than ten bump filters of any kind !!!!
                new_key = key_orig+"."+str(i)
                if new_key not in self.bumps:
                    break
        else:
            new_key = my_key
        self.bumps[new_key] = np.array([np.float(position),np.float(amplitude),np.float(width)])
        return
    
    def set_prior_from_param_func(self,pft):
        """

        sets up k-space parameter function based on default or 
        specified paramters

        Inputs:
           pft = paramter space function type

        Outputs:
           sets up k-space parameter function

        """
        if np.isscalar(self.init_image):
            print("Error: The Transform Space parameter function can only be set after loading an initial image")
        else:
            # test ##########
            # The following doesn't work and I don't know why
            # It seems that the prior function (in particulr the decay
            # function) should be normalized by the zero k-space value
            # of the image (i.e. the total "power") rather than an arbitrary
            # value like 500
            # self.bvec[1] = self.mod_image[0,0]
            # print("Zero k-space value:",self.mod_image[0,0])
            # self.bvec[1] = 50.0
            # test ##########
            if self.basis == "Fourier": # the usual re. elses for other tx spaces
                # Try setting bvec[1] = self.mod_image[0]
                # if self.imdim == 1:
                #    self.bvec[1] = self.mod_image[0]
                # elif self.imdim == 2:
                #    self.bvec[1] = self.mod_image[0,0]
                # elif self.imdim == 3:
                #    self.bvec[1] = self.mod_image[0,0,0]
                # else:
                #    pass
                
                if pft == "Inverse Power Decay":
                    self.prior_mean = self.bas.ixsc(self.bvec,self.kdist,self.decay)
                    self.prior_std = self.prior_scale*self.prior_mean # default for now
                                                        # i.e. strong prior
                elif pft == "Banded Inverse Power Decay":
                    self.prior_mean = self.bas.ixscbanded(self.bvec,self.kdist,self.decay,self.banded_cutoff,self.bumps)
                    self.prior_std = self.prior_scale*self.prior_mean # default for now
                                                        # i.e. strong prior
                elif pft == "Linear Decay":
                    self.prior_mean = self.bas.linsc(self.bvec,self.kdist)
                    self.prior_std = self.prior_scale*self.prior_mean # default for now
                                                        # i.e. strong prior
                elif pft == "Empirical":
                    # should already have loaded prior_mean and prior_std
                    self.prior_std *= self.prior_scale
                    # ugly: this is done elsewhere for all other types
                    # which is weird since we reset the sd for some of the other types just above
                    self.prior_std2 = self.prior_std*self.prior_std
                else:
                    print("Please specify recognized transform space parameter function, one of:")
                    for pf in self.param_func_choices:
                        print(pf)
                    return
            # Add bumps if there are any
            if self.bumps:
                # print("Adding some bumps ",self.bumps)
                self.prior_mean += self.bas.add_bumps_to_pf(self.bvec,self.kdist,self.bumps,np.int(np.max(self.kdist)))
            # "Zero" center of transform space prior and
            # set large std at origin
            if self.imdim == 1:
                self.prior_mean_init = self.prior_mean[0]
                self.prior_mean[0] = 0.
                # self.prior_mean[0] = 1.
                self.prior_std[0] = self.prior_scale_orig
            elif self.imdim == 2:
                self.prior_mean_init = self.prior_mean[0,0]
                self.prior_mean[0,0] = 0.
                # self.prior_mean[0,0] = 1.
                self.prior_std[0,0] = self.prior_scale_orig
            elif self.imdim == 3:
                self.prior_mean_init = self.prior_mean[0,0,0]
                self.prior_mean[0,0,0] = 0.
                # self.prior_mean[0,0,0] = 1.
                self.prior_std[0,0,0] = self.prior_scale_orig
            else:
                pass
        return

    def load_empirical(self, fname):
        """Load empirical prior from named file and set mode to Empirical
        Also sets prior scale to 1, since the default value is very small.
        You can and probably should make it larger via the Gaussian gui specification.
        """
        x = np.load(fname)
        self.prior_mean = x["mean"]
        self.prior_std = x["sd"]
        self.param_func_type = "Empirical"
        self.prior_scale = 1.0
     
    def bifs_map_gauss_gauss(self,pm,ps2,mf,ds2):
        """

        returns MAP estimate for posterior function using 
        gaussian prior and likelihood and uses analytic 
        conjugate function

        Inputs:
           pm - prior mean
           ps2 - square of prior std
           mf - k-space moduls (from data)
           ds2 - square of likelihood std

        Outputs:
           MAP estimate

        """
        # return((pm/ps2 + mf/ms2) / (1./ps2 + 1./ms2))
        # more efficient:
        return((pm*ds2 + mf*ps2) / (ps2 + ds2))

    # Don't think the arguments to the prior and likelihood are
    # right yet, e.g. maybe need self.prior_mean_orig as argument
    # for scale factor of Gaussian prior - need to estimate the
    # the scale factor for the Rician from the signal to noise
    # in the image... KY
    
    def bifs_map_gauss_rice(self,pm,ps,mf,ds):
        """
        returns MAP estimate for posterior function using 
        gaussian prior and rician likelihood and uses uses analytic
        function derived using sympy and the modified Bessel
        function approximation from the paper cited above
        
        Inputs:
           pm - prior mean
           ps - prior std
           mf - k-space moduls (from data)
           ds - data std
        Outputs:
           MAP estimate
        """
        # Kludge - slows things down but we use the Gaussian conjugate
        # prior calculation to estimate where the MAP value will be and
        # ergo estimate which of the exponential coefficient values to use
        # for the modified Bessel function approximation
        #
        conj = self.bifs_map_gauss_gauss(self.prior_mean,self.prior_std2,self.mod_image,self.ksd2)
        d = (2*ds**2 - 2*ps**2)
        dabs = np.abs(d)
        # Estimate which Bessel approximation coefficients
        # to use based on the Gaussian conjugate MAP value

        sn = np.zeros(d.shape)
        sn = where(np.abs(d) > 0,np.sign(d),-1.0)
        denom = sn*np.maximum(np.abs(d),self.rice_denom_cutoff)
        b = mf/ds
        result = 0.0
        for i in range(self.bsa.shape[-1]):
            if self.imdim == 1:
                ba = self.bsa[:,i]
            elif self.imdim == 2:
                ba = self.bsa[:,:,i]
            elif self.imdim == 3:
                ba = self.bsa[:,:,:,i]
            else:
                pass
            num = (-(b*ba*ds*ps**2 - mf*ds**2 + 2*mf*ps**2 - ds**2*pm + ds*np.sqrt(b**2*ba**2*ps**4 + 2*b*ba*mf*ds*ps**2 - 2*b*ba*ds*pm*ps**2 + mf**2*ds**2 - 2*mf*ds**2*pm + ds**2*pm**2 - 4*ds**2*ps**2 + 4*ps**4)))
            ### Note - the closer the scales get to each other
            # (ergo the pathalogical blow up), the MAP estimate
            # should actually get closer to the average of the
            # the prior mean and likelihood mean - that's the
            # logic for:
            # denom = np.where(dabs < self.rice_denom_cutoff,(2*num)/(mf+pm),dabs)
            ### Doesn't seem to work so back to above denom setting with
            ### cutoff
                
            result += (num/denom)/len(ba.shape)
        # Need to include the sign part here in case the setting
        # for the denominator screwed up the sign
        return np.sign(result)*result
        
    def BIFS_MAP(self):
        """

        performs MAP otpimization individually at each k-space 
        point using multiprocessing package and the above specific 
        MAP estimation functions, recombines the results, performs
        inverse transform and returns final image 

        Inputs:

        Outputs:
          Sets final_image based on performing k-space MAP estimation 

        """
        # Break up image(s) into linear array(s) for sending off
        # to multiprocessing - this is stupid and unecessarily
        # time consuming but I still need to figure out how to send
        # chunks of multivariate arrays to a multiprocessing pool
        if np.isscalar(self.init_image):
            print ("Error: Need to load an image before running MAP")
            return
        else:
            # In case prior scale was reset:
            if self.param_func_type != "Empirical":
                self.prior_std[np.where(self.prior_mean != 0.0)] = self.prior_scale*self.prior_mean[np.where(self.prior_mean != 0.0)]
                self.prior_std[np.where(self.prior_mean == 0.0)] = self.prior_scale*self.prior_mean_init
                self.prior_std2 = self.prior_std*self.prior_std

            # In case parameter function or decay value were reset:
            self.set_prior_from_param_func(self.param_func_type)

            # Square of likelihood
            self.ksd2 = self.likelihood_scale**2
            # self.data_std2 = self.data_std**2
            # Rician parameter
            # self.rice_arg = self.mod_image/self.likelihood_scale
            if self.prior == "Gaussian" and self.likelihood == "Gaussian":
                self.bifsk_image = self.bifs_map_gauss_gauss(self.prior_mean,self.prior_std2,self.mod_image,self.ksd2)
            elif self.prior == "Gaussian" and self.likelihood == "Rician":
                conj = self.bifs_map_gauss_gauss(self.prior_mean,self.prior_std2,self.mod_image,self.ksd2)
                besind = np.zeros(self.init_image.shape,dtype=int)
                besind[np.where(conj > self.bessel_approx_lims[1])] = 1
                besind[np.where(conj > self.bessel_approx_lims[2])] = 2
                besind[np.where(conj > self.bessel_approx_lims[3])] = 3
                self.bsa = self.bessel_approx_array[besind,:]
                self.bifsk_image = self.bifs_map_gauss_rice(self.prior_mean,self.prior_std,self.mod_image,self.likelihood_scale)
            else:
                pass
            # Send back to image space
            if self.basis == "Fourier": # usual add else for other tx 
                if self.imdim == 1:
                    self.final_image = np.real(self.bas.itx1(self.bifsk_image*np.exp(1j*self.phase_image)))
                elif self.imdim == 2:
                    self.final_image = np.real(self.bas.itx2(self.bifsk_image*np.exp(1j*self.phase_image)))
                elif self.imdim == 3:
                    self.final_image = np.real(self.bas.itxn(self.bifsk_image*np.exp(1j*self.phase_image)))
            return
