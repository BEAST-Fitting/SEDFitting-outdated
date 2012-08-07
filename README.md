# Stellar Nut Cracker Code

Bayesian stellar SED fitting for PHAT.

### Author:

Karl Gordon (STScI)

### Acknowledgements:

Help on the stellar models for this code has come from:
- Danny Lennon (STScI)

Help on Bayesian techniques for this code has come from:
- David Hogg (NYU)
- Dan Weisz, Morgan Fouesneau (UofWa)

Support for this work from the Panchromatic Hubble Andromeda Treasury (PHAT)
- Julianne Dalcanton (PI, UofWa)

### License:

Copyright 2011, 2012 the author.  All rights reserved.

If you have interest in using or re-using any of this content, get in touch with Karl Gordon.

### Revision history:

7 Aug 2012: Basic testing and bugfixing done.  Code runs and produces
            reasonable output.  More extensive testing still needed.
            Enhancements:
              1) Major changes to the output files
                 output now down in sets of stars (chunks) to save creating a new file for every stars
                 full likelihood output done as a sparse matrix - [indices of 4D matric,log likelihood values]
                 output now much smaller
              2) using log likelihoods instead of linear likelihoods
              3) major update of variable names to be correct and clearer
              4) streamlining of code to cache when possible to improve speed
              5) output catalog has most probable values of model parameters instead of best fit (max likelihood)

6 Aug 2012: Finished cleanup and enhancement.  Now to test the code.  
            Commiting to make sure nothing gets lost.

2  Aug 2012: First batch of cleaning and documentation
             Not finished yet and probably is currently broken
             some changes not backward compatible (testing needed)

25 Jul 2012: Commit of current code (needs cleaning and documentation)
             Prior work done outside of github.

