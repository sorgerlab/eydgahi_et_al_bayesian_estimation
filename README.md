Eydgahi et al. 2013 supplementary materials
===========================================

Code and data to support Eydgahi, H. et al. Properties of cell death models
calibrated and compared using Bayesian approaches. Mol Syst Biol 9, (2013).
doi:10.1038/msb.2012.69

The code supplied will perform one MCMC random walk for Bayesian
estimation of the parameters of a previously published ODE-based model 
of receptor-mediated apoptosis (EARM 1.3). (Albeck et al, 2008).

The master .m file to run the code is: BAYESIAN_OPTIONS.M
To run it, you will need to provide it values for the following variables:

TIME_MAX: this is an integer; it specifices the number of steps in the MCMC
algorithm.

NUM_NODES: this is an integer; it is the number of nodes on the cluster
that should be used to calculate the hessian matrix.
Typically a value of 70 is chosen.

HESS_OPTION: this is a string; it takes two values: 'hess' or 'nohess'
'hess' indicates that a hessian-guided random walk should be performed.
'nohess' indicates that a non-guided random walk should be performed.

START_OPTION: this is a string; it takes three values: 'random', 'same',
or 'continue'
'random' indicates that the first position in parameter space is chosen
randomly.
'same' indicates that the same position is used to start the MCMC random
walk; this position is set to be the nominal values for the model as
reported in Albeck et al. (2008).
'continue' is used if the cluster for some reason crashes. since the
algorithm saves everything every 10000 steps, choosing continue will
revert you to the last saved position.
Typically 'random' or 'continue' is chosen.

BOUNDARY_OPTION: this is a string; it takes two values: 'boundary' or
'noboundary'
'boundary' is used to restrict the algorithm to parameter values +/- 2
around the mean values reported by Albeck et al (2008).
'noboundary' is enables the algorithm to sample all regions in parameter
space.
Typically 'noboundary' is chosen.

PRIOR_OPTION: this is a string; it takes in two values: 'albeck' or
'general'. it essentially determines the values for the mean and variance
that should be used to construct the prior.
'albeck' uses the mean parameter values as reported in Albeck et al
(2008) and a variance of 2 in log space on all parameters.
'general' uses the mean and variance values obtained from a literature
search and as reported in Table 1 in the manuscript.

TRAJ_OPTION: this is an integer; it takes in either a value or 1 or 2.
1 indicates that the algorithm should be calibrated to the EC-RP
trajectory only.
2 indicates that the algorithm should be calibrated to both EC-RP and
IC-RP trajectories.

DATA_FILE: this is a call to function that has the biological data.
If the value for TRAJ_OPTION is 1, we use either @albeck_data_2
@sabrina_slope2
if the value for TRAJ_OPTION is 2, nothing needs to be inserted

START_FILE: if the option for START_OPTION is set to 'continue', then a
value for START_FILE must be provided. this is a string that has the name
of the .mat file that you would like the algorithm to load in order to
complete the MCMC algorithm.

An example of how to start the algorithm is as follows:
bayesian_options(1000000,70,'hess','random','noboundary','general',1,@albeck_data_2)





This section of the README file describes the purpose of the .m files
necessary to run the algorithm. It is important that all of these files
must be downloaded into the same directory in order for the algorithm to
properly run.

ACTUAL_INIT_DATA.M: this is an .m file that has the parameter values for
k_f, k_r, and k_c (forward, reverse, and catalysis, respectively) as
reported in Albeck et al. (2008).

ALBECK_DATA_2.M: this is an .m file that contains data for calibrating to
one signal (EC-RP) as reported in Albeck et al. (2009)

BAYESIAN_OPTIONS.M: this is the master .m file that must be run.

HESS_GRADIENT_CALC.M: it uses a central-difference approximation to
calculate the value of hessian for the elements of the matrix as provided
by HESSIAN_PARALLEL_SHORT.M

HESS_PARALLEL_NODES.M: based on the values for NUM_NODES, it creates an
appropriate number of jobs to run on the cluster.

HESSIAN_PARALLEL_SHORT.M: it determines which elements of the hessian
matrix should be calculated by a particular node on the cluster.

JACOBIAN.M: it solves the system of 69 ODEs that describes the model.

LIKELIHOOD_GEN.M: calculates the -ln(likelihood) value for a given
position in parameter space.

ODES.M: it provides the 69 ODEs of the model.

PRIOR_GEN.M: calculates the -ln(prior) value for a given position in
parameter space.

PRIOR_STATS.M: it determines the values for the mean and variance
that should be used to construct the prior; it uses the mean parameter
values as reported in Albeck et al (2008) and a variance of 2 in log 
space for all parameters.

PRIOR_STATS_OLD.M: it determines the values for the mean and variance
that should be used to construct the prior; it uses the mean and variance
values obtained from a literature. These values are reported in Table 1 
of the manuscript

RAND_GAUSS.M: it generates a new step by taking a step of size 
NORM_STEP_SIZE (set in BAYESIAN_OPTIONS to 0.75) in log space (this 
function is used for non-guided portions of the random walk).

SABRINA_SLOPE_EC.M: this is an .m file that is used when calibrating to 
two trajectories(EC-RP and IC-RP); it provides the EC-RP data.

SABRINA_SLOPE_IC.M: this is an .m file that is used when calibrating to 
two trajectories(EC-RP and IC-RP); it provides the IC-RP data.

SABRINA_SLOPE2.M: this is an .m file that contains data for calibrating to
one signal (EC-RP) as reported in Spencer et al. (2009)

STEP_GENERATOR.M: it uses the eigenvalues and eigenvectors of the hessian
matrix to generate is new position in parameter space (this function is 
used for hessian-guided portions of the random walk).

VALUES.M: it provides the values for k_deg, k_c, k_f, and k_r
(degradation, catalysis, forward, and reverse rates, respectively) as
reported in Albeck et al (2008).

VALUES_CONC.M: it provides intial protein concetrations and k_s
(synthesis rate) as reported in Albeck et al (2008).
