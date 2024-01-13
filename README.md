# R-package-trajMSM
The package trajMSM is based on the paper Marginal Structural Models with Latent Class Growth
Analysis of Treatment Trajectories: https://doi.org/10.48550/arXiv.2105.12720. Latent class growth
analysis (LCGA) are increasingly proposed as a solution to summarize the observed longitudinal
treatment in a few distinct groups. When combined with standard approaches like Cox proportional
hazards models, LCGA can fail to control time-dependent confounding bias because of timevarying
covariates that have a double role of confounders and mediators. We propose to use LCGA
to classify individuals into a few latent classes based on their medication adherence pattern, then
choose a working marginal structural model (MSM) that relates the outcome to these groups. The
parameter of interest is nonparametrically defined as the projection of the true MSM onto the chosen
working model. The combination of LCGA with MSM (LCGA-MSM) is a convenient way
to describe treatment adherence and can effectively control time-dependent confounding. Several
approaches exist to estimate the parameters of a MSM and one of the most popular is the inverse
probability weighting (IPW). The IPW mimics a random assignment of the treatment by creating
a pseudo-population where the treated and the untreated groups are comparable. In longitudinal
settings, IPW can appropriately adjust for time-varying covariates affected by prior exposure and
selection bias. In this first version, we proposed to estimate parameters of the LCGA-MSM using
the IPW, g-computation and pooled LTMLE. We proposed an extension of the LCGA-MSM to a time-dependent outcome.
We called this approach LCGA-HRMSM for LCGA and history-rectricted HRMSM. The same three estimators: IPW, g-computation
and pooled LTMLE are proposed to estimate parameters of LCGA-HRMSMs. 

To access the R codes: https://github.com/awamaeva/R-package-trajMSM.
