Lab Meeting
========================================================
author: Meridith L. Bartley 



Data Summary
========================================================
- two hours of observations
- high density (1 chamber) & low density (4 chambers)
- 60 ants interacting 
- 232 interactions in high density
  - 10 with queen
- 288 interactions in low density
  - 17 with queen

Interactions Per Ant
========================================================
$\\$

![plot of chunk unnamed-chunk-1](LabMeeting-figure/unnamed-chunk-1-1.png)
- 60 individuals in each

***

$\\$

![plot of chunk unnamed-chunk-2](LabMeeting-figure/unnamed-chunk-2-1.png)

Interactions per Ant by Location
========================================================
$\\$

![plot of chunk unnamed-chunk-3](LabMeeting-figure/unnamed-chunk-3-1.png)
- 37 and 32 individuals

***
$\\$

![plot of chunk unnamed-chunk-4](LabMeeting-figure/unnamed-chunk-4-1.png)



Total Interactions over Time
========================================================
$\\$

![plot of chunk unnamed-chunk-5](LabMeeting-figure/unnamed-chunk-5-1.png)

***

$\\$

![plot of chunk unnamed-chunk-6](LabMeeting-figure/unnamed-chunk-6-1.png)

Low Density Interactions by Location
========================================================
$\\$


<img src="LabMeeting-figure/unnamed-chunk-7-1.png" title="plot of chunk unnamed-chunk-7" alt="plot of chunk unnamed-chunk-7" style="display: block; margin: auto;" />

***
$\\$

<img src="LabMeeting-figure/unnamed-chunk-8-1.png" title="plot of chunk unnamed-chunk-8" alt="plot of chunk unnamed-chunk-8" style="display: block; margin: auto;" />


Model for Interactions
========================================================

Observed Interactions: $\text{ N}_t \sim \text{Pois}(\lambda_{X_t})$

  Unobserved States: $\text{ X}_t | \text{ X}_{t-1} \sim \text{Multinom} (1, \underline{p}_{X_{t-1}})$

  where, 
  
  $$\mathbf{P} = \begin{bmatrix} p_{ 11 } & \dots  & p_{ 1n } \\ \vdots  & \ddots  & \vdots  \\ p_{ n1 } & \dots  & p_{ nn } \end{bmatrix}$$
    
$\lambda_{X_t} \sim \text{Gamma}(\alpha, \beta)$

$\underline{p}_{X_t} \sim \text{Dirichlet}(\underline{\theta})$

    
    
MCMC - High Density
========================================================

![plot of chunk unnamed-chunk-9](LabMeeting-figure/unnamed-chunk-9-1.png)

2 State Model - High Density
=========================================================

<img src="LabMeeting-figure/unnamed-chunk-10-1.png" title="plot of chunk unnamed-chunk-10" alt="plot of chunk unnamed-chunk-10" style="display: block; margin: auto;" />


MCMC - Low Density
========================================================

![plot of chunk unnamed-chunk-11](LabMeeting-figure/unnamed-chunk-11-1.png)

2 State Model - Low Density
=========================================================

<img src="LabMeeting-figure/unnamed-chunk-12-1.png" title="plot of chunk unnamed-chunk-12" alt="plot of chunk unnamed-chunk-12" style="display: block; margin: auto;" />

MCMC - Low Density Location 1
========================================================

![plot of chunk unnamed-chunk-13](LabMeeting-figure/unnamed-chunk-13-1.png)

2 State Model - Low Density Location 1
=========================================================

<img src="LabMeeting-figure/unnamed-chunk-14-1.png" title="plot of chunk unnamed-chunk-14" alt="plot of chunk unnamed-chunk-14" style="display: block; margin: auto;" />

MCMC - Low Density Location 4
========================================================

![plot of chunk unnamed-chunk-15](LabMeeting-figure/unnamed-chunk-15-1.png)

2 State Model - Low Density Location 4
=========================================================

<img src="LabMeeting-figure/unnamed-chunk-16-1.png" title="plot of chunk unnamed-chunk-16" alt="plot of chunk unnamed-chunk-16" style="display: block; margin: auto;" />
