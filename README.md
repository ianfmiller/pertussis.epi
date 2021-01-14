# pertussis.epi
## Questions: How do the effects of acellular pertussis vaccines compare to those of whole-cell pertussis vaccines? Do these differences explain the resurgence of pertussis incidence? Are inferences consistent across regions and scales of inference?
## Approach: 
#### 1) Code biweekly pertussis incidence from CDC data
#### 2) Fit POMP models w/ different vaccine effects to 
##### a) identify differences/similarities between vaccine effects
##### b) identify which set of differences/similarities best explains resurgence
#### 3) Replicate model fitting across scales of inference (US, CDC reporting region, state) and regions (CDC reporting regions, states)
## Modeling Approach
### Because questions deal with differences between pertussis vaccines, I do not seek to identify similarities/differences between the effects of vaccinal and natural immunity
### Split effects of immunity into
#### 1) waning protection--rate at which protection is lost
#### 2) failure in protection--rate at which protection fails to prevent against infection
#### 3) protection against symptoms--rate at which immune individuals display symptoms if they are infected
#### 4) reduction in asymptomatic transmission--rate at which asymptomatic immune individuals transmit relative to symptomatic individuals
### Assume:
#### 1) Any infection (asymptomatic or symptomatic) results in natural immunity
#### 2) All symptomatic individuals are equally infectious
#### 3) Latent period, recovery time identical for all infected classes