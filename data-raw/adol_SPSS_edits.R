# This file analyses  adolescent data collected using the adolescent version of the AqoL-6D. The algorithm produces a utility score for the instrument based on adult weights which is corrected for age by a factor of 1.84.
# The final output contains ADULT values for each of the dimensions which are not adjusted for age, and a final instrument utility score which is adjusted for ADOLESCENT age.
# The algorithm produces a utility score for the instrument based on adult weights which is corrected for age by a factor of 1.84.
# The AQoL-6D ADOLESCENT utility scores are scaled such that the:
# "AQoL-6D worst health state" according to the algorithm output is as follows for the different countries:
# Australia = 0.072, Fiji= 0.094, NZ= 0,053, Tonga= 0,068 (where Death = 0.00).
# "AQoL-6D best health state" = 1.00
# These algorithm 'worst' scores are unreliable and the lowest adolescent's score should be used as the worst.
# Scores below this are increasingly unreliable. The further away from the mean, the more likely the error and the authors have no confidence in the algorithm 'worst' scores.
# Dimensions with 3-4 items will allow for 1 missing value to be imputed.
# However, if more item responses in the dimensions are missing the observations will be dropped and there will not be an instrument score for the individual.
# ‘v’ is a value score, ‘dv’ is a  disvalue score – these terms apply to item and dimension scores which have been scaled on a Best-Worst scale.
# ‘u’ is an instrument utility score, ‘du’ is an instrument disutility score  - these have been scaled on a Life-Death scale.

### CHANGE 1: Missing Values
## Verify impute algorithm
## Drop if missing > 2 for each dimension.

### CHANGE 2: Item Disvalues
# D1IT1 if (Q1=6) dvQ1=0.073.# DIFF
# D1IT18 if (Q18=4) dvQ18=0.622.#0.621

### CHANGE 3: Dimension scores
# DIMENSION 5 - PAIN. DIMENSION SCALING CONSTANT kD5 = -0.962. #-0.96

### CHANGE 4: OVERALL SCORE ON A 0-1 DISVALUE SCALE***
# DIMENSION SCALING CONSTANT kA = -0.965 # NEW
# DIMENSION WORST WEIGHTS (wDi) # NEW
# wD1=0.472 # NEW
# wD2=0.448 # NEW
# wD3=0.479 # NEW
# wD4=0.345 # NEW
# wD5=0.592 # NEW
# wD6=0.637 # NEW
# scaling factor for wDi = 1/1.132181 = 0.883 # NEW

# dimension formula # NEW
#duaqol = 1/kA)*[(1+(kA*wD1*0.883*dvD1))*(1+(kA*wD2*0.883*dvD2))*(1+(kA*wD3*0.883*dvD3))*(1+(kA*wD4*0.883*dvD4))*(1+(kA*wD5*0.883*dvD5))* (1+(kA*wD6*0.883*dvD6))-1]
#Compute duaqol =(1/-0.965)*((1+(-0.965*0.472*0.883*dvD1))*(1+(-0.965*0.448*0.883*dvD2))*
#(1+(-0.965*0.479*0.883*dvD3))*(1+(-0.965*0.345*0.883*dvD4))*(1+(-0.965*0.592*0.883*dvD5))*
#(1+(-0.965*0.637*0.883*dvD6))-1).

### CHANGE 5:  4. OVERALL INSTRUMENT SCORE ON A LIFE-DEATH DISUTILITY SCALE***
#Scaling constant for LIFE – DEATH scale, wld = 1.132
# duaqolld = wld*duaqol
# Compute duaqolld =1.132*duaqol. # NEW

### CHANGE 5: AQoL-6D MULTIPLICATIVE MODE- ECONOMETRIC CORRECTION FOR ADOLESCENT AGE ***
# Compute duAQoL6D = duaqolld**1.8407651. # NEW
# If duAQoL6D > 1  duAQoL6D = (1 + (duAQoL6D-1)*.25773196). # NEW
# Variable duAQoL6D = " AQoL-6D AdolescentDisutility Score".  ***

### CHANGE 6: INSTRUMENT UTILITY SCORE
# Verify:
# Compute uAQoL6D = (1- duAQoL6D). # NEW
# Add:
# Compute uAQoL6DRotated = uAQoL6D*.8807+.0889
# produces best .97 worst .03
# hence add .03 to constant = .0889 + .03 = .1189
# Compute uAQoL6DRotated = uAQoL6D*.8807+.1189. # NEW

### CHANGE 7:  ECONOMETRIC CORRECTION
# Compute uaqol=1-(1-uAQoL6DRotated)**1.19. # NEW
# VARIABLE LABELS uaqol 'AQOL 6D ADOLESCENT Utility Score'.
