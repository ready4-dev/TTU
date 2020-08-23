*************************************************************************************
  *** SPSS ALGORITHM FOR AQoL-6D ADOLESCENT ( (CHE version 11 dated 14 January 2014).***********************
                                                *** REM THIS ALGORITHM (CHE version 11 dated 14 January 2014) IS AN INTERIM RELEASE
                                              ***AND MAY BE CHANGED WITHOUT NOTICE.
                                              *** REM RESEARCHERS SHOULD CHECK WITH THE AQOL GROUP AT
                                              ***MONASH UNIVERSITY FOR ANY MODIFICATION  www.aqol.com.au
                                              *******************************************************************************************
                                                ***This file analyses  adolescent data collected using the adolescent version of the AqoL-6D. The algorithm produces a utility score for the instrument based on adult weights which is corrected for age by a factor of 1.84.
                                              ***The final output contains ADULT values for each of the dimensions which are not adjusted for age, and a final instrument utility score which is adjusted for ADOLESCENT age.
                                              ***The algorithm produces a utility score for the instrument based on adult weights which is corrected for age by a factor of 1.84.

                                              *******************************************************************************************
                                                * The AQoL-6D ADOLESCENT utility scores are scaled such that the:
                                                * "AQoL-6D worst health state" according to the algorithm output is as follows for the different countries: Australia = 0.072, Fiji= 0.094, NZ= 0,053, Tonga= 0,068 (where Death = 0.00).
                                              * "AQoL-6D best health state" = 1.00

                                              *These algorithm 'worst' scores are unreliable and the lowest adolescent's score should be used as the worst.
** Scores below this are increasingly unreliable. The further away from the mean, the more likely the error and the authors have no confidence in the algorithm 'worst' scores.
**Missing Values: Note that missing data are represented by a blank or dot and are handled by imputing values within each dimension.
*Dimensions with 3-4 items will allow for 1 missing value to be imputed.
*However, if more item responses in the dimensions are missing the observations will be dropped and there will not be an instrument score for the individual.
*** (Dimension scores are not utility values as they have not been evaluated on a life-death scale)
***The dimensions are scaled on a  "Dimension Worst Health State - Dimension Best Health State" scale, where DWHS = 0.00 and DBHS = 1.00.
**********************************************************************************
* aqol# are item responses in your data
* ‘v’ is a value score, ‘dv’ is a  disvalue score – these terms apply to item and dimension scores which have been scaled on a Best-Worst scale.
* ‘u’ is an instrument utility score, ‘du’ is an instrument disutility score  - these have been scaled on a Life-Death scale.
* Missing values represented by a blank or dot.

compute Q1  = aqol1.
compute Q2  = aqol2.
compute Q3  = aqol3.
compute Q4  = aqol4.
compute Q5  = aqol5.
compute Q6  = aqol6.
compute Q7  = aqol7.
compute Q8  = aqol8.
compute Q9  = aqol9.
compute Q10 = aqol10.
compute Q11 = aqol11.
compute Q12 = aqol12.
compute Q13 = aqol13.
compute Q14 = aqol14.
compute Q15 = aqol15.
compute Q16 =aqol16.
compute Q17 = aqol17.
compute Q18 = aqol18.
compute Q19 = aqol19.
compute Q20 = aqol20.

*****************************************************************
 *********    Imputing Missing Values in Database     *********
*****************************************************************
** Independent Living - Dimension 1**

Compute ILmiss = Nmiss (Q1,Q2,Q3,Q4).
Do if ILmiss < 2.
Do repeat
  A = Q1,Q2,Q3,Q4.
If (Missing (A)) A = RND(Mean (Q1,Q2,Q3,Q4)).
End repeat.
End if.

** Relationships - Dimension 2**

Compute RELmiss = Nmiss (Q5,Q6,Q7).
Do if RELmiss < 2.
Do repeat
  A = Q5,Q6,Q7.
If (Missing (A)) A = RND(Mean (Q5,Q6,Q7)).
End repeat.
End if.

** Mental Health - Dimension 3**

Compute MENmiss = Nmiss (Q8,Q9,Q10,Q11).
Do if MENmiss < 2.
Do repeat
  A = Q8,Q9,Q10,Q11.
If (Missing (A)) A = RND(Mean (Q8,Q9,Q10,Q11)).
End repeat.
End if.

** Coping - Dimension 4**

Compute COPmiss = Nmiss (Q12,Q13,Q14).
Do if COPmiss < 2.
Do repeat
  A = Q12,Q13,Q14.
If (Missing (A)) A = RND(Mean (Q12,Q13,Q14)).
End repeat.
End if.

** Pain - Dimension 5**

Compute PAINmiss = Nmiss (Q15,Q16,Q17).
Do if PAINmiss < 2.
Do repeat
  A = Q15,Q16,Q17.
If (Missing (A)) A = RND(Mean (Q15,Q16,Q17)).
End repeat.
End if.

** Senses - Dimension 6**

Compute SENmiss = Nmiss (Q18,Q19,Q20).
Do if SENmiss < 2.
Do repeat
  A = Q18,Q19,Q20.
If (Missing (A)) A = RND(Mean (Q18,Q19,Q20)).
End repeat.
End if.
Execute.

********************************************************************************************************
*** ITEM DISVALUES ***
********************************************************************************************************

***Dimension 1. Independent living*

***1.  Household Help

if (Q1=1) dvQ1 = 0.
if (Q1=2) dvQ1=0.073.
if (Q1=6) dvQ1=0.073.# DIFF
if (Q1=3) dvQ1=0.435.
if (Q1=4) dvQ1=0.820.
if (Q1=5) dvQ1=1.

***2. Getting Around Outside

if (Q2=1) dvQ2 = 0.
if (Q2=2) dvQ2=0.033.
if (Q2=3) dvQ2=0.240.
if (Q2=4) dvQ2=0.471.
if (Q2=5) dvQ2=0.840.
if (Q2=6) dvQ2=1.

***3. Mobility

if (Q3=1) dvQ3 = 0.
if (Q3=2) dvQ3=0.041.
if (Q3=3) dvQ3=0.251.
if (Q3=4) dvQ3=0.570.
if (Q3=5) dvQ3=0.827.
if (Q3=6) dvQ3=1.

***4. Personal Care

if (Q4=1) dvQ4 = 0.
if (Q4=2) dvQ4=0.040.
if (Q4=3) dvQ4=0.297.
if (Q4=4) dvQ4=0.797.
if (Q4=5) dvQ4=1.

*****Dimension 2. Relationships***************

***5.Intimate

if (Q5=1) dvQ5 = 0.
if (Q5=2) dvQ5=0.074.
if (Q5=3) dvQ5=0.461.
if (Q5=4) dvQ5=0.841.
if (Q5=5) dvQ5=1.

***6. Family Role

if (Q6=1) dvQ6=0.
if (Q6=2) dvQ6=0.193.
if (Q6=3) dvQ6=0.759.
if (Q6=4) dvQ6=1.

***7. Community Role

if (Q7=1) dvQ7=0.
if (Q7=2) dvQ7=0.197.
if (Q7=3) dvQ7=0.648.
if (Q7=4) dvQ7=1.

***Dimension 3. Mental Health

***8. Despair

if (Q8=1) dvQ8=0.
if (Q8=2) dvQ8=0.133.
if (Q8=3) dvQ8=0.392.
if (Q8=4) dvQ8=0.838.
if (Q8=5) dvQ8=1.

***9. Worried

if (Q9=1) dvQ9=0.
if (Q9=2) dvQ9=0.142.
if (Q9=3) dvQ9=0.392.
if (Q9=4) dvQ9=0.824.
if (Q9=5) dvQ9=1.

***10. Sad

if (Q10=1) dvQ10=0.
if (Q10=2) dvQ10=0.097.
if (Q10=3) dvQ10=0.330.
if (Q10=4) dvQ10=0.784.
if (Q10=5) dvQ10=1.

***11. Calm

if (Q11=1) dvQ11=0.
if (Q11=2) dvQ11=0.064.
if (Q11=3) dvQ11=0.368.
if (Q11=4) dvQ11=0.837.
if (Q11=5) dvQ11=1.

***Dimension 4. Coping

***12. Energy

if (Q12=1) dvQ12=0.
if (Q12=2) dvQ12=0.056.
if (Q12=3) dvQ12=0.338.
if (Q12=4) dvQ12=0.722.
if (Q12=5) dvQ12=1.

***13. Control

if (Q13=1) dvQ13=0.
if (Q13=2) dvQ13=0.055.
if (Q13=3) dvQ13=0.382.
if (Q13=4) dvQ13=0.774.
if (Q13=5) dvQ13=1.

***14. Coping

if (Q14=1) dvQ14=0.
if (Q14=2) dvQ14=0.057.
if (Q14=3) dvQ14=0.423.
if (Q14=4) dvQ14=0.826.
if (Q14=5) dvQ14=1.

***Dimension 5. Pain

***15. Serious pain

if (Q15=1) dvQ15=0.
if (Q15=2) dvQ15=0.133.
if (Q15=3) dvQ15=0.642.
if (Q15=4) dvQ15=1.

***16. Pain

if (Q16=1) dvQ16=0.
if (Q16=2) dvQ16=0.200.
if (Q16=3) dvQ16=0.758.
if (Q16=4) dvQ16=1.

***17. Pain interferes

if (Q17=1) dvQ17=0.
if (Q17=2) dvQ17=0.072.
if (Q17=3) dvQ17=0.338.
if (Q17=4) dvQ17=0.752.
if (Q17=5) dvQ17=1.

***Dimension 6. Senses

***18. Vision

if (Q18=1) dvQ18=0.
if (Q18=2) dvQ18=0.033.
if (Q18=3) dvQ18=0.223.
if (Q18=4) dvQ18=0.622.#0.521
if (Q18=5) dvQ18=0.843.
if (Q18=6) dvQ18=1.

***19. Hearing

if (Q19=1) dvQ19=0.
if (Q19=2) dvQ19=0.024.
if (Q19=3) dvQ19=0.205.
if (Q19=4) dvQ19=0.586.
if (Q19=5) dvQ19=0.826.
if (Q19=6) dvQ19=1.

***20. Communicate

if (Q20=1) dvQ20=0.
if (Q20=2) dvQ20=0.187.
if (Q20=3) dvQ20=0.695.
if (Q20=4) dvQ20=1.

*********************************************************************************************************
***2. DIMENSION SCORES ***
***********************************************************************************************The dimension scores in the Adolescent version of the AQoL-6D are based on adult weights and these scores are not econometrically transformed at the end as the instrument utility score is. They are included in the output.

*** (Dimension scores are not strict utility values as they have not been evaluated on a life-death scale)

*** The dimensions are scaled on a "Dimension Worst Health State - Dimension Best Health State" scale where DWHS = 0.00 and DBHS = 1.00.***

***DIMENSION 1 - IND LIV.

***DIMENSION SCALING CONSTANT kD1=-0.978.
***IND LIV HAS 4 ITEMS.

***ITEM WORST WEIGHTS (Wi).
**This model uses w1=0.385412.
**This model uses w2=0.593819.
**This model uses w3=0.630323.
**This model uses w4=0.794888.
**4 item formula.
**dvD1=(1/kD1)*[(1+(kD1*w1*dvQ1))*(1+(kD1*w2*dvQ2))*(1+(kD1*w3*dvQ3))*(1+(kD1*w4*dvQ4))-1].

Compute dvD1=(1/-0.978)*((1+(-0.978*0.385412*dvQ1))*(1+(-0.978*0.593819*dvQ2))*(1+(-0.978*0.630323*dvQ3))*(1+(-0.978*0.794888*dvQ4))-1).

**Variable dvD1 = "Disvalue Score for Dimension 1 - Independent Living".

***DIMENSION 2 - REL.

***DIMENSION SCALING CONSTANT kD2 = -0.923.
***REL HAS 3 ITEMS.

***ITEM WORST WEIGHTS (Wi).
**This model uses w5=0.64303.
**This model uses w6=0.697742.
**This model uses w7=0.508658.

**3 item formula
**dvD2=(1/kD2)*[(1+(kD2*w5*dvQ5))*(1+(kD2*w6*dvQ6))*(1+(kD2*w7*dvQ7))-1].

Compute dvD2=(1/-0.923)*((1+(-0.923*0.64303*dvQ5))*(1+(-0.923*0.697742*dvQ6))*(1+(-0.923*0.508658*dvQ7))-1).

**Variable dvD2 = "Disvalue Score for Dimension 2 - Relationships". ***

***DIMENSION 3 - MEN. ***

***DIMENSION SCALING CONSTANT kD3 = -0.983.  ***
***MEN HAS 4 ITEMS.

***ITEM WORST WEIGHTS (Wi).  ***
**This model uses w8=0.640377.	 ***
**This model uses w9=0.588422.	 ***
**This model uses w10=0.648748.	 ***
**This model uses w11=0.71122. ***

Compute dvD3=(1/-0.983)*((1+(-0.983*0.640377*dvQ8))*(1+(-0.983*0.588422*dvQ9))*(1+(-0.983*0.648748*dvQ10))*(1+(-0.983*0.71122*dvQ11))-1).

**Variable dvD3 = "Disvalue Score for Dimension 3 - Mental Health".  ***

***DIMENSION 4 - COPING.

***DIMENSION SCALING CONSTANT kD4 = -0.930.  ***
***COPI HAS 3 ITEMS.  ***

***ITEM WORST WEIGHTS (Wi).  ***
**This model uses w12=0.415694.	 ***
**This model uses w13=0.636994.	 ***
**This model uses w14=0.773296.	 ***

Compute dvD4=(1/-0.930)*((1+(-0.930*0.415694*dvQ12))*(1+(-0.930*0.636994*dvQ13))*(1+(-0.930*0.773296*dvQ14))-1).

**Variable dvD4 = "Disvalue Score for Dimension 4 - Coping".

***DIMENSION 5 - PAIN.

***DIMENSION SCALING CONSTANT kD5 = -0.962. #-0.96
***PAIN HAS 3 ITEMS.
***ITEM WORST WEIGHTS (Wi).
**	w15=0.631833.
**	w16=0.767573.
**	w17=0.652241.

Compute dvD5=(1/-0.962)*((1+(-0.962*0.631833*dvQ15))*(1+(-0.962*0.767573*dvQ16))*(1+(-0.962*0.652241*dvQ17))-1).

**Variable dvD5 = "Disvalue Score for Dimension 5 - Pain".

***DIMENSION 6 - SEN.
***DIMENSION SCALING CONSTANT kD6 = -0.851.
***SEN HAS 3 ITEMS.
***ITEM WORST WEIGHTS (Wi).
**	w18=0.580696.
**	w19=0.463022.
**	w20=0.604613.

Compute dvD6=(1/-0.851)*((1+(-0.851*0.580696*dvQ18))*(1+(-0.851*0.463022*dvQ19))*(1+(-0.851*0.604613*dvQ20))-1).

**Variable dvD6 = "Disvalue Score for Dimension 6 - Senses".

Compute vD1 =1-dvD1.
Compute vD2 =1-dvD2.
Compute vD3 =1-dvD3.
Compute vD4 =1-dvD4.
Compute vD5 =1-dvD5.
Compute vD6 =1-dvD6.
Execute.

VARIABLE LABELS vD1  "Adult Score Dimension 1 - Independent Living".
VARIABLE LABELS vD2  "Adult Score Dimension 2 - Relationships".
VARIABLE LABELS vD3  "Adult Score Dimension 3 - Mental Health".
VARIABLE LABELS vD4  "Adult Score Dimension 4 - Coping".
VARIABLE LABELS vD5  "Adult Score Dimension 5 - Pain".
VARIABLE LABELS vD6  "Adult Score Dimension 6 - Senses".
EXECUTE.

***3. OVERALL SCORE ON A 0-1 DISVALUE SCALE***
*** On a scale, Utility/Disutility involves preferences, Value/Disvalue does not involve preferences***
**********************************************************************************************************
***DIMENSION SCALING CONSTANT kA = -0.965 # NEW
***DIMENSION WORST WEIGHTS (wDi) # NEW
* wD1=0.472 # NEW
* wD2=0.448 # NEW
* wD3=0.479 # NEW
* wD4=0.345 # NEW
* wD5=0.592 # NEW
* wD6=0.637 # NEW

*scaling factor for wDi = 1/1.132181 = 0.883 # NEW

**6 dimension formula # NEW
**duaqol = 1/kA)*[(1+(kA*wD1*0.883*dvD1))*(1+(kA*wD2*0.883*dvD2))*(1+(kA*wD3*0.883*dvD3))*(1+(kA*wD4*0.883*dvD4))*(1+(kA*wD5*0.883*dvD5))* (1+(kA*wD6*0.883*dvD6))-1]

Compute duaqol =(1/-0.965)*((1+(-0.965*0.472*0.883*dvD1))*(1+(-0.965*0.448*0.883*dvD2))*
(1+(-0.965*0.479*0.883*dvD3))*(1+(-0.965*0.345*0.883*dvD4))*(1+(-0.965*0.592*0.883*dvD5))*
(1+(-0.965*0.637*0.883*dvD6))-1).


*** 4. OVERALL INSTRUMENT SCORE ON A LIFE-DEATH DISUTILITY SCALE***
*************************************************************************************************

**Scaling constant for LIFE – DEATH scale, wld = 1.132
** duaqolld = wld*duaqol

Compute duaqolld =1.132*duaqol. # NEW


**************************************************************************************************
***AQoL-6D MULTIPLICATIVE MODE- ECONOMETRIC CORRECTION FOR ADOLESCENT AGE ***
**************************************************************************************************

Compute duAQoL6D = duaqolld**1.8407651. # NEW

If duAQoL6D > 1  duAQoL6D = (1 + (duAQoL6D-1)*.25773196). # NEW


**Variable duAQoL6D = " AQoL-6D AdolescentDisutility Score".  ***

***INSTRUMENT UTILITY SCORE
****************************************************

Compute uAQoL6D = (1- duAQoL6D). # NEW


***Compute uAQoL6DRotated = uAQoL6D*.8807+.0889
***produces best .97 worst .03
***hence add .03 to constant = .0889 + .03 = .1189

Compute uAQoL6DRotated = uAQoL6D*.8807+.1189. # NEW

*********************************************************************************************************************************************
*** ECONOMETRIC CORRECTION
*********************************************************************************************************************************************

Compute uaqol=1-(1-uAQoL6DRotated)**1.19. # NEW
EXECUTE.

VARIABLE LABELS uaqol 'AQOL 6D ADOLESCENT Utility Score'.
EXECUTE.

DESCRIPTIVES VARIABLES=vD1 vD2 vD3 vD4 vD5 vD6 uaqol
     /STATISTICS=MEAN STDDEV RANGE MIN MAX SEMEAN.
EXECUTE.

Delete Variables Q1 Q2 Q3 Q4 Q5 Q6 Q7 Q8 Q9 Q10 Q11 Q12 Q13 Q14 Q15 Q16 Q17 Q18 Q19 Q20
ILmiss  RELmiss MENmiss COPmiss PAINmiss SENmiss
dvQ1 dvQ2 dvQ3 dvQ4 dvQ5 dvQ6 dvQ7 dvQ8 dvQ9 dvQ10 dvQ11 dvQ12 dvQ13 dvQ14 dvQ15 dvQ16 dvQ17 dvQ18 dvQ19 dvQ20
duaqol duaqolld dvD1 dvD2 dvD3 dvD4 dvD5 dvD6 duAQoL6D uAQoL6D uAQoL6DRotated.
Execute.
