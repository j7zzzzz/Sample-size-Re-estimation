# README #
**1. all_function.R**: Contains all the functions used in this paper.

**2. sample_size.R**: Contains the code for the sample size re-estimation formula proposed in this paper, which is used to calculate the required number of patients for the trial.

**3. power_calculation.R**: Used in the simulation studies and redesign clinical trials sections to validate the performance of our proposed sample size estimation method in trials without interim analysis.

**4. power_calculation_interim.R**: Used in the simulation studies and redesign clinical trials sections to validate the performance of our proposed sample size re-estimation method in trials with interim analysis.

---
**5. Additional Information**:

In power_calculation.R and power_calculation_interim.R, the different cases correspond to the following sections in the paper:



**SIMULATION STUDIES**

(1) case1: Corresponds to case 1 in Section 5.1.<br>
(2) case2: Corresponds to case 2 in Section 5.1.<br>
(3) case3_K=3: Corresponds to case 3 in Section 5.1, with $Z_1$ discretized into three groups.<br>
(4) case3_K=4: Corresponds to case 3 in Section 5.1, with $Z_1$ discretized into four groups.<br>
(5) case4_K=3: Corresponds to case 4 in Section 5.1, with $Z_1$ following a three-category multinomial distribution.<br>
(6) case4_K=4: Corresponds to case 4 in Section 5.1, with $Z_1$ following a four-category multinomial distribution.<br>
(7) case5_K=3: Corresponds to case 5 in Section 5.1, with $Z_1$ following a three-category multinomial distribution.<br>
(8) case5_K=4: Corresponds to case 5 in Section 5.1, with $Z_1$ following a four-category multinomial distribution.<br>


**REDESIGN CLINICAL TRIALS**

(1) caseR01: Uses the parameters of the RE01 trial from Section 6.1, with a design similar to case 1 in Section 5.1.<br>
(2) case3_R01: Uses the parameters of the RE01 trial from Section 6.1, with a design similar to case 3 in Section 5.1, and $Z_1$ discretized into three groups.<br>
(3) case4_R01: Uses the parameters of the RE01 trial from Section 6.1, with a design similar to case 4 in Section 5.1, and $Z_1$ following a three-category multinomial distribution.<br>
(4) case5_R01: Uses the parameters of the RE01 trial from Section 6.1, with a design similar to case 5 in Section 5.1, and $Z_1$ following a three-category multinomial distribution.<br>

(5) caseL: Uses the parameters of the Lung Cancer trial from Section 6.2, with a design similar to case 2 in Section 5.1.<br>
(6) case3_L: Uses the parameters of the Lung Cancer trial from Section 6.2, with a design similar to case 3 in Section 5.1, and $Z_1$ discretized into four groups.<br>
(7) case4_L: Uses the parameters of the Lung Cancer trial from Section 6.2, with a design similar to case 4 in Section 5.1, and $Z_1$ following a four-category multinomial distribution.<br>
(8) case5_L: Uses the parameters of the Lung Cancer trial from Section 6.2, with a design similar to case 5 in Section 5.1, and $Z_1$ following a four-category multinomial distribution.
