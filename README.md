# magnetostriction-spline
The attached MATLAB scripts produce the results presented in the manuscript “Flexible Identification Procedure for Thermodynamic Constitutive Models for Magnetostrictive Materials”. They have been tested with MATLAB R2017b. The scripts can be directly run from the MATLAB Editor or Command Window. The content of the scripts is briefly described below. 

- Script1_Produce_SMS_Results.m
  
  The script produces the reference data with the simplified multiscale model, against which the splines are fitted. It uses the methods described in Section 2.5.

- Script2_Fit2_Hsig_SMS.m

  The script fits the bivariate spline against the data produced with the simplified multiscale (SMS) model using the magnetic field strength **_H_** and stress **_σ_** as the state variables. It uses the methods described in Section 2.2 of the paper and produces the results presented in Section 3.1 and Figure 4.

- Script3_Fit3_Hsig_SMS.m

  The script fits the trivariate spline against the data produced with the SMS model using **_H_** and **_σ_** as the state variables. It uses the methods described in Section 2.2 of the paper and produces the results presented in Section 3.1 and Figure 5.

- Script4_Fit2_Beps_SMS.m

  The script fits the bivariate spline against the data produced with the SMS model using the magnetic flux density **_B_** and strain **_ε_** as the state variables. It uses the methods described in Section 2.3 and produces the results presented in Section 3.1 and Figure 6.

- Script5_Fit3_Beps_SMS.m

  The script fits the trivariate spline against the data produced with the SMS model using **_B_** and **_ε_** as the state variables. It uses the methods described in Section 2.3 and produces the results presented in Section 3.1 and Figure 7.

- Script6_Comparison1_Hsig.m

  The script compares the spline-based thermodynamic model and the SMS model under different multiaxial stresses. It uses the methods described in Sections 2.5 and 3.3 and produces the results presented in Section 3.3 and Figure 9.

- Script7_Comparison2_Hsig.m

  The script compares the spline-based thermodynamic model and the SMS model under shear stresses. It uses the methods described in Sections 2.5 and 3.4 and produces the results presented in Section 3.4 and Figure 10.

- Script8_Fit2_Beps_Measurement.m

  The script fits the bivariate spline against measurement data using **_B_** and **_ε_** as the state variables. It uses the methods described in Section 2.3 and produces the results presented in Section 3.5 and Figure 11.

- Script9_Fit2_Hsig_SMS_error.m

  The script fits the bivariate spline against the data produced with the SMS model using the magnetic field strength **_H_** and stress **_σ_** as the state variables similarly to Script2_Fit2_Hsig_SMS.m., but adding some artificial numerical error to the SMS results. It produces the results presented in Section 3.2 and Figure 8 (a).

- Script10_Fit3_Hsig_SMS_error.m

  The script fits the trivariate spline against the data produced with the SMS model using **_H_** and **_σ_** as the state variables similarly to Script3_Fit3_Hsig_SMS.m., but adding some artificial numerical error to the SMS results. It produces the results presented in Section 3.2 and Figure 8 (b).
  