# LongitudinalDynamicsOfHarrier
A study of longitudinal dynamics of AV-8A Harrier Aircraft 

Keep both files on the same working directory:
1. f_CLTFM
2. Harrier

# Contents of the study

1. Modelling of plant which includes Actuator and Engine dynamics of the motors.
2. Model of the plant - Ap, Bp, Cp, Dp, After including actuator there are total 6 states and 2 inputs and 2 outputs.
3. Modal Analysis - Eigenvalue-EigenVector function
4. Transmission Zeros
5. SYSTEM TRANSFER FUNCTIONS: From u_i to y_j
6. Controllability and Observability
7. Modal analysis
8. Frequency Response: Singular Values of Plant
9. Augment Plant with Integrators at Plant Input and Plot Singular Values
10. Bilinear transformation
11. LQR Design
12. Kalman filter
13. LQG Design
14. LQG-LTR Design
15. H-infinity Robust Constrol Design
16. CLOSED LOOP TIME RESPONSE 

# Background

The McDonnell Douglas(now Boeing) AV-8A Harrier is a single engine ground attack aircraft which is capable of vertical or short takeoff and landing. It was developed in the 1960s and formed the first generation of the Harrier series of aircrafts. It is powered by a single Pegasus turbofan engine mounted in the fuselage. The engine is fitted with four vectoring nozzles for directing the thrust generated(two for the bypass flow and two for the jet exhaust) and two air intakes. The aircraft also has several smaller reaction nozzles in the nose , tail and wingtips for the purpose of balancing during vertical flight. The aircraft is capable of forward flight like a fixed wing aircraft. It is also capable of doing VTOL and STOL manoeuvres where the lift and control surfaces are useless. The harrier also has two control elements namely the thrust vector and the reaction control system which is not found in conventional fixed-wing aircraft.
