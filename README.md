inmed2014
=========

##International Conference on Innovation in Medicine and Healthcare 2014

This repository contains both the both the article and the presentation slides for the work entitled: **Gated Sensor Fusion: A way to Improve the Precision of Ambulatory Human Body Motion Estimation** which was accepted for presentation at the INMED2014 *(San Sebastian, Spain, 9-11 July)*.

**ABSTRACT:** Human body motion is usually variable in terms of intensity and, therefore, any Inertial Measurement Unit attached to a subject will measure both low and high angular rate and accelerations. This can be a problem for the accuracy of orientation estimation algorithms based on adaptive filters such as the Kalman filter, since both the variances of the process noise and the measurement noise are set at the beginning of the algorithm and remain constant during its execution. Setting fixed noise parameters burdens the adaptation capability of the filter if the intensity of the motion changes rapidly. In this work we present a conjoint novel algorithm which uses a motion intensity detector to dynamically vary the noise statistical parameters of different approaches of the Kalman filter. Results show that the precision of the estimated orientation in terms of the RMSE can be improved up to 29% with respect to the standard fixed-parameters approaches.

**KEYWORDS:** Human motion, Orientation estimation, MARG sensors, Kalman filtering, Intensity detection.

- To compile both the presentation slides and the article use DVI-PS-PDF chain in conjunction with LaTeX.
