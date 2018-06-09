# Implementation

- `odefun_control1.m`: Geometric control on SE(3)xR^2 quadrotor suspended with point mass. which illustrates the Matlab simulation and implementation of 'Geometric Control and Differential Flatness of a Quadrotor UAV with a Cable-Suspended Load - Koushil Sreenath, Taeyoung Lee, Vijay Kumar' link https://ieeexplore.ieee.org/document/6760219/

- `odefun_control2.m`: Apply previous geometric control on SE(3)xR^2 quadrotor suspended with point mass into the offset dynamical model, where we can observe a final configuration error

# Robust Analysis of geometric control into outset dynamical model
![Old Control Performance](/ControlPerformanceComparaison/OldModel.png )
where we observe significant final configuration error, not convergent to zero.

![Improved Control Performance](/ControlPerformanceComparaison/ImprovedModel.png)
where the configuration error converges perfectly.
