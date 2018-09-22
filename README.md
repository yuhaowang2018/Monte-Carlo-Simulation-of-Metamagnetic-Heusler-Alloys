# Monte Carlo Simulation of Heusler Alloys
## Usage:
   The Monte Carlo simulation's implementation is inside "simulation" folder. To calibrate the model parameters, compile the c++ program and put the excutable inside the "parameter_tuning" folder and run the Python program.

## Introduction:
* Implemented Monte Carlo simulation with Metropolis algorithm using C++.
* Potts model was used based on Hamiltonian with interchanges between elastic, chemical and magnetic free energy.
* Gaussian process and efficient global optimization were used to efficiently calibrate model parameters.

## Sample output for Ni-Co-Mn-In
* Strain distortion:

![alt text][plot1]

* Magnetic transformation:

![alt text][plot2]

[plot1]: https://github.com/yuhaowang2018/Monte-Carlo-Simulation-of-Metamagnetic-Heusler-Alloys/blob/master/strain1.png
[plot2]: https://github.com/yuhaowang2018/Monte-Carlo-Simulation-of-Metamagnetic-Heusler-Alloys/blob/master/Full_H_2.png
