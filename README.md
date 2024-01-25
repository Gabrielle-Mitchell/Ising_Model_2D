# Ising_Model_2D
The Ising model of a ferromagnet is a lattice model used to interpret ferromagnetism in materials. In this simulation, the ferromagnet was modelled as a two dimensional square lattice with dipole spins at each lattice point, pointing either up or down on the y-axis. Various observables, including the magnetization, energy, and heat capacity were calculated as a function of temperature to show the temperature at which the phase transition occurs. This is known as the Curie temperature, and it is the temperature at which a ferromagnet becomes paramagnetic and its spins become randomly aligned. This was accomplished using the Metropolis-Hastings Algorithm, which uses a Metropolis loop to decide whether or not a randomly chosen spin should be flipped. The Monte Carlo loop, which contains the Metropolis loop, was used to simulate 10,000 different lattice representative states per temperature, and the average of the accumulated change in each observable was taken. The heat capacity graphs were used to estimate the Curie temperature for a variety of lattice sizes. The code is primarily adapted from [2].


References

[1] Dieterle and Witthauer. (2007). The Phase Transition of the 2D-Ising Model. Retrieved from http://quantumtheory.physik.unibas.ch/people/bruder/Semesterprojekte2007/p1/Ising.pdf

[2] Jacques Kotze. (2008). Introduction to Monte Carlo Methods for an Ising Model of a Ferromagnet. Cornell University Library. Retrieved from https://arxiv.org/abs/0803.0217

[3] Ramaswamy and Shrinivas. (2013). Metropolis Monte Carlo Simulation of the Ising Model. Massachusetts Institute of Technology. Retrieved from http://web.mit.edu/krish_s/www/files/ising_Model.pdf

[4] Ben Hammel. (2014). Monte Carlo Simulation of the Ising Model using Python. Retrived from http://www.thebrokendesk.com/post/monte-carlo-simulation-of-the-ising-model-using- python/

[5] Sedgewick and Wayne. (2005). Modular Programming. Princeton University. Retrieved from http://introcs.cs.princeton.edu/java/lectures/98ising.pdf

[6] Fraden, Kosowsky, Steinberg. (2013). Simulations: The Ising Model. Brandies University, Waltham. Retrieved from http://fraden.brandeis.edu/courses/phys39/simulations/AsherIsingModelReport.pdfhttp://fraden.brandeis.edu/courses/phys39/simulations/AsherIsingModelReport.pdf

[7] Yu, Sun, Alber. (2015), Monte Carlo Simulation of Ising Model and Phase Transition Studies. University of Notre Dame. Retrieved from http://slideplayer.com/slide/5191426/
