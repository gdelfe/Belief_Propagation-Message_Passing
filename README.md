# belief-propagation-message-passing
This repo contains the codes to compute the magnetization at a given temperature for a Ising spin model on a Erdős–Rényi random graph (ERRG) via **message-passing (belief-propagatation/cavity method)**. </br>
Results are then compared with the magnetization obtained via Monte Carlo simulations.

Some reference to the methods used:</br>

The algorithm is known as "cavity method" in the physics community and "message-passing" or "belief-propagation" in the computer science community. </br>
1. <a href="https://arxiv.org/abs/1409.3048">Cavity Method: Message Passing from a Physics Perspective</a> </br>
Gino Del Ferraro, Chuang Wang, Dani Martí, Marc Mézard
2. <a href=https://www.diva-portal.org/smash/get/diva2:957675/FULLTEXT01.pdf>Equilibrium and dynamics on complex networks</a> </br>
Gino Del Ferraro  - (Chapter 2. - see ref. therein)

To run the code:
1. Generate a Erdos Renyi random graph with Generate_ERRG_deg.c, follow the instruction at the top of the code, comment section
2. Run static_BP.c (follow the instruction at the top part of the code) to compute the magnetization of an Ising model with arbitrary connectivity as a function of the temperature.
3. If you only want to run the code in 2. for one given temeperature, you can run static_BP_single_beta.c

All the codes run for a given instance of the ERRG Ising model. If you want to get results for multiple realization of the ERRG you will have to average the results over different run, each made on a different ERRG. 

Please cite this work if you use any of these codes: 

Del Ferraro, G. (2016). Equilibrium and Dynamics on Complex Networkds. PhD dissertation, KTH Royal Institute of Technology.

Del Ferraro, G., Wang, C., Martí, D., & Mézard, M. (2014). Cavity method: Message passing from a physics perspective. arXiv preprint arXiv:1409.3048.

