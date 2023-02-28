# Notes on Particle Optimization

## References
[pyswarms collection of resources](http://particleswarm.info/index.html)
[github main repo](https://github.com/ljvmiranda921/pyswarms)


**Optimization**

[TWD article: PSO visually explained](https://towardsdatascience.com/particle-swarm-optimization-visually-explained-46289eeb2e14)
[pyswarms - options handler tutorial](https://pyswarms.readthedocs.io/en/latest/examples/tutorials/options_handler.html)

> As you might have noticed, I have not yet talked about the *inertia*, *cognitive* and *social coefficients*. These coefficients control the levels of exploration and exploitation. Exploitation is the ability of particles to target the best solutions found so far. Exploration, on the other hand, is the ability of particles to evaluate the entire research space. The challenge of the remaining part of the article will be to determine the impact of these coefficients to find a good balance between exploration and exploitation.

**Inertia *w***
> The inertia weight w thus makes a balance between the exploration and the exploitation of the best solutions found so far.

---
> According to the paper by M. Clerc and J. Kennedy to define a standard for Particle Swarm Optimization, the best static parameters are w=0.72984 and c1 + c2 > 4. More exactly c1 = c2 = 2.05. Additionally, the linear decay of the parameter w was initially proposed by Yuhui and Russ Y. H. Shi and R. C. Eberhart.



![](/assets/images/san-juan-mountains.jpg "")

**Visualization**
[pyswarms - plotters](https://pyswarms.readthedocs.io/en/latest/api/pyswarms.utils.plotters.html#module-pyswarms.utils.plotters)
[pyswarms - plotting in 2D space](https://pyswarms.readthedocs.io/en/latest/examples/tutorials/options_handler.html#Plotting-in-2-D-space)

**Features: benchmark functions, plotters, search options for hyperparams**
[pyswarms - features](https://pyswarms.readthedocs.io/en/latest/features.html)
