# c-sensor-registration
* <a href="https://arxiv.org/pdf/2001.04286.pdf" target="_blank">Nonparametric Continuous Sensor Registration.</a>

This repository contains examples of sensor registration using different manifolds and Lie groups. For the RGB-D visual odometry case, i.e., R^3 and SE(3), see:
<a href="https://github.com/MaaniGhaffari/cvo-rgbd" target="_blank">Continuous Direct Sparse Visual Odometry from RGB-D Images.</a>

**Continuous sensor registration** is a new mathematical framework that enables nonparametric joint semantic/appearance and geometric representation of continuous functions using data. The joint semantic and geometric embedding is modeled by representing the processes in a reproducing kernel Hilbert space. The framework allows the functions to be defined on arbitrary smooth manifolds where the action of a Lie group is used to align them. The continuous functions allow the registration to be independent of a specific signal resolution and the framework is fully analytical with a closed-form derivation of the Riemannian gradient and Hessian.

## Citations
* William Clark, Maani Ghaffari, Anthony Bloch. "Nonparametric Continuous Sensor Registration." Journal of Machine Learning Research 22.271 (2021): 1-50. http://jmlr.org/papers/v22/20-1468.html
```
@article{JMLR:v22:20-1468,
  author  = {William Clark and Maani Ghaffari and Anthony Bloch},
  title   = {Nonparametric Continuous Sensor Registration},
  journal = {Journal of Machine Learning Research},
  year    = {2021},
  volume  = {22},
  number  = {271},
  pages   = {1-50},
  url     = {http://jmlr.org/papers/v22/20-1468.html}
}
```
