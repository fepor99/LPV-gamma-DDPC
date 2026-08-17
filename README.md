# Reduced LPV $\gamma$-DDPC
Code of the reduced LPV $\gamma$-DDPC controller implemented in the paper
> Federico Porcari, Chris Verhoek, Roland Tóth, Valentina Breschi, Simone Formentin.
> "A subspace approach to data-driven predictive control for linear parameter varying systems".
> _arXiv preprint arXiv:2607.28490_, 2026

Released under the BSD 3-Clause License.

Please cite us if you use this code for any work development:
```bibtex
@article{Porcari2026:LPV-gamma-DDPC,
  title     = {A subspace approach to data-driven predictive control for linear parameter varying systems},
  author    = {Porcari, Federico and Verhoek, Chris and T{\'o}th, Roland and Breschi, Valentina and Formentin, Simone},
  note      = {arXiv preprint arXiv:2607.28490},
  year      = {2026},
  journal   = {Submitted to Automatica}
}
```

For any questions related to the code, you can contact me via email at <federico.porcari@polimi.it>

### Usage
To run the code, one needs MATLAB (tested in 2020b) with the `Signal Processing Toolbox`, `YALMIP`, and the `GUROBI` solver. The latter can be changed for another SDP solver, but it may obtain slightly different results compared to the ones presented in the paper. 

Although optional, the offline LQ decomposition is computed through C functions if the `MATLAB Coder Toolbox` is available. 

### Acknowledgements
The code for the LPV-IO-DPC controller has been taken from [its original GitLab repository](https://gitlab.com/releases-c-verhoek/lpvdpc) by Chris Verhoek.
