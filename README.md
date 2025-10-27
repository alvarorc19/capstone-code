# Characterisation of the XY model

## Ising model
Defined by Hamiltonian\
$$
H = E = -J \sum_{\langle i j \rangle}s_i s_j - h \sum_i s_i
.$$ 

## TODO
- [ ] Implement table of neighbours in 2D
    - Use the 1d array that where for an index i, for 2d neighbours are stored at $4\dot i + \textrm{dir}$, where dir is direction and can be $[0,3]$, up, down, right, left. The index i is defined as index $= i\dot L +j$ for a position in lattice $(i,j)$.
- [ ] Change computation of energy 
    - Currently $\Delta E = E_f-E_i$, change it to a more general scenario for the Potts' model (might need change for the XY model due to dot product)\
    $$
    \Delta E = -J (s_f-s_i)\sum_{j}s_j - h (s_f-s_i)    
    .$$ 
    Where the sum over $j$ is over the neighbours
- [ ] Add Hamiltonian for Potts' model
- [ ] Add Hamiltonian for XY model
- [ ] Code cluster algorithm
