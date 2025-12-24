# Characterisation of the XY model

## Ising model
Defined by Hamiltonian\
$$H = E = -J \sum_{\langle i j \rangle}s_i s_j - h \sum_i s_i.$$ 

## Simulation steps
1. Make a parser and have hold the relevant data for the simulation
2. Initialise lattice
3. Run simulation
    i. Measure current state
    ii. Perform algorithm
    iii. Write new state (tell how verbose in an option, also log)
    iv. Compute observables or measurables (also tell if do this in config)
4. Quit simulation and end program

## Issues or code choices to address
- Make one lattice class?
- For the models make a general model class and have the other ones inherit?
    - Could make an $O(1)$ case and then an $O(n)$ for $n\geq2$ because of the implementation and usage of vectors
- Think how to do the magnetic field, difference between $O(1)$ and $O(2)$. Make sure to include assertions. 
- How am I going to store the values for the XY model?
  - Have an array way bigger s.t. they are next to each other?
  - Have it represented as a complex number?
  - Have only the number and then compute the sine and cosine every time?
    - I feel like this one is non-efficient!
    - However, if Hamiltonian is w.r.t. angle...

## TODO
- Implement table of neighbours in 2D
    - Use the 1d array that where for an index i, for 2d neighbours are stored at $4\dot i + \textrm{dir}$, where dir is direction and can be $[0,3]$, up, down, right, left. The index i is defined as index $= i\dot L +j$ for a position in lattice $(i,j)$.
- Change computation of energy 
    - Currently $\Delta E = E(s_f)-E(s_i)$, for a spin flip $s_i \rightarrow s_f$ change it to a more general scenario for the Potts' model (might need change for the XY model due to dot product)\
    $$\Delta E = -J (s_f-s_i)\sum_{j}s_j - h (s_f-s_i).$$ \
    Where the sum over $j$ is over the neighbours
- Add Hamiltonian for Potts' model
- Add Hamiltonian for XY model
- Code cluster algorithm
- Make parser that reads from the terminal and then from a toml file
    - This is in conjunction to creating a folder where the hdf5 file would be dumped as well as the log file
