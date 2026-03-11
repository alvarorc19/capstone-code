# Code implementation
- Implement RG in C++ code
- Make code saving the last 100 lattice steps? Have done it for 
    1000 but might occupy too much storage. Let me try simulation with only 100.
    - Create a separate dataset that its only job is to store the last 10 states
    - This means create a new function to write the lattice at the end of simulation

- Newman and Barkema page 102 mira que hay una manera de calcular la susceptibility con el dinamical exponent gracias a lo relacionado que esta la correlation length con esto.

## Major rework
- Haz todo el recorrido mirando que declaraciones son innecesarias, ocupan memoria que no deberia cuando se puede hacer un acceso simple mediante pointer que hace lo mismo.
- Implementa el cálculo de observables a poca escala para el método de renormalización.
- Incorpora método híbrido, a la que haya una cantidad de clusters estándar, cambia al método de Metropolis y solo cambia de vuelta si todos los spins son cambiados $\implies$ faster update.
 

# Observables implementation
- Add thermalisation check
    - Save all of the magnetisations to see when it actually thermalises
    - Add some standard deviation check or something.

## H5 file structure

├── observables/
│   ├── energy
│   ├── susceptibility
│   ├── correlation function?(t and r)
│   └── magnetisation (x and y column)
├── renormalised_observables/ (do all of the metrics but include different values of $b$). If you do this one not do the other one.
│   ├── energy
│   ├── susceptibility
│   ├── correlation function?(t and r)
│   └── magnetisation (x and y column)
├── clusters/
│   ├── number of clusters
│   └── size of clusters (could do an N by time matrix and store there how big the clusters made are)
└── lattice (just keep last sweep

Keep this here:
│

└── 

├── 
