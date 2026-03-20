# Code implementation
## Implement RG in C++ code
- Make function that calculates the reduced observables
1. Calculate neighbour table for the necessary $b$ s
l1. Execute 512 sweeps
l2. Compute rg matrix (from XYModel)
l3. Compute magnetisation (from XYModel)
l4. Compute energy (from XYModel)
l5. Write it into h5 (make function for this)
2. Do last sweep < 512 


- I did the code but it seems like the thing does not work. The energies it obtains are all over the place. Need to further look into algorithm
- Instead of doing only blocking, mix finite size and blocking and plot everything together to observe the crossings. Do that for simple case of $16\times16$ and $32\times32$ blocking and actually reducing the lattice size, and plot them all together. See Newman and Barkema RG for magnetisation.

- Mira de hacer un loop con las diferentes $b$. Guardas los pointers o datasets en un array y los nombres van puestos segun el array de las $b$.

- Quiza hacer Hoshen-Kopelman para encontrar tamaño de los clusters para rg.

## Restructuration
- Make it so that array of $b$ is more interactive and not just copy paste.


## Hybrid method
- Incorpora método híbrido, a la que haya una cantidad de clusters estándar, cambia al método de Metropolis y solo cambia de vuelta si todos los spins son cambiados $\implies$ faster update.

# Observables implementation
- Add thermalisation check
    - Add some standard deviation check or something.**??**
- Bear in mind that correlation time goes up as lattice size goes up by Wolff method (see Newman page 101).

## Python
- Compute correlation function of the lattices I get at the end
- With average number of spins in cluster compute correlation length and susceptibility
    - Compare with the normal one
    - Do it in the file observables_plots, inside the big function.
- Plot thermalisation but only first steps to have more resolution and see where it changes
- Extract $\tau_{int}$ values to see where the error comes from, especially in specific heat and susceptibility (I think there's something related to $N$).
    - Maybe deal with autocorrelation function. Look at Newman page 60 something. Also here [here](../../../../../Downloads/BF01022990-1.pdf). This is the one used in `pyerrors`.
    - Really look in depth in Obs class object to see what stuff it does.
- **Error analysis**
    - See how jacknife and bootstrap are implemented, and compare it to them.


## H5 file structure

├── observables/
│   ├── energy
<!-- │   ├── susceptibility -->
│   └── magnetisation (x and y column)
├── renormalisation/ (do all of the metrics but include different values of $b$). If you do this one not do the other one.
│   ├── energy
<!-- │   ├── susceptibility ? necessary -->
│   └── magnetisation (x and y column)
├── clusters/
<!-- │   ├── number of clusters -->
│   └── average size of clusters
└── lattice (just keep last 100 steps)

Keep this here:
│

└── 

├── 
