# Code implementation
## Implement RG in C++ code

- FTODO: Mira de hacer un loop con las diferentes $b$. Guardas los pointers o datasets en un array y los nombres van puestos segun el array de las $b$.

- FTODO: Quiza hacer Hoshen-Kopelman para encontrar tamaño de los clusters para rg.

## Vortices
- TODO: Add how to compute the vorticity (winding number over plaquettes)

## Restructuration
- Make it so that array of $b$ is more interactive and not just copy paste.


## Hybrid method
- Incorpora método híbrido, a la que haya una cantidad de clusters estándar, cambia al método de Metropolis y solo cambia de vuelta si todos los spins son cambiados $\implies$ faster update.

# Observables implementation
- Add thermalisation check
    - Add some standard deviation check or something.**??**
- Bear in mind that correlation time goes up as lattice size goes up by Wolff method (see Newman page 101).
- Read fully finite size difference stuff to see if that method can work.

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
- Look at test.ipynb in the Documents/test folder to see how to implement the method to find the inflection point and the error to make the extrapolation.


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
