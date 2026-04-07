# Code implementation
## Implement RG in C++ code

- FTODO: Mira de hacer un loop con las diferentes $b$. Guardas los pointers o datasets en un array y los nombres van puestos segun el array de las $b$.

- FTODO: Quiza hacer Hoshen-Kopelman para encontrar tamaño de los clusters para rg.

## Vortices
- FTODO: Add how to compute the vorticity (winding number over plaquettes)

## Restructuration
- FTODO: Make it so that array of $b$ is more interactive and not just copy paste.


# Observables implementation
- Bear in mind that correlation time goes up as lattice size goes up by Wolff method (see Newman page 101).
- Read fully finite size difference stuff to see if that method can work.

## Python
- FTODO: Compute correlation function of the lattices I get at the end

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
