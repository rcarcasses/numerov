
plotQuantumGraph <- function(nEigen) {
    getEnergies();
}

quantumAnalysis <- function(px = seq(-5, 5, by = 0.01), 
                            py = seq(-5, 5, by = 1),
                            nEigen = 6,
                            scale = 60,
                            dE = 0.1, 
                            tol = 1e-9) {
    if(length(px) != length(py)) {
        setPotential();
    } else {
        setPotential(px, py);
    }
    
    computeSpectrum(nEigen, dE, tol);
    wfs      = getWavefunctions();
    energies = getEnergies();
    potential = getPotential();
    # make a nice quantum graph
    plot.new();
    plot(potential$x, potential$y, xlab = 'x', ylab = 'V(x)', type = 'l', lwd = 3);
    
    cl <- rainbow(length(wfs));
    for(i in 1:length(wfs)) {
        lines(wfs[[i]]$x, scale * wfs[[i]]$y + energies[i], col = cl[i], lwd = 2);
    }
}