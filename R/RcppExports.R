# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rcpp_numerov <- function() {
    .Call('numerov_rcpp_numerov', PACKAGE = 'numerov')
}

getEnergiesAndIndices <- function() {
    .Call('numerov_getEnergiesAndIndices', PACKAGE = 'numerov')
}

setPotential <- function(px = numeric(), py = numeric()) {
    invisible(.Call('numerov_setPotential', PACKAGE = 'numerov', px, py))
}

getPotential <- function() {
    .Call('numerov_getPotential', PACKAGE = 'numerov')
}

computeSpectrum <- function(nEigen, dE = 0.1, tol = 1e-9) {
    invisible(.Call('numerov_computeSpectrum', PACKAGE = 'numerov', nEigen, dE, tol))
}

getEnergies <- function() {
    .Call('numerov_getEnergies', PACKAGE = 'numerov')
}

getWavefunctions <- function() {
    .Call('numerov_getWavefunctions', PACKAGE = 'numerov')
}

