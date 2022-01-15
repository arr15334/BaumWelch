# Baum-Welch

R package of manual implementation of Baum-Welch algorithm

## Installation

Run the following command from RStudio console

`devtools::install_github('arr15334/BaumWelch')`

And then load the library using

`library(baumw)`

## Usage

`hmm = BaumWelch(visible, transitionMatrix, emissionMatrix, initial_distribution, iterations)`

-`visible` is a vector of the emitted symbols
-`transitionMatrix` is the first estimation of the transition matrix A. This can be randomly initialized if there is no prior information
-`emissionMatrix` is the first estimation of the emission matrix A. This can be randomly initialized if there is no prior information
-`initial_distribution` is the first estimation of the initial distribution pi.
-`iterations` is the number of maximum iterations that the algorithm should perform.

This will return a list containing the estimations of the transition and emission matrices. 