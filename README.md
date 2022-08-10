# Blim

Estimation of Blim for ICES stocks following guidelines from ICES 2021 (https://doi.org/10.17895/ices.advice.7891)

Can be installed by typing: 

```R
devtools::install_github("https://github.com/vtrijoulet/Blim/Blim")
```

## Example

A simple example can be run with the following code:

```R
library(Blim)
example("Blim_estim")
```

The example code is as follows:

```R
library(Blim)

s <- seq(100, 2000, length=50)
r <- exp((2+log(s)-0.001*s)+rnorm(length(s)))
data <- data.frame(Rec=r, SSB=s, row.names=1971:2020)

# The stock has spasmodic recruitment (type 1):
Blim_estim(data, Rage=0, is.spasmodic=TRUE, doplot=TRUE) # only type 1 is assumed

# The algorithm tests for ICES stock types 2-6 and chooses the best:
Blim_estim(data, Rage=0, doplot=TRUE) 

# Type 2 is not considered and the algorithm only tests for types 3-6:
Blim_estim(data, Rage=0, types=3:6, doplot=TRUE) 

```

If the stock is believed to have sporadic recruitment, the argument "is.spasmodic" should be set to TRUE, and Blim will be estimated according to the ICES type 1.

If the stock does not have sporadic recruitment, Blim can be estimated by running the "Blim_estim" function with the default arguments. In this case types 2-6 will be tested and the best type and corresponding Blim will be reported.

It is possible to avoid testing for type 2 (fitting the hockey-stick stock-recruitment relationship) by modifying the "types" argument so it does not contain the number 2.


More details available by running:
```R
library(Blim)
?Blim_estim
```
