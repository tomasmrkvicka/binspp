# binspp
Bayesian inference for Neyman-Scott point processes (R package)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

R software with VGAM, spatstat, FNN and cluster libraries installed. Additional required packages are: Rcpp, RcppArmadillo and RcppEigen.

```
In R shell write:

> install.packages(c("VGAM", "spatstat", "FNN", "cluster"))
> install.packages(c("Rcpp", "RcppArmadillo", "RcppEigen"))
```

### Installing

Install package by downloading from CRAN:
```
install.packages("binspp")
```
You can also download binspp.tar.gz package and install it to your R software.

```
install.packages("C:/path/to/directory/binspp.tar.gz", 
    repos = NULL, 
    lib = "C:/path/to/libraryDirectory")
```

## Running the tests

Load data dataset_N4.Rdata, run example scripts to test package functionality.

## Built With

R Studio or any other R software.

* [RStudio](https://rstudio.com/products/rstudio/download/) - The R Studio

## Versioning

We use [GitHub](http://github.com/) for versioning. For the versions available, see the [binspp](https://github.com/tomasmrkvicka/binspp). 
You can also get binspp package on the [CRAN](https://cran.r-project.org/package=binspp).

## Authors

* **Tomas Mrkvicka** - *creator*, *author* - [ResearchGate](https://www.researchgate.net/profile/Tomas_Mrkvicka)
* **Jiri Dvorak** - *author* - [ResearchGate](https://www.researchgate.net/profile/Jiri-Dvorak-5)
* **Ladislav Beranek** - *author* - [GitHub](https://github.com/lberanek)
* **Radim Remes** - *author*, *maintainer* - [GitHub](https://github.com/radimremes)

See also the list of [contributors](https://github.com/tomasmrkvicka/binspp/contributors) who participated in this project.

## License

This project is licensed under the GNU GPL 3 License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

Anderson, C. Mrkvička T. (2020). Inference for cluster point processes with over- or under-dispersed cluster sizes, *Statistics and computing* **30**, 1573–1590. https://doi.org/10.1007/s11222-020-09960-8
