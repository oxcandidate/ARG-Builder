# Arg-Builder and KwARG

Software implementing parsimony-based algorithms for reconstructing ancestral recombination graphs (ARGs) with recurrent mutation. 
There are two algorithms: KwARG and ARG-Builder.
The KwARG algorithm here is an updated version to that available at https://github.com/a-ignatieva/kwarg.git

## Citation
KwARG:
Ignatieva, A., Lyngs&oslash;, R. B., Jenkins, P. A., and Hein, J. KwARG: Parsimonious reconstruction of ancestral recombination graphs with recurrent       mutation. [[arXiv]](http://arxiv.org/abs/2012.09562) [[bioRxiv]](https://www.biorxiv.org/content/10.1101/2020.12.17.423233v1)

## License information
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Installation

This repository contains the programs KwARG and ARG-Builder as well as several python scripts for calculating lower bounds on the number of recombinations required for a dataset. 


## Compile from source

Clone the repository using git via command line:
```sh
git clone https://github.com/oxcandidate/ARG-Builder.git
```
### Linux and Mac OS

Open Terminal and navigate to the source folder for either kwarg of arg_builder:

The following command will compile:
```sh
make all
```

The following command will delete the compiled objects:
```sh
make clean
```
