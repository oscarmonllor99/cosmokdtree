# RAPID-kdtree

## About

Developed at the Departament d'Astronomia i Astrofísica of Universitat de València by Óscar Monllor-Berbegal in collaboration with David Vallés-Pérez, Susana Planelles and Vicent Quilis. This work has been supported by the European Union NextGenerationEU (PRTR-C17.I1), the Spanish Ministerio de Ciencia e Innovación(ASFAE/2022/001 and PID2022-138855NB-C33), the Generalitat Valenciana (CIPROM/2022/49), and Óscar Monllor-Berbegal acknowledges support from Universitat de València through an Atracció de Talent fellowship.

## Index of contents 

1. [Brief description](#brief-description)

2. [Repository organisation](#repository-organisation)

3. [Installation](#installation)

4. [Running the code](#running-the-code)

## Brief description

Fortran 2003/2008 $k$-d tree implementation with OpenMP directives for parallel tree construction. Designed for efficient spatial indexing and nearest-neighbour searches in large simulation datasets. Below, we provide a brief description of the different techniques used to build the method:

* The maximum variance splitting is used to choose the axis across which the input dataset is split at each tree depth.
* *Quickselect* is employed to find the median point dividing the data into two parts across the splitting axis. A *lazy* approach is followed, as the algorithm is not required to find the exact median point, but a close candidate inside the (45%, 55%) range.
* An adaptive leaf size is used, scaling with the number of input points. Faster tree construction and queries are obtained with this technique.
* Tree construction is parallelised leveraging *OpenMP task* constructs, until a maximum tree depth ($d_\text{max}$) is reached and all physical cores ($N_\text{CPU}$) are being used.
* The *Max Heap* structure is utilised to efficiently save and update the $k$-nearest neighbours, as the tree is traversed.
* In case of needing sorted results, a *Quicksort* implementation is leveraged.
* The algorithm is fully customisable and versatile, as the user can specify the integer size, floating point arithmetics precision, periodic boundary conditions, or dimensionality of the input space at compilation time.

A detailed description of the algorithm, together with several tests demonstrating its performance, can be found in _Monllor-Berbegal et al. 2025 in prep._

## Repository organisation

The source code and Makefile can be found inside the `src` folder. All demos utilised to test the tree construction and query routines are inside the `demos` folder, together with proper compilation and executable options. Lastly, all Jupyter Notebooks used to save the test results and plot the benchmarks can be found in the `benchmarks` folder, together with a friends-of-friends implementation testing the $k$-d tree efficiency when applied to this particular issue.

## Installation

### Download

Clone the repository into the desired directory by running:

```
git clone https://github.com/oscarmonllor99/cosmo_kdtree.git
```

Then access the root directory:

```
cd cosmo_kdtree-main
```

### Dependencies

Our code needs few dependencies: a Fortran compiler and OpenMP for shared-memory parallelisation, as it is practically self-contained.

* [gfortran](https://gcc.gnu.org/wiki/GFortran)
* [OpenMP](https://www.openmp.org/)

### Make (Compilation)

We already provide the user with an example Makefile to compile the code inside `src`. A call to make should look like this:

```
make PERIODIC=A DIMEN=B LONGINT=C DOUBLEPRECISION=D
```

`PERIODIC` tells the code whether to use periodic boundary conditions or not (`0` deactivates them, while `1` will switch them on). Defaults to `0`.

`DIMEN` sets the dimensionality of the points provided as input. Hence `DIMEN` should be an integer equal to or greater than `1`. Defaults to `3`.

`LONGINT` specifies the size of integers employed to index points in the $k$-d tree. If more than $N = 2,147,483,648$ are used, `LONGINT=1` is necessary. Defaults to `0`.

`DOUBLEPRECISION` sets the floating point arithmetic precision to `REAL*4` if it is `0`, or to `REAL*8` if `1`. If high precision is needed, it should be activated. Defaults to `0`.


## Running the code inside Fortran

As *cosmokdtree* builds the $k$-d tree and saves query results using `pointers` and `derived types` (KDTreeNode and KDTreeResult), its usage is quite intuitive and straightforward. Inside the `demos` folders we provided many examples displaying the wide range of possibilities in which the $k$-d tree could be applied, depending on the compilation options. The tree must be declared using the following directive:

```
type(KDTreeNode), pointer :: tree
```

Then, calling the building function is as easy as:

```
tree => build_kdtree(points)
```

After that, carrying out a query on the tree is simple. Firstly, the query variable must be declared:

```
type(KDTreeResult) :: query
```

Then, when can call the query function and recover the results:

```
query = knn_search(tree, target_point, k, sorted)
indices = query%idx
distances = query%dist
```

Several query functions are available: `knn_search`, `ball_search` and `box_search`. The first finds the $k$-nearest neighbours from a point, the second finds all neighbours within a $R$ distance from a point and, the last one finds all points inside a query box. In the first two cases, we provide the option (`sorted = .true.`) to sort results by ascending distance to the query point. In all cases, if `PERIODIC=1`, periodic boundary conditions will be taken into account.






