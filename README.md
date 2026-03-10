
## About

 This work has been supported by the Agencia Estatal de Investigación Española (AEI; grant PID2022-138855NB-C33), by the Ministerio de Ciencia e Innovación (MCIN) within the Plan de Recuperación, Transformación y Resiliencia del Gobierno de España through the project ASFAE/2022/001, with funding from European Union NextGenerationEU (PRTR-C17.I1), and by the Generalitat Valenciana (grant PROMETEO CIPROM/2022/49). OM and DV acknowledge support from Universitat de València through Atracció de Talent fellowships. DV acknowledges additional support from the ERC CoG $\vec{B}$ELOVED, GA n. 101169773. Simulations have been carried out using the supercomputer Lluís Vives at the Servei d'Informàtica of the Universitat de València.

## Index of contents 

1. [Brief description](#brief-description)

2. [Repository organisation](#repository-organisation)

3. [Installation](#installation)

4. [Running the code](#running-the-code)

## Brief description

Fortran 2003/2008 $k$-d tree implementation with OpenMP directives for parallel tree construction. Designed for efficient spatial indexing and nearest-neighbour searches of large datasets. Below, we provide a brief description of the different techniques used to build the method:

* The sliding-midpoint splitting is used to choose the axis across which the input dataset is split at each tree depth.
* A configurable leaf size is used, defaulting to $N_\text{leaf} = 16$.
* Tree top-down construction is parallelized leveraging *OpenMP task* constructs, until a maximum tree depth is reached and all physical cores ($N_\text{CPU}$) are being used.
* The *Max Heap* and *Min Heap* structures are utilized to efficiently traverse the tree in nearest-neighbor queries. 
* In case of needing sorted results, a *Quicksort* implementation is leveraged.
* The algorithm is fully customisable and versatile, as the user can specify the integer size, floating point arithmetics precision, periodic boundary conditions, or dimensionality of the input space at compilation time.
* Python bindings are also included.

A detailed description of the algorithm, together with several tests demonstrating its performance, can be found in _Monllor-Berbegal et al. 2026 in prep._

## Repository organisation

The source code and Makefile can be found inside the `src` folder. All demos used to test the tree construction and query routines are inside the `demos` folder, together with suitable compilation and executable options. All scripts used to carry out Python benchmarks can be found in the `benchmarks` folder.

## Installation

### Download

Clone the repository into the desired directory by running:

```
git clone https://github.com/oscarmonllor99/cosmokdtree.git
```

Then access the root directory:

```
cd cosmokdtree-main
```

### Dependencies

Our code needs few dependencies: a Fortran compiler and OpenMP for shared-memory parallelization, as it is practically self-contained.

* [gfortran](https://gcc.gnu.org/wiki/GFortran)
* [OpenMP](https://www.openmp.org/)

### Make (Compilation)

We already provide the user with an example Makefile to compile the code inside `src`. A call to make should look like this:

```
make PERIODIC=A DIMEN=B LONGINT=C DOUBLEPRECISION=D PYTHON=E
```

`PERIODIC` tells the code whether to use periodic boundary conditions or not (`0` deactivates them, while `1` will switch them on). Defaults to `0`.

`DIMEN` sets the dimensionality of the points provided as input. Hence `DIMEN` should be an integer equal to or greater than `1`. Defaults to `3`.

`LONGINT` specifies the size of integers employed to index points in the $k$-d tree. If more than $N = 2,147,483,648$ are used, `LONGINT=1` is necessary. Defaults to `0`.

`DOUBLEPRECISION` sets the floating point arithmetic precision to `REAL*4` if it is `0`, or to `REAL*8` if `1`. If high precision is needed, it should be activated. Defaults to `0`.

`PYTHON` tells the code whether to create the Python module. If `1`, it fixes LONGINT=1 and DOUBLEPRECISION=1 and calls F2PY to create the module file that is imported by `kdtree.py` inside `/src/cosmokdtree`.


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

Then, one can call the query function and recover the results:

```
query = knn_search(tree, target_point, k, sorted)
indices = query%idx
distances = query%dist
```

Several query functions are available: `knn_search`, `ball_search` and `box_search`. The first finds the $k$-nearest neighbours from a point, the second finds all neighbours within a $R$ distance from a point and, the last one finds all points inside a query box. In the first two cases, we provide the option (`sorted = .true.`) to sort results by ascending distance to the query point. In all cases, if `PERIODIC=1`, periodic boundary conditions will be taken into account.

## Running the code inside Python

First of all, the code must be compile using `PYTHON=1` in the `make` call. After that, a `pycosmokdtree.cpython` file is created inside `src`, which is the module object called by `kdtree.py` inside the `src/cosmokdtree` Python module folder. Note that, if the user wants to use the Python module in an arbitrary folder, this folder must be added to the system `PATH`. 

Inside `kdtree.py` the user will find the `build_kdtree`, `knn_search`, `ball_search` and `box_search` functions working exactly as the ones in the Fortran module:

```
build_kdtree(points, leaf=None, boxsize=None)
``` 

```
knn_search(tree, target, k, sorted = True)
```

```
ball_search(tree, target, radius, sorted = False)
```

```
box_search(tree, box)
```

## Deallocation of variables (avoid memory leaks)

### Tree

To deallocate the tree in Fortran:

```
call deallocate_kdtree(node)
```

which will deallocate node and all its children.

To deallocate the tree in Python (with the same effect):

```
deallocate_kdtree(tree)
```

### Queries

To deallocate query data in Fortran:

```
call deallocate_query(query)
```

which will deallocate both `dist` and `idx` arrays of the query object.

In Python, query data is automatically handled by the Python/Fortran wrapper after search results are retrieved.




