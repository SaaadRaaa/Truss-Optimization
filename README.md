# Truss-Optimization
Analysis and Optimization of Truss Structures using Genetic Algorithms

This repository contains a Python program that optimizes truss designs for maximum strength while minimizing the total weight. The program uses genetic algorithms to find the best truss design given a set of constraints.

## Table of Contents
* Prerequisites
* Installation
* Usage
* Inputs
* Outputs
* Contributing

## Prerequisites
Before running this program, you will need to have the following software installed on your computer:
* Python 3
* NumPy
* SciPy

## Installation
To use this program, you will need to have Python 3 installed on your computer. You can download Python from the official website: https://www.python.org/downloads/

Once you have Python installed, you can clone this repository to your computer using Git:

```bash
git clone https://github.com/SaaadRaaa/Truss-Optimization.git
```

To install the necessary packages, you can use pip:

```terminal
pip install numpy scipy
```

## Usage
To run the program, navigate to the directory where the repository is located and run the following command:

```terminal
python main.py
```
By default, the program will generate a truss structure with 10 nodes and optimize it for maximum strength while minimizing the total weight. You can customize the number of nodes and other parameters in the program by modifying the values at the loadcase file. You will also need to specify the load and support conditions.

The program will then use genetic algorithms to optimize the truss design for minimum weight, subject to the given constraints.

## Inputs
The program takes the following inputs:

* Node coordinates
* Element connectivity
* Load conditions (magnitude and direction)
* Support/Boundary conditions
* Element properties (density, cross-section area and elastic modulus)

## Outputs
The program outputs the following:

* Optimal truss design (member cross-section areas)
* Total weight of the truss
* Visualization of the truss design

## Contributing
Contributions are welcome! If you have any ideas for improvements or new features, feel free to submit a pull request.
