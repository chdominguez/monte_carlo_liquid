
# Requirements
* A C compiler (CLANG or GCC)

# Compiling
1.  Clone this repository 
`git clone https://github.com/chdominguez/monte_carlo_liquid`
2. Compile with
`make`

# Usage
After compiling, a **bin** folder should appear with the compiled mc.exe executable. Then, run: 

`mc.exe <input>` 

Where *input* is the configuration for the simulation. Two sample lennard.mc and stillinger.mc files is provided in the root of this directory.

# Theoretical background
## Monte Carlo Simulations
Monte Carlo simulations allow to generate molecular ensembles at a given temperature. The main idea is to perform random trial changes in the molecular structre and accept or reject themn using a Boltzmann distribution criterion. Considering two structures, $r^{(i)}$ and $r^{(j)}$ where $r$ is a vector notation referring to a particular arrangement of the molecular structure. At a given temperature T, the probability of observing two different structures is given by: 

$$\dfrac{P(r^{(i)})}{P(r^{(j)})}=exp(-\dfrac{V(r^{(i)})-V(r^{(j)})}{k_bT})$$

Where we can define $V(r^{(i)})$ as the potential energy for that particular structure.

This program satisfies the above relation by changing the structure of the system randomly such that probability for changing from structure $r^{(i)}$ to structure $r^{(j)}$ and for the reverse are equal. 

## Periodic boundary conditions
This software is aimed at liquid simulations. Because liquids have no definite structure on long lenghtscales but only on short distances, the simulation follows the minimum-image convention to implement periodic boundary conditions. The main idea of periodic boundary conditions is to have the system surrounded by an infinite number of replicas of itself in each direction. The lenght of the repeating box must be chosen such that the characteristics of the liquid can be well simulated.

Periodic boundary conditions can easily be implemented in code following these expressions:

$$
\begin{gathered}
\mathbf{r}_{i j}=\mathbf{r}_j-\mathbf{r}_i=\left\lbrace\begin{array}{c}
\Delta x_{i j} \\
\Delta y_{i j} \\
\Delta z_{i j}
\end{array}\right\rbrace=\left\lbrace\begin{array}{c}
x_j-x_i \\
y_j-y_i \\
z_j-z_i
\end{array}\right\rbrace \\
\mathbf{r}_{i j}^{\mathrm{MI}}=\left(\mathbf{r}_j-\mathbf{r}_i\right)^{\mathrm{MI}}=\left\lbrace\begin{array}{c}
\Delta x_{i j}^{\mathrm{MI}} \\
\Delta y_{i j}^{\mathrm{MI}} \\
\Delta z_{i j}^{\mathrm{MI}}
\end{array}\right\rbrace=\left\lbrace\begin{array}{c}
x_j-x_i+n_x L \\
y_j-y_i+n_y L \\
z_j-z_i+n_z L
\end{array}\right\rbrace
\end{gathered}
$$

Then, one has to choose the shortest distance between the regular $r_{ij}$ and the mirror image $r_{ij}^{MI}$. Furthermore, each time an atom is moved "outside the box" one can simply add or substract L (the box lenght) to the position to re-center the given atom.

## Potential Energy Models for Atomic Liquids

This software utilitzes twoclassical methods for accounting the potential energy of the system.

### Lennard-Jones Fluid
This potential is aimed to simulate liquids comprised of rare gas atoms such as argon. Its the first type implemented in the software.

$$
V_{i j}^{\mathrm{LJ}}=4 \epsilon\left\lbrace\left(\frac{\sigma}{r_{i j}}\right)^{12}-\left(\frac{\sigma}{r_{i j}}\right)^6\right\rbrace
$$

### The Stillinger Model
This potential is aimed to model liquid silicon, the second type of simulation implemented in this software.

$$
V^{\text {Stillinger }}=V_2+V_3=\epsilon \sum_{i=1}^{N-1} \sum_{j=i+1}^N f_2\left(\frac{r_{i j}}{\sigma}\right)+\epsilon \sum_{i=1}^{N-2} \sum_{j=i+1}^{N-1} \sum_{k=j+1}^N f_3\left(\frac{\mathbf{r}_{i j}}{\sigma}, \frac{\mathbf{r}_{i k}}{\sigma}, \frac{\mathbf{r}_{j k}}{\sigma}\right)
$$

## The radial distribution function
In order to compare the outcome of simulations with experiment or with other simulations, the radial distribution function can be used. The output of this function is saved in a .gdist file.

$$g(r)=\dfrac{\rho(r)}{\rho}$$

Which in this program is implemented as follows:

$$
g\left(r_i\right)=\frac{N_i}{N_{\text {accum }} \times N / 2 \times 4 \pi\left(\left(r_i+\Delta r\right)^3-r_i^3\right) \rho}
$$

Where $N_i$ is the number of times a given distance between a pair of atoms happens during the simulation, $N_{accum}$ the number of times the list of distances has been updated, $N$ the total number of atoms, $r_i$ the given distance, $\Delta r$ the discretization of the distances and $\rho$ the density.

## Neighbour lists
In order to save computational time, expensive loops are reduced by computing the terms only for the atoms present in a given range. This range is defined as the cutoff plus a distance $r_{skin}$. Every few Monte-Carlo steps, the neighbour list of each atom is updated. Thus, the program has to evaluate the full atom-atom distance fewer times during all of the simulation.

# License
MIT License

Copyright (c) 2022 Christian Domínguez Dalmases

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
