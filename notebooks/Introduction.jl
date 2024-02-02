### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# ╔═╡ 6bd4ad0a-c089-11ee-1e74-01a1b24bbd2c
md""" 
# Castawave -- an interactive guide
###### Written by: Aidan Blaser
---

In this Julia notebook we first break down the methodology behind the Castawave solver and then test a few examples.

### Governing Equation
In the Eulerian frame, irrotational and incompressible flow is governed by the equation

$\quad \mathbf{u} = \mathbf{\nabla}\phi \, , \quad \nabla \cdot \mathbf{u} = 0 \quad \Rightarrow \quad \nabla^2 \phi = 0 \, ,$

where $\phi$ is the velocity potential and $\nabla^2$ the Laplacian operator. This is known as [Laplace's equation](https://en.wikipedia.org/wiki/Laplace%27s_equation) which occurs frequently in nature. Solutions to Laplace's equation are solely determined by boundary conditions.

### Boundary Conditions

Because the fluid cannot penetrate the (potentially infinitely deep) bottom, we require that at $v(y = -h) = \phi_y(y = -h) = 0$. An infinite bottom takes $-h \rightarrow -\infty$. 

The first surface boundary condition, called the kinematic boundary condition, is the mathematical statement that fluid particles at the surface stay at the surface. If we define the surface as $y = \eta(x,t)$, then we require that the surface only change if the fluid velocity pushes it up or down, i.e.

$\frac{d \eta}{d t} = \eta_t + \eta_x \phi_x = \phi_y \Big|_{y = \eta}$

However, because we've introduced a new unknown (the surface $\eta(x,t)$), we need to have another boundary condition to close the system. This is done by invoking Bernoulli's equation, which for unforced, irrotational time-dependent flows is

$\phi_t + \tfrac{1}{2}(\nabla \phi)^2 + g \eta = 0 \Big|_{y = \eta}$
"""

# ╔═╡ 516e6a59-e31c-47d4-a51f-5d30d148d39d


# ╔═╡ Cell order:
# ╟─6bd4ad0a-c089-11ee-1e74-01a1b24bbd2c
# ╠═516e6a59-e31c-47d4-a51f-5d30d148d39d
