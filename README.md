# 🔥 Finite Element Method for the Heat Equation | MATLAB Project

## 📌 Description

This project provides a **numerical solution of the heat equation** in both stationary and time-dependent settings using the **Finite Element Method (FEM)** in MATLAB.

It focuses on modeling heat propagation through different spatial domains over time and includes:
- Dirichlet Boundary Conditions
- Neumann Boundary Conditions
- Fourier (Robin) Boundary Conditions

This project was conducted as part of an academic module on numerical analysis and partial differential equations.

---

## 🚀 Features

✅ Implementation of FEM for 1D and 2D heat propagation  
✅ Treatment of multiple types of boundary conditions  
✅ Time integration using implicit schemes  
✅ Modular and reusable MATLAB code  
✅ Visualizations for temperature evolution over space and time  

---

## 🧠 Theory Overview

The **heat equation** in its time-dependent form is:

\[
$\frac{\partial u}{\partial t} - \nabla \cdot (k \nabla u) = f \quad \text{in } \Omega \times (0, T)$
\]


Where:
- \( u \) is the temperature
- \( k \) is the thermal conductivity
- \( f \) is the heat source
- $\( \Omega \)$ is the spatial domain

The problem is solved using:
- **Galerkin Finite Element Method**
- **Mass and stiffness matrices construction**
- **Time discretization** using backward Euler

---

## 📂 Project Structure

```
/finite_element_method_project
│── principal_chaleur.m      # Main simulation program for heat equation
│── domaine.geo              # Geometry definition file (1)
│── domaine.msh              # Mesh file generated from geometry (1)
│── geomRectangle.geo        # Geometry definition file (2)
│── geomRectangle.msh        # Mesh file generated from geometry (2)
│── lecture_msh.m            # Mesh reading and processing function
│── affichemaillage.m        # Mesh visualization function
│── matK_elem.m              # Elementary stiffness matrix computation
│── matM_elem.m              # Elementary mass matrix computation
│── mat_elem_surface.m       # Surface element matrix computation
│── sigma_1.m                # Thermal conductivity function (domain 1)
│── sigma_2.m                # Thermal conductivity function (domain 2)
│── condition_initiale.m     # Initial condition definition
│── f.m                      # Heat source function
│── f_t.m                    # Time-dependent heat source function
│── T_Gamma.m                # Boundary temperature function
│── elimine.m                # Function to apply Dirichlet conditions
│── affiche.m                # Results visualization function
│── Compte_rendu_TP2_ANN201.pdf  # Detailed project report
```
