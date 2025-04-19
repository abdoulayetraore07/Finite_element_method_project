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
\frac{\partial u}{\partial t} - \nabla \cdot (k \nabla u) = f \quad \text{in } \Omega \times (0, T)
\]

Where:
- \( u \) is the temperature
- \( k \) is the thermal conductivity
- \( f \) is the heat source
- \( \Omega \) is the spatial domain

The problem is solved using:
- **Galerkin Finite Element Method**
- **Mass and stiffness matrices construction**
- **Time discretization** using backward Euler

---

## 📂 Project Structure

