# Metabolic Pathway Simulation: Allosteric Feedback Inhibition

This repository contains a Python-based kinetic simulation of a 4-reaction metabolic pathway. The project was developed to model the dynamic behavior of a microbial host optimized for high-value biochemical production, specifically focusing on the impact of allosteric feedback inhibition and byproduct diversion.

## 📌 Overview

The simulation models a structural topology where:
- **Substrate X** is converted to **Intermediate A** (Commitment Step $v_1$).
- **Intermediate A** can either proceed to **Product P** via **Intermediate B** ($v_2, v_3$) or be diverted to an unwanted **Byproduct** ($v_4$).
- **Product P** exerts non-competitive allosteric inhibition on the first reaction ($v_1$).

## 🧬 Mathematical Model

The system is governed by the following Ordinary Differential Equations (ODEs):

$$ \frac{d[A]}{dt} = v_1 - v_2 - v_4 $$
$$ \frac{d[B]}{dt} = v_2 - v_3 $$
$$ \frac{d[P]}{dt} = v_3 $$

Where the reaction rates are defined as:
- $v_1 = \frac{V_{max} \cdot [X]}{(K_{m1} + [X]) \cdot (1 + \frac{[P]}{K_i})}$ (Inhibited Michaelis-Menten)
- $v_2 = k_2 \cdot [A]$
- $v_3 = k_3 \cdot [B]$
- $v_4 = k_4 \cdot [A]$

## 🚀 Installation & Usage

1. [Hasil](Hasil%20Programming.png)
2. [Script](Script%20Programming.py)
