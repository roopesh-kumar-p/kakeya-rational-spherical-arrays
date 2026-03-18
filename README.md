# kakeya-rational-spherical-arrays

MATLAB code for Kakeya-inspired rational spherical arrays (beamforming, sum rate capacity, mutual coupling, TTFC analysis)







\# Kakeya-Inspired Rational Spherical Arrays (Reference Research Code)



\## Title

\*\*Kakeya-Inspired Rational Spherical Arrays: Achieving the Lower Bound for Full-Sphere Coverage and Capacity Optimization\*\*



\## Authors

\- Roopesh Kumar Polaganga  

\- Qilian Liang  



\## Affiliation

Department of Electrical Engineering  

The University of Texas at Arlington  



\## Corresponding Author

Roopesh Kumar Polaganga  

Email: RoopeshKumar.Polaganga@mavs.uta.edu  



\## Repository / DOI

https://github.com/roopesh-kumar-p/kakeya-rational-spherical-arrays



\---



\## Description

This repository provides MATLAB code to reproduce the main results presented in the manuscript, including:



\- Integer spherical array baseline (Fig-2)  

\- Rational-KLB array construction (Algorithm-1)  

\- Array factor comparison (Fig-4)  

\- Beam quality metrics: PSLL, ISLR, Directivity (Fig-5)  

\- Sum-rate capacity comparison (Fig-6)  

\- Mutual coupling analysis (Fig-7 / Fig-8)  

\- TTFC (Time-To-Full-Coverage) evaluation (Fig-9)  



\---



\## Algorithm Mapping

The implementation follows \*\*Algorithm-1\*\* in the manuscript:



\- Step 1–9   : Integer array + directivity → N\_KLB estimation  

\- Step 10–17 : Coprime selection (A, B) and rational spherical construction  

\- Step 18–25 : Beamforming and array factor evaluation  

\- Step 26–29 : Coverage and TTFC computation  



\---



\## Usage

Run the script directly in MATLAB:



```matlab

main

