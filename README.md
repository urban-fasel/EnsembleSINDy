# EnsembleSINDy
 
Ensemble-SINDy: Robust sparse identification of nonlinear dynamics for ordinary- and partial differential equations with uncertainty quantification

Urban Fasel^1, J. Nathan Kutz^2, Bingni W. Brunton^3, Steven L. Brunton^1

^1 Department of Mechanical Engineering, University of Washington, USA

^2 Department of Applied Mathematics, University of Washington, USA

^3 Department of Biology, University of Washington, USA 


Matlab code to reproduce the results: 

## Ensemble-SINDy

### Lorenz system

SINDy/main_runEnsembleSINDy_heatmap.m
 * compare different ensemble SINDy methods in terms of model error (error of SINDy coefficient) and success rate (probability to get correct model structure) over a range of noise levels and data length
 * plot heatmap
    
SINDy/main_runEnsembleSINDy_UQ.m
 * for one particular noise level and data length, run Enesmble-SINDy and get plots of:
    * SINDy coefficient uncertainties and inclusion probabilities
    * ensemble dynamics reconstruction and forecast 

### Lotka Volterra Lynx and Hare polpulation data 1900-1920 Hudson Bay Company

SINDy/main_runLotkaVolterra.m
 * SINDy coefficient uncertainties and inclusion probabilities
 * ensemble dynamics reconstruction


## PDE-FIND

PDE-FIND/main_weakEnsemblePDE.m
 * run comparison between weak PDE-FIND and weak-ensemble PDE-FIND to identify five different PDEs fromvery noisy data

PDE-FIND/main_get_FiguresTable.m
 * plot figures shown in table below 
