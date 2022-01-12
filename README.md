# Ensemble SINDy and PDE-FIND
 
This repository contains Matlab code to reproduce the results of the paper:
  
**Ensemble-SINDy: Robust sparse model discovery in the low-data, high-noise limit, with active learning and control**

Urban Fasel, J. Nathan Kutz, Bingni W. Brunton, Steven L. Brunton.  

https://arxiv.org/pdf/2111.10992.pdf


## Ensemble-SINDy

### Lorenz system: Comparing different ensemble-SINDy methods 

[SINDy/main_runEnsembleSINDy_heatmap.m](/SINDY/main_runEnsembleSINDy_heatmap.m)  
 * compare different ensemble SINDy methods in terms of model error (error of SINDy coefficients) and success rate (rate of identifying the correct non-zero and zero terms in the library) over a range of noise levels and data length
 * plot heatmap
    
![ensembleSINDy_heatmap](/SINDY/results/ensembleSINDy_heatmap.png)  


### Lorenz system: UQ and ensemble forecasting
    
[SINDy/main_runEnsembleSINDy_UQ.m](/SINDY/main_runEnsembleSINDy_UQ.m)   
 * for one particular noise level and data length, run Enesmble-SINDy and get plots of:
    * SINDy coefficient uncertainties and inclusion probabilities
    * ensemble dynamics reconstruction and forecast 

### Lotka Volterra: Lynx and Hare polpulation data 1900-1920 Hudson Bay Company

[SINDy/main_runLotkaVolterra.m](/SINDY/main_runLotkaVolterra.m)  
 * SINDy coefficient uncertainties and inclusion probabilities
 * ensemble dynamics reconstruction


## Ensemble-PDE-FIND

Here, we use the weak formulation [code](https://github.com/dm973/WSINDy_PDE) from [Messenger and Bortz](https://arxiv.org/abs/2007.02848) as a baseline.

[PDE-FIND/main_weakEnsemblePDE.m](/PDE-FIND/main_weakEnsemblePDE.m)  
 * run comparison (model error and success rate) between weak PDE-FIND and weak-ensemble PDE-FIND to identify five different PDEs from very noisy data

[PDE-FIND/main_get_FiguresTable.m](/PDE-FIND/main_get_FiguresTable.m)  
 * plot figures shown in table below 

[PDE-FIND/main_weakEnsemblePDE_success09.m](/PDE-FIND/main_weakEnsemblePDE_success09.m)  
 * compare noise levels for different PDEs where weak PDE-FIND and weak-ensemble PDE-FIND success rate drops below 0.90.

![table_ensemble-PDEFIND](/PDE-FIND/PlotsPaper/table_ensemble-PDEFIND.png)


## Active E-SINDy

[Active_and_control/Active/main_activePlot1and2.m](/Active_and_control/Active/main_activePlot1and2.m)  
 * plots active E-SINDy left and center figure 

[Active_and_control/Active/main_activePlot3.m](/Active_and_control/Active/main_activePlot3.m)  
 * plots active E-SINDy right figure 

![ensembleSINDyActive](/Active_and_control/Active/results/ensembleSINDyActive.png)  


## MPC E-SINDy

[Active_and_control/MPC/main_MPC.m](/Active_and_control/MPC/main_MPC.m)  
 * compare SINDy with E-SINDy for model predictive control

![ensembleSINDyMPC](/Active_and_control/MPC/Results/ensembleSINDyMPC.png)  

