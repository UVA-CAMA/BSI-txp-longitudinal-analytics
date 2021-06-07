# Longitudinal modeling for transplant Bloodstream Infection
Language: R </br>
Modeling Strategy: </br>
1. Preprocess data: slice 48 hours time windows for pre-defined positive / negative / baseline episodes around blood culture acquisitions.
2. Regress physiological and clinical variables against restricted cubic spline trasnformed time effect using Generalized Least Squares(GLS). 
3. Adjust within-subject dependence (repeated measures) using AR(1) random effect.
4. Interact with 4 primary groups: transplant recipients' positive infection episodes / transplant recipients baseline and negative episodes / no-transplant patients’ positive episodes / no-transplannt patients’ baseline and negative episodes.
5. Statistical inference on subtraction of fitted effects to contrast across groups.
