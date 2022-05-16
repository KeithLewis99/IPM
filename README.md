# IPM
This project contains the code and data required to develop an Integrated Population Model (IPM) in support of a Limit Reference Point (LRP) for 2J3KL capelin (Mallotus villosus).  

Integrated Population models are defined as a model with a state-space model as a sub-component.....see Schuab and Kerry (pg).

Data
For complete metadata and details on the various data components, see the 2K3KL Data Dictionary (in development).  The current iteration of the model contains 
1. the age-disagregated data from the spring acoustic survey (ages 2-4) from 1999-2019,
2. capelin condition from 1999 to present
3. the percent mature capelin from 1999 to present.  This is a fixed amount for ages 3 and 4 but time varying for age 2.
3. larval density from Bellevue Beach 2001 to present.
4. Timing of ice retreat from 1999 to present

Model information
The current iteration also utilizes the capelin forecast model (Lewis et al. 2019) to model abundance of age 2 fish.  At each age, it is assumed that mature fish migrate to the inshore, spawn and die (may relax this last assumption later).  The immature fish stay offshore until the next year.  The model without larval density, i.e., timing of ice retreat and condition, then "kills off" fish between ages 2 and 3 as well as ages 3 and 4.  These fish are then available to the spring acoustic index.

There are three models in the current iteration.
1. A model with an independent parameter for each of the variables for each age
2. A model with an independent parameter for each of the variables for age-2 but the same parameters for age-3 & age-4
3. A "null" model that focuses only on the demographics of age-3 and age-4.

Project Information
The current project is organized in two streams, both starting with the data file, IPM_dat.R.  The first branch is for exploratory data analysis (EDA).  The file IPM_EDA creates figures from the data and these are presented in a dashboard.  The second branch is for the IPM.  The IPM_out.R file sources the data file, the model file with the JAGS code (IPM_mod.R), some Zuur code for the model diagnostics.  These are then presented in a dashboard.

Finally, there is a "scratch pad" where I develop new code and try things out.


