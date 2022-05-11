# IPM
This project contains the code and data required to develop an Integrated Population model in support of a Limit Reference Point (LRP) for 2J3KL capelin (Mallotus villosus).  

Integrated Population models are defined as a model with a state-space model as a sub-component.....see Schuab and Kerry (pg).

For complete metadata and details on the various data components, see the Data Dictionary (in development).  The current interation contains 
1. the age-disagregated data from the spring acoustic survey (ages 2-4) from 1999-2019,
2. capelin condition from 1999 to present
3. the percent mature capelin from 1999 to present.  This is a fixed amount for ages 3 and 4 but time varying for age 2.
3. larval density from Bellevue Beach 2001 to present.
4. Timing of ice retreat from 1999 to present

The current interation also utilizes the capelin forecast model (Lewis et al. 2019) to model abundance of age 2 fish.  At each age, it is assumed that mature fish migrate to the inshore, spawn and die (may relax this last assumption later).  The immature fish stay offshore until the next year.  The model without larval density, i.e., timing of ice retreat and condition, then "kills off" fish between ages 2 and 3 as well as ages 3 and 4.  These fish are then available to the spring acoustic index.


