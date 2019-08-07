#!/bin/bash

rsync --delete -rvu --progress --exclude={"GrIS/Atmosphere_Forcing/aSMB_observed/v0","GrIS/Ocean_Forcing/Melt_Implementation/v0","GrIS/Ocean_Forcing/Melt_Implementation/v1","GrIS/Ocean_Forcing/Melt_Implementation/v2","GrIS/Ocean_Forcing/Melt_Implementation/v3","GrIS/Ocean_Forcing/Retreat_Implementation","GrIS/Ocean_Forcing/Retreat_Processing","GrIS/output"} aaschwanden@transfer.ccr.buffalo.edu:/projects/grid/ghub/ISMIP6/Projections/GrIS $1
