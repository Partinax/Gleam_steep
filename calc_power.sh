#!/bin/bash

# Use the stilts built-in functions to turn a GLEAM flux density into a 325MHz and 1.4GHz radio power

z=0.11
flux=0.5 #Jy
freq=200. #MHz
alpha=-1.5

DL=`stilts calc "luminosityDistance($z,73,0.27,0.73)"`

flux325=`stilts calc "${flux}*pow((325./$freq),$alpha)"`
flux1400=`stilts calc "${flux}*pow((1400./$freq),$alpha)"`
P325=`stilts calc "pow(10,-26)*fluxToLuminosity($flux325,MpcToM($DL))"`
P1400=`stilts calc "pow(10,-26)*fluxToLuminosity($flux1400,MpcToM($DL))"`

echo $DL $flux325 $flux1400 $P325 $P1400


