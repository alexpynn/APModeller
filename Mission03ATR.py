
# =============================================================================
#  Mission Profile
#
# ATR Based on 275Kt Cruise
#
# Predicted fuel burn for this mission is within 1% of the value quoted in the
# ATR fuel efficency manual
# =============================================================================

import matplotlib.pyplot as plt
from APModeller import Aircraft, atmosphere

# Create the atmospere
atm = atmosphere()

# Set HERA Properties
ATRProps = {'TOMass': 22000,
             'Mass': 22000,
             'Weight': 215820.0,
             'WingArea': 61.7,
             'Cd0': 0.027403,
             'Cdi': 0.031,
             'GATOREff': 0.98,
             'PropEff': 0.85,
             'StallSpeed': 47.6}

TurboGenerator = {'PSFC_Constant': None, # Assuming a Constant Power specific fuel consumption (not realistic) 
                  'PSFC_Map': {'PowerFactor': [0.628, -1.2146, 1.5826],
                               'VelocityFactor': [1, 0], #3.0769231
                               'AltitudeFactor': [1],
                               'NormalPower': 2.2e6,
                               'NormalPSFC': 0.3}
                  } # A map constaining PSFC for Airspeed and engine speed

# Creat an Aircraft
ATR = Aircraft(ATRProps, atm, TurboGenerator)

# Start route
ATR.NewMission()

#Initial Taxi
Taxi_Distance = 3500 # meters (m)
Taxi_Velocity = 8 # meters per second (ms-1)
ATR.Taxi(Taxi_Distance, Taxi_Velocity)#

# Take Off
ATR.TakeOff()
 
# Climb
FL1 = 0 # meters (m)
FL2 = 7620 # meters (m)
V1 = 60 # meters per second
V2 = 87.5 # meters per second
Accel = 0.06 # meters per second squared
ClimbTime = 0.3*3600

ATR.Climb(FL1,FL2,V1,V2,Accel,ClimbTime)

ATR.Climb(FL2,FL2,87.5,141,0.066,815)

#Cruise
CruiseTime=1100
ATR.Cruise(FL2,141,CruiseTime)

#Descent
FL3 = 3000
V3 = 123
DescentAngle = 2.5 #deg
ATR.Descent_Slow(FL2,FL3,141,V3,DescentAngle)

#Hold
HoldTime = 120
ATR.Cruise(FL3,V3,HoldTime)


#Descent
V4 = 58
DescentAngle = 4 #deg
ATR.Descent_Slow(FL3,0,V3,V4,DescentAngle)

#Land
StoppingDistance = 1000
ATR.Land(StoppingDistance,V4,Taxi_Velocity)

#Taxi
ATR.Taxi(Taxi_Distance,Taxi_Velocity)
