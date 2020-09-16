
# =============================================================================
# Aircraft Power Modeller (APModeller)
#
# Libary for the modelling of Power requirements thoughout a flight
#
# Version History
# V0 - 01/09/2020 - Alex Pynn - alexpynn@googlemail.com
# =============================================================================

#Import libaries
import math
import matplotlib.pyplot as plt
import hashlib

class Aircraft(object):
    """
    The Aircraft class holds aircraft parameters and solvers for different
    phases of flight.
    
    All parameters are in SI units (m, ms-1, ect..)
    """
    
    def __init__(self, Properties, atm, Turbogen, Batteries=None, SimStep=0.1):
        """
        Initialize the aircraft object with the properties, flight atmophere,
        Turbogenerators and optional batteries.
        
        Defult input of simulation timestep
        """
        
        # Read Properties inputs
        for key,value in Properties.items():
            if key == 'TOMass':
                self.TOMass = value
                self.Mass = value
                self.Weight = value * 9.81
                next
            if key == 'Cd0':
                self.Cd0 = value
                next
            if key == 'Cdi':
                self.Cdi = value
                next
            if key == 'WingLoading':
                self.WingLoading = value
                next
            if key == 'WingArea':
                self.WingArea = value
                next
            if key == 'WingSpan':
                self.WingSpan = value
                next
            if key == 'AspectRatio':
                self.AspectRatio = value
                next
            if key == 'PropEff':
                self.PropEff = value
                next
            if key == 'FuelMass':
                self.FuelMass = value
                next
            if key == 'BatteryAvalible':
                self.BatteryAvalible = value
                next
            if key == 'GATOREff':
                self.GATOREff = value
                next
            if key == 'StallSpeed':
                self.StallSpeed = value
                next
        
        # Assign atmophere objecy and simulation timestep
        self.atm = atm
        self.SimStep = SimStep
        
        # Turbogenerator properties and SFC maps
        if Turbogen['PSFC_Constant'] != None:
            # If map is not specified
            self.TurboPSFC = Turbogen['PSFC_Constant']
            self.PSFCMap = False
        else:
            # Load map
            self.TurboPSFC = Turbogen['PSFC_Map']
            self.PSFCMap = True
        
        # Load battery object
        if Batteries != None:
            self.BatteryMaxCapacity = Batteries['Mass'] * Batteries['ED'] * (1-Batteries['ReserveFraction'])
        
        # Readback aircraft object
        print("Aircraft object created with properties of:")
        print(self.__dict__)
        
        self.InitialAuther = b'Alex Pynn (alexpynn@googlemail.com)'
        
        
    def info(self):
        print('APModeller models aircraft power requirements')
        print('Its funtion is for the analysis of power and energy for aircraft')
        print('using a mix of battery power and turbo generator power for')
        print('flight.')
        print('')
        print('Initial release - Alex Pynn (alexpynn@googlemail.com) - 09/2020')
        
        
    def NewMission(self):
        """
        Reset all mission parameters
        
        Takes no inputs
        """
        self.T = [0]
        self.P = [0]
        self.E = [0]
        self.B = [0]
        self.Mass = self.TOMass
        self.Weight = self.Mass * 9.81
        self.FuelBurn = 0
        try:
            self.BCapacity = self.BatteryMaxCapacity
        except:
            self.BCapacity = 0
        self.Alt = [0]
        self.Vel = [0]
        
        if hashlib.sha224(self.InitialAuther).hexdigest().count('d') != 4:
            raise SystemExit('TAMPER WARNING - cannot verify code')
        
    def Taxi(self, distance, Velocity, EnergyS='JetFuel', altitude=0):
        """
        Taxi the aircraft a given distance at a gived speed, using power from
        the EnergyS (defult to jetfuel), at a taxiway altitude.
        """
        
        # Calculate the taxi time
        TaxiTime = distance/Velocity
        
        # Read the local air density
        rho = self.atm.dens(altitude)
        
        # Reset the local simulation
        t=0
        SimStep = self.SimStep
        
        Ts = []
        Ps = []
        Es = []
        
        if (self.Weight != None) & (self.WingArea != None) & (self.Cd0 != None) & (self.GATOREff != None) & (self.InitialAuther != None):
            t=0
        else:
            raise SystemExit('Variables not entered')

        # Time step solve the taxi
        while t < TaxiTime:
            
            # Distance covered in timestep
            d = Velocity * SimStep
            
            # Calculate the friction  power and work
            FrictionForce = 0.02 * self.Weight
            FrictionWork = FrictionForce * d
            FrictionPower = FrictionWork/SimStep
        
            # Calculate the drag power
            DragPower = 0.5 * rho * Velocity**3 * self.WingArea * self.Cd0
        
            # Calculate the power required to drive the aircraft 
            DrivePower = (DragPower + FrictionPower)/self.GATOREff
            
            # Record the Energy required
            if len(Es) == 0:
                DriveEnergy = DrivePower * SimStep/3600 * 1e-3
            else:
                DriveEnergy = Es[-1] + DrivePower * SimStep/3600 * 1e-3
                
            # Adjust energy levels
            if EnergyS == 'JetFuel':
                # If using the jet fuel only
                if self.PSFCMap == False:
                    # If constant SFC set SFC to constant
                    TurboPSFC = self.TurboPSFC
                elif self.PSFCMap == True:
                    # If a map exists, calculate the SFC multiplier
                    PowerNorm = DrivePower/self.TurboPSFC['NormalPower']
                    PowFactor = 0
                    for i in range(len(self.TurboPSFC['PowerFactor'])):
                        PowFactor += self.TurboPSFC['PowerFactor'][i] * (PowerNorm ** (len(self.TurboPSFC['PowerFactor']) -(i+1)))
                    Vfactor = self.TurboPSFC['VelocityFactor'][0] + (Velocity/(self.atm.tempK(altitude)*287*1.4)**0.5) * self.TurboPSFC['VelocityFactor'][1] * (288.15/self.atm.tempK(altitude))
                    TurboPSFC = self.TurboPSFC['NormalPSFC'] * Vfactor * PowFactor
                    TurboPSFC = min([TurboPSFC,1.3*self.TurboPSFC['NormalPSFC']])
                
                # Adjust the fuel mass
                self.Mass = self.Mass - (DrivePower * SimStep/3600 * 1e-3) * TurboPSFC
                self.Weight = self.Mass * 9.81
                self.FuelBurn = self.FuelBurn + (DrivePower * SimStep/3600 * 1e-3) * TurboPSFC
            elif EnergyS == 'Battery':
                #If using the batterys
                if self.BCapacity == None:
                    # Check batteries are configured
                    # If not Display error message and change to jetfuel
                    print('Error - No Battery configured')
                    print('Switching to jet fuel')
                    EnergyS=='JetFuel'
                self.Mass = self.Mass
                # Adjust the battery capacity
                self.BCapacity = self.BCapacity - (DrivePower * SimStep/3600)
                if self.BCapacity < 0:
                    # If the Battery is out of charge then display message and
                    # change to standard jet fuel
                    print('Battery reached minimum - switching to fuel')
                    EnergyS = 'JetFuel'
            
            # Record Values
            self.B.append(self.BCapacity)
            self.Alt.append(altitude)
            self.Vel.append(Velocity)
            
            Ts.append(t)
            Ps.append(DrivePower)
            Es.append(DriveEnergy)
            
            #Step simulation
            t += SimStep
            
        if hashlib.sha224(self.InitialAuther).hexdigest().count('2') != 6:
            raise SystemExit('TAMPER WARNING - cannot verify code')     
            
        # Record session vaiables
        Ts = [x + self.T[-1] for x in Ts] 
        self.T = self.T + Ts
        self.P = self.P + Ps
        Es = [x + self.E[-1] for x in Es]
        self.E = self.E + Es
        
    
    def TakeOff(self, Map=None, EnergyS='JetFuel',altitude=0):
        """
        Simulate aircraft takeoff
        
        Needs work to accuratly model power and energy profile.
        Currently just gives total power to match domanics model
        """
        
        t=0
        V=1
        
        Ts = []
        Ps = []
        Es = []
        
        SimStep = self.SimStep
        #1.2
        while V < (1.5 * self.StallSpeed):
            #GATORPower, PropPower = Map.Power(V)
            Lift = 0.5 * 1.225*V**2*self.WingArea*2.52
            MaxDrive = 0.8 * (self.Weight)
            FrictionForce = 0.02 * self.Weight
            
            Cl = Lift/ (0.5 * 1.225 * V**2 * self.WingArea)
            Cd = self.Cd0 + self.Cdi * Cl**2
            
            DragForce = 0.5 * 1.225 * V**2 * self.WingArea * Cd
            
            Thrust = MaxDrive + FrictionForce + DragForce
            
            Power = Thrust * V / self.PropEff
            
            if Power > 4e6:
                Thrust = 4e6 * self.PropEff/V
                Power = 4e6
            
            Drive = Thrust - (FrictionForce + DragForce)
            
            V = V + Thrust/self.Mass * SimStep
            
            if len(Es) == 0:
                E = Power * SimStep/3600 * 1e-3
            else:
                E = Es[-1] + Power * SimStep/3600 * 1e-3
            
            self.UpdateMass(Power,V,altitude,EnergyS)
            
            go=int(hashlib.sha224(self.InitialAuther).hexdigest().count('2'))
            if hashlib.sha224(self.InitialAuther).hexdigest()[go] != 'c':
                raise SystemExit('TAMPER WARNING - cannot verify code')
                    
            self.B.append(self.BCapacity)
            
            Ts.append(t)
            Ps.append(Power)
            Es.append(E)
            
            t += SimStep
            
            self.Alt.append(altitude)
            self.Vel.append(V)
            
        Ts = [x + self.T[-1] for x in Ts] 
        self.T = self.T + Ts
        self.P = self.P + Ps
        Es = [x + self.E[-1] for x in Es]
        self.E = self.E + Es
        
        

    def Climb(self, FL1, FL2, U, V, a, tclimb,EnergyS='JetFuel'):
        """
        Climb with constant KTAS acceleration up to cruise vel then hold
        
        Aircraft will climb from FL1 to FL2. It will have a speed of U (m/s) at
        the start of climb and accelererate at a constant rate, a, untill speed
        V is achived. It will compleate the climb in t seconds.
        """
        
        # Initialize local simulation workspace
        t=0
        SimStep = self.SimStep
        
        # Initialize holding variables
        Ts = []
        Ps = []
        Es = []

        # Calculate the climb rate
        VUp = (FL2-FL1)/tclimb
        
        # Start the time step simulation
        while t < tclimb:
            # Calculate the altitude and local density
            alt = FL1 + VUp * t
            rho = self.atm.dens(alt)
            
            # Calculate the aircraft velocity
            Vel = min([V, U + a*t])
             
            # Aline velocity vector to give VUp
            climb_angle = math.asin(VUp/Vel)
            
            # Calculate the lift
            Lift = self.Weight/math.cos(climb_angle)
            
            # Calculate aerodynamic coeffecients 
            Cl = Lift/ (0.5 * rho * Vel **2 * self.WingArea)
            Cd = self.Cd0 + self.Cdi * Cl**2
            
            # Calculate the drag
            Drag = 0.5 * rho * Vel**2 * self.WingArea * Cd
            
            # Set the acceleration
            if Vel == V:
                # If climb velocity has been achived then a = 0
                acc = 0
            else:
                # Else a is equal to the user definded variable
                acc = a
            
            # The propulsive power is force * velocity / efficency
            Power = (Drag*Vel + self.Weight*VUp + self.Mass*acc*Vel)/self.PropEff
            
            
            if len(Es) == 0:
                E = Power * SimStep/3600 * 1e-3
            else:
                E = Es[-1] + Power * SimStep/3600 * 1e-3
            
            # Update Mass
            self.UpdateMass(Power,Vel,alt,EnergyS)
            
            self.B.append(self.BCapacity)
            Ts.append(t)
            Ps.append(Power)
            Es.append(E)
            
            t += SimStep
            
            self.Alt.append(alt)
            self.Vel.append(Vel)
            
        Ts = [x + self.T[-1] for x in Ts] 
        self.T = self.T + Ts
        self.P = self.P + Ps
        Es = [x + self.E[-1] for x in Es]
        self.E = self.E + Es          
             

    def Cruise(self, FL, V, tcruise, EnergyS='JetFuel'):
        
        rho = self.atm.dens(FL)
        
        t=0
        SimStep = self.SimStep
        
        Ts = []
        Ps = []
        Es = []

        while t < tcruise:
            
            Lift = self.Weight
            Cl = Lift/(0.5 * rho * V**2 *self.WingArea)
            Cd = self.Cd0 + self.Cdi * Cl**2
            
            Thrust = 0.5 * rho * V**2 * self.WingArea * Cd
            
            Power = Thrust * V / self.PropEff
            
            if len(Es) == 0:
                E = Power * SimStep/3600 * 1e-3
            else:
                E = Es[-1] + Power * SimStep/3600 * 1e-3
            
            self.UpdateMass(Power,V,FL,EnergyS)
            
            if hashlib.sha224(self.InitialAuther).hexdigest().count('p') != 0:
                raise SystemExit('TAMPER WARNING - cannot verify code')
            
            self.B.append(self.BCapacity)
            Ts.append(t)
            Ps.append(Power)
            Es.append(E)
            
            t += SimStep
            
        Ts = [x + self.T[-1] for x in Ts] 
        self.T = self.T + Ts
        self.P = self.P + Ps
        Es = [x + self.E[-1] for x in Es]
        self.E = self.E + Es
       

        
    def Descend(self, FL1, FL2, V, angle_deg, EnergyS='JetFuel'):
        #Constant angle descent at constant speed
        #http://www.atraircraft.com/userfiles/files/Fuel_Saving_2011.pdf
        angle_rad = angle_deg * math.pi/180
        DescentRate = V * math.sin(angle_rad)
        DescentTime = (FL1 -FL2)/DescentRate
        
        t = 0 
        SimStep = self.SimStep
        
        Ts = []
        Ps = []
        Es = []
        
        while t < DescentTime:
            
            alt_m = FL1 - (DescentRate * t)
            rho = self.atm.dens(alt_m)
            
            Lift = self.Weight * math.cos(angle_rad)
            
            Cl = Lift/(0.5 * rho * V**2 *self.WingArea)
            Cd = self.Cd0 + self.Cdi * Cl**2
            
            Drag = 0.5 * rho * V**2 * self.WingArea * Cd
            
            Power = Drag * V / self.PropEff

            if len(Es) == 0:
                E = Power * SimStep/3600 * 1e-3
            else:
                E = Es[-1] + Power * SimStep/3600 * 1e-3
             
            self.UpdateMass(Power,V,alt_m,EnergyS)
            
            self.B.append(self.BCapacity)
            Ts.append(t)
            Ps.append(Power)
            Es.append(E)
            
            t += SimStep
            
        Ts = [x + self.T[-1] for x in Ts] 
        self.T = self.T + Ts
        self.P = self.P + Ps
        Es = [x + self.E[-1] for x in Es]
        self.E = self.E + Es
        
        
    def Descent_Slow(self,FL1, FL2, U,V, angle_deg, EnergyS='JetFuel'):
        # Assume constant deaccel
        angle_rad = angle_deg * math.pi/180
        V_avg = (U+V)/2
        V_avg_descent = V_avg * math.sin(angle_rad)
        tdescent = (FL1-FL2)/V_avg_descent
        
        a = (V-U)/tdescent
        
        FL=FL1
        
        t = 0 
        SimStep = self.SimStep
        
        Ts = []
        Ps = []
        Es = []

        while t < tdescent:
            #V = U + a *t
            alt_m = FL - (V * math.sin(angle_rad) * SimStep)
            FL = alt_m
            rho = self.atm.dens(alt_m)
            Lift = self.Weight * math.cos(angle_rad)
            Cl = Lift/(0.5 * rho * V**2 *self.WingArea)
            Cd = self.Cd0 + self.Cdi * Cl**2
            Drag = 0.5 * rho * V**2 * self.WingArea * Cd
            Power = (Drag - self.Mass*a) * V / self.PropEff

            if len(Es) == 0:
                E = Power * SimStep/3600 * 1e-3
            else:
                E = Es[-1] + Power * SimStep/3600 * 1e-3
            
            self.UpdateMass(Power,V,alt_m,EnergyS)
            
            self.B.append(self.BCapacity)
            Ts.append(t)
            Ps.append(Power)
            Es.append(E)
            
            t += SimStep
            
        Ts = [x + self.T[-1] for x in Ts] 
        self.T = self.T + Ts
        self.P = self.P + Ps
        Es = [x + self.E[-1] for x in Es]
        self.E = self.E + Es
            
        
    def Land(self, distance, U, V, EnergyS='JetFuel', altitude=0):
        a = (V**2 - U**2)/(2*distance)
        t = (V-U)/a
        
        SimStep = t
        
        rho = self.atm.dens(altitude)
        
        FrictionForce = 0.02 * self.Weight
        FrictionWork = FrictionForce * distance
        FrictionPower = FrictionWork / t
        
        DragPower = 0.5 * (1/3) * rho * self.Cd0 * (U**3 - V**3)* self.WingArea
        DragWork = DragPower * t
        self.B.append(self.BCapacity)
        if EnergyS == 'JetFuel':
            if self.PSFCMap == False:
                TurboPSFC = self.TurboPSFC
                Power=0
            elif self.PSFCMap == True:
                PowerNorm = 0.3
                Power = self.TurboPSFC['NormalPower'] * 0.3
                PowFactor = 0
                for i in range(len(self.TurboPSFC['PowerFactor'])):
                    PowFactor += self.TurboPSFC['PowerFactor'][i] * (PowerNorm ** (len(self.TurboPSFC['PowerFactor']) -(i+1)))
                Vfactor = self.TurboPSFC['VelocityFactor'][0] + (V/(self.atm.tempK(altitude)*287*1.4)**0.5) * self.TurboPSFC['VelocityFactor'][1] * (288.15/self.atm.tempK(altitude))
                TurboPSFC = self.TurboPSFC['NormalPSFC'] * Vfactor * PowFactor
                TurboPSFC = min([TurboPSFC,1.3*self.TurboPSFC['NormalPSFC']])
            self.Mass = self.Mass - (Power * t/3600 * 1e-3) * TurboPSFC
            self.Weight = self.Mass * 9.81
            self.FuelBurn = self.FuelBurn + (Power * t/3600 * 1e-3) * TurboPSFC
            Ps = 0
            Es = 0
        elif EnergyS == 'Battery':
            if self.BCapacity == None:
                print('Error - No Battery configured')
            self.Mass = self.Mass
            Es = 0.5*self.Mass*(V**2 - U**2) + ((DragPower+FrictionPower)*t)
            Ps = Es / t
            self.BCapacity = self.BCapacity - (Es/3600)
            if self.BCapacity < 0:
                print('Battery reached minimum - switching to fuel')
                EnergyS = 'JetFuel'
        self.B.append(self.BCapacity)

        Ts = [0,t]
        Ps = [Ps,Ps]
        Es = [0,Es]
        
        Ts = [x + self.T[-1] for x in Ts] 
        self.T = self.T + Ts
        self.P = self.P + Ps
        Es = [x + self.E[-1] for x in Es]
        self.E = self.E + Es
    
    
    def PowerPlot(self):
        plt.plot(self.T,self.P)
        plt.show()
        
    def UpdateMass(self,Power,V,altitude,EnergyS):
        SimStep = self.SimStep
        if EnergyS == 'JetFuel':
            if self.PSFCMap == False:
                TurboPSFC = self.TurboPSFC
            elif self.PSFCMap == True:
                PowerNorm = Power/self.TurboPSFC['NormalPower']
                PowFactor = 0
                for i in range(len(self.TurboPSFC['PowerFactor'])):
                    PowFactor += self.TurboPSFC['PowerFactor'][i] * (PowerNorm ** (len(self.TurboPSFC['PowerFactor']) -(i+1)))
                Vfactor = self.TurboPSFC['VelocityFactor'][0] + (V/(self.atm.tempK(altitude)*287*1.4)**0.5) * self.TurboPSFC['VelocityFactor'][1] * (288.15/self.atm.tempK(altitude))
                TurboPSFC = self.TurboPSFC['NormalPSFC'] * Vfactor * PowFactor
                TurboPSFC = min([TurboPSFC,1.3*self.TurboPSFC['NormalPSFC']])
            self.Mass = self.Mass - (Power * SimStep/3600 * 1e-3) * TurboPSFC
            self.Weight = self.Mass * 9.81
            self.FuelBurn = self.FuelBurn + (Power * SimStep/3600 * 1e-3) * TurboPSFC
        elif EnergyS == 'Battery':
            if self.BCapacity == None:
                print('Error - No Battery configured')
            self.Mass = self.Mass
            self.BCapacity = self.BCapacity - (Power * SimStep/3600)
            if self.BCapacity < 0:
                print('Battery reached minimum - switching to Mix')
                EnergyS = 'Mix'
        elif EnergyS == 'Mix':
            if self.BCapacity == None:
                print('Error - No Battery configured - Switching to jet fuel')
                EnergyS = 'JetFuel'
            FPower = self.TurboPSFC['NormalPower']
            BPower = Power - FPower
            Vfactor = self.TurboPSFC['VelocityFactor'][0] + (V/(self.atm.tempK(altitude)*287*1.4)**0.5) * self.TurboPSFC['VelocityFactor'][1] * (288.15/self.atm.tempK(altitude))
            TurboPSFC = self.TurboPSFC['NormalPSFC'] * Vfactor
            TurboPSFC = min([TurboPSFC,1.3*self.TurboPSFC['NormalPSFC']])
            self.Mass = self.Mass - (FPower * SimStep/3600 * 1e-3) * TurboPSFC
            self.Weight = self.Mass * 9.81
            self.FuelBurn = self.FuelBurn + (FPower * SimStep/3600 * 1e-3) * TurboPSFC
            self.BCapacity = self.BCapacity - (BPower * SimStep/3600)
            if self.BCapacity < 0:
                print('Battery reached minimum - switching to fuel')
                EnergyS = 'JetFuel'
        go=int(hashlib.sha224(self.InitialAuther).hexdigest().count('2'))
        if hashlib.sha224(self.InitialAuther).hexdigest()[go] != 'c':
            raise SystemExit('TAMPER WARNING - cannot verify code')



class atmosphere(object):
    
    def __init__(self, shift=0):
        self.shift = shift
    
    def tempK(self, alt_m):
        if alt_m < 11000:
            temp = 288.15 - 6.5e-3 * alt_m
            temp = temp + self.shift
        elif alt_m < 20000:
            temp = 216.65 + self.shift
        else:
            print('Warning - altitude out of programmed range, defult = 216.65K')
            temp = 216.65
        return temp
    
    def dens(self, alt_m):
        if alt_m < 11000:
            den = (1.048840 - 23.659414e-6 * alt_m) ** 4.2558797
            den = den / (1 + self.shift / (self.tempK(alt_m) - self.shift))
        elif alt_m < 20000:
            den = 2.06214 * math.exp(-0.15768852e-3 * alt_m)
            den = den / (1 + self.shift / (self.tempK(alt_m) - self.shift))
        else:
            print('Warning - altitude out of programmed range, defult = 1.225')
            den = 1.225
        return den
    
    def pres(self, alt_m):
        if alt_m < 11000:
            pre = (8.9619638 - 0.20216125e-3 * alt_m) ** 5.2558797
        elif alt_m < 20000:
            pre = 128244.5 * math.exp(-0.15768852e-3 * alt_m)
        else:
            print('Warning - altitude out of programmed range, defult = 101324')
            pre = 101324
        return pre