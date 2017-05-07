# ISA Calculator

import math

print "What do you want to do?"
print "\n"                                      
print "1)Enter an altitude in meters"
print "2)Enter an altitude in feet"
print "3)Enter pressure in Pa"
print "4)Calculate pressure, temperarure and density given the atmospheric layer and altitude"
print "5)Quit"
print "\n"
print "\n"


c = input("Your choice:")

if c == 1:
    n = input ("Enter an altitude in meters:")
    h = n

elif c == 5:
        quit()

elif c == 2:
    n = input ("Enter an altitude in feet:")
    h = 0.3048*n    

    #feet to meters conversion

# question a

# (if h>11000:
#    print "Values cannot be calculated"
#    quit())

# question a

    if h<=11000:

        print "Value in the required domain. h =", (h)
    
        #collect all constants
        To = (288.15)
        a = (-0.0065)
        g0 = 9.80665
        R = 287.00
        p0 = 101325.0
        rho0 = 1.225

        #calculation of temperature, pressure and density 

        T = (To + a*h)
        print "T =",(T)

        P = p0*((T/To)**(-g0/(a*R)))
        print "P =",(P)

        rho = rho0*((T/To)**((-g0/(a*R))-1))
        print "rho =",(rho)

    elif h>11000 and h<20000:

    #include stratosphere isotherm layer (h>=20km)
    
           #temperature at the start of stratosphere isotherm layer

        #collect all constants
        To = (288.15)
        a = (-0.0065)
        g0 = 9.80665
        R = 287.00
        p0 = 101325.0
        rho0 = 1.225
        T1 = 216.65

    #old calculations
    
        T = (To + a*h)
        P = p0*((T/To)**(-g0/(a*R)))
        rho = rho0*((T/To)**((-g0/(a*R))-1))

    #new calculations
        T2 = T1 
        P1 = P*(math.e)**((-g0/(R*T))*(h-11000))
        rho1 = rho*(math.e)**(((-g0/(R*T))*(h-11000)))

        print (T2)
        print (P1)
        print (rho1)

    elif h>20000:

    #collect all constants
        To = (288.15)
        a1 = (1)
        g0 = 9.80665
        R = 287.00
        p3 = 5474.9
        rho2 = 0.29626
        T2 = 216.65

        T3 = (T2 + a1*(h-20000))
        P4 = p3*((T3/T2)**(-g0/(a1*R)))
        rho3 = rho2*((T3/T2)**((-g0/(a1*R))-1))
    
        print (T3)
        print (P4)
        print (rho3)

elif c == 3:
    Pre = input ("Enter an pressure in Pa:")

    if Pre<22625:
        T1 = (P1/P0)**((-a*R)/g0)         

#collect all constants
        To = (288.15)
        a1 = (1)
        g0 = 9.80665
        R = 287.00
        p3 = 5474.9
        rho2 = 0.29626
        Tlayer = 216.65
       
        print ("Altitude in troposphere, with h = ")

    if Pre<5471.74:

#collect all constants
        To = (288.15)
        a1 = (1)
        g0 = 9.80665
        R = 287.00
        p3 = 5474.9
        rho2 = 0.29626
        Tlayer = 216.65        

        print ("Altitude in isotherm layer of strarosphere, with h = ")


elif c == 4:
    
    def Troposphere(h):
     #temperature gradient#
        a = -0.0065
        
    #collect all constants
        T0 = (288.15)
        g0 = 9.80665
        R = 287.00
        P = 101325.0
        rho0 = 1.225

        T = 288.15 + a*h
        P = P*((T/T0)**(-g0/(a*R)))
        rho = P/(T*R)
        return P,T,rho

    def Tropopause(h):

    #collect all constants
        T0 = (288.15)
        g0 = 9.80665
        R = 287.00
        P = 101325.0
        rho0 = 1.225
        
        P1 = 101325.0*math.exp((9.80665)*(h-11000)/(287*T0))
        rho1 = P/(R*T0)
        return P1, T0, rho1

    def Stratosphere1(h):

    #collect all constants
        T0 = (288.15)
        g0 = 9.80665
        R = 287.00
        P = 101325.0
        rho0 = 1.225
        
        a = 0.001
        T2 = T +a*(h-20000)
        P2 = P1*((T2/T1)**(-g0/(a*R)))
        rho2 = P2/(R*T2)
        return P2, T2, rho2
    
    def Stratosphere2(h):

    #collect all constants
        T0 = (288.15)
        g0 = 9.80665
        R = 287.00
        P = 101325.0
        rho0 = 1.225
        
        a = 0.0028
        T3 = T2 +a*(h-32000)
        P3 = P2*((T3/T2)**(-g0/(a*R)))
        rho3 = P3/(R*T3)
        return P3, T3, rho3

    def Stratopause(h):

    #collect all constants
        T0 = (288.15)
        g0 = 9.80665
        R = 287.00
        P = 101325.0
        rho0 = 1.225
        
        a = 0.00
        T4 = T3
        P4 = P3*math.exp((g0)*(h-47000)/(R*T4))
        rho4 = P4/(R*T4)
        return P4,T4, rho4

    def Mesosphere1(h):

    #collect all constants
        T0 = (288.15)
        g0 = 9.80665
        R = 287.00
        P = 101325.0
        rho0 = 1.225
        
        a = -0.0028
        T5 = T4 +a*(h-51000)
        P5 = P4*((T5/T4)**(-g0/(a*R)))
        rho5 = P5/(R*T5)
        return P5,T5, rho5

    def Mesosphere2(h):

    #collect all constants
        T0 = (288.15)
        g0 = 9.80665
        R = 287.00
        P = 101325.0
        rho0 = 1.225
        
        a = -0.002
        T6 = T5 +a*(h-51000)
        P6 = P5*((T6/T5)**(-g0/(a*R)))
        rho6 = P6/(R*T6)
        return P6,T6, rho6

    def Mesopause(h):

    #collect all constants
        T0 = (288.15)
        g0 = 9.80665
        R = 287.00
        P = 101325.0
        rho0 = 1.225
        
        a = 0.00
        T7 = T6
        P7 = P*math.exp((g0)*(h-84552)/(R*T7))
        rho7 = P7/(R*T7)
        return P7, T7, rho7


    import math

    h = input("Enter altitude in meters: ")

    #base altitudes list

    hbase = [0.0, 11000., 20000., 32000., 47000., 51000., 71000., 84852.]
    htop = [11000., 20000., 32000., 47000., 51000., 71000., 84852., 99999.]
    H = [0, 11000, 20000, 32000, 47000, 51000, 71000, 84852, 99999]

    #temperature gradients list
    a = [-0.0065, 0, 0.001, 0, -0.0028, -0.002, 0]

if h>H[len(H)-1]:
    print "The required altitude is not in the atmosphere"

else:
    for i in range(len(H)):

        if h <= H[i + 1]:
            if i == 0 :
                Pressure, Temperature, Density = Troposphere(h)
                print "Altitude is in the Troposphere. Pressure, Temperature and density at this altitude :"
                print "Pressure (in Pa): ", (Pressure)
                print "Temperature (in K): ", (Temperature)
                print "Density (in kg/m**3): ", (Density)
                break

        elif i == 1:
                Pressure, Temperature, Density = Tropopause(h)
                print "Altitude is in the Tropopause, the isotherm layer of the Troposphere. Pressure, Temperature and density at this altitude :"
                print "Pressure (in Pa): ", (Pressure)
                print "Temperature (in K): ", (Temperature)
                print "Density (in kg/m**3): ", (Density)
                break

        elif i == 2:
                Pressure, Temperature, Density = Stratosphere1(h)
                print "Altitude is in the first layer of the Stratosphere. Pressure, Temperature and density at this altitude :"
                print "Pressure (in Pa): ", (Pressure)
                print "Temperature (in K): ", (Temperature)
                print "Density (in kg/m**3): ", (Density)
                break
            
        elif i == 3:
                Pressure, Temperature, Density = Stratosphere2(h)
                print "Altitude is in the second layer of the Stratosphere. Pressure, Temperature and density at this altitude :"
                print "Pressure (in Pa): ", (Pressure)
                print "Temperature (in K): ", (Temperature)
                print "Density (in kg/m**3): ", (Density)
                break

        elif i == 4:
                Pressure, Temperature, Density = Stratopause(h)
                print "Altitude is in the third and isotherm layer of the Stratosphere, the Stratopause. Pressure, Temperature and density at this altitude :"
                print "Pressure (in Pa): ", (Pressure)
                print "Temperature (in K): ", (Temperature)
                print "Density (in kg/m**3): ", (Density)
                break

        elif i == 5:
                Pressure, Temperature, Density = Mesosphere1(h)
                print "Altitude is in the first layer of the Mesosphere. Pressure, Temperature and density at this altitude :"
                print "Pressure (in Pa): ", (Pressure)
                print "Temperature (in K): ", (Temperature)
                print "Density (in kg/m**3): ", (Density)
                break

        elif i == 6:
                Pressure, Temperature, Density = Mesosphere2(h)
                print "Altitude is in the second layer of the Mesosphere. Pressure, Temperature and density at this altitude :"
                print "Pressure (in Pa): ", (Pressure)
                print "Temperature (in K): ", (Temperature)
                print "Density (in kg/m**3): ", (Density)
                break

        elif i == 7:
                Pressure, Temperature, Density = Stratopause(h)
                print "Altitude is in the third and isotherm layer of the Mesosphere, the Mesopause. Pressure, Temperature and density at this altitude :"
                print "Pressure (in Pa): ", (Pressure)
                print "Temperature (in K): ", (Temperature)
                print "Density (in kg/m**3): ", (Density)
                break 
