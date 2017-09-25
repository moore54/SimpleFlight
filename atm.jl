fileLoc = splitdir(@__FILE__)
push!(LOAD_PATH,"$(fileLoc[1])/../hale-module-miscellany/")
using StandardAtmosphere
using PyPlot

#h in km
#T in K

function stdatm(h)
    TSL = 288.15
    PSL = 1.01325e5
    usl = 1.79e-5
    g = 9.80665
    R = 287.053
    T = TSL - 71.5 + 2.0 *log(1 + exp(35.75 - 3.25*h) + exp(-3.0 + 0.0003*h^3))
    P = PSL*exp(-0.118*h- 0.0015*h^2 /(1-0.018*h+0.0011*h^2))
    rho = P/(R*T)

    y = 1.4
    RoM = 8.314/0.0289644 #287.058
    a=sqrt(y*RoM*T)
    return rho, P,T,a
end

N = 1000
alt = linspace(0,47,N)
rho = zeros(alt)
P = zeros(alt)
T = zeros(alt)
a= zeros(alt)
rho2 = zeros(alt)
P2 = zeros(alt)
T2 = zeros(alt)
a2= zeros(alt)
for i = 1:N
    rho[i], P[i],T[i],a[i] = stdatm(alt[i])
    # rho2[i], mu,q,a2[i],P2[i],T2[i] = atmosphere(alt[i]*1000,"Fit")
end


# PyPlot.figure()
# PyPlot.plot(rho,alt,label = "myatm rho")
# PyPlot.plot(rho2,alt,label = "judd atm rho")
# PyPlot.legend(loc = "best")
# PyPlot.ylabel("Altitude (km)")
# PyPlot.xlabel("Rho (K)")
#
# PyPlot.figure()
# PyPlot.plot(a,alt,label = "myatm a")
# PyPlot.plot(a2,alt,label = "judd atm a")
# PyPlot.legend(loc = "best")
# PyPlot.ylabel("Altitude (km)")
# PyPlot.xlabel("a (m/s)")
#
# PyPlot.figure()
# PyPlot.plot(P,alt,label = "myatm P")
# PyPlot.plot(P2*1000,alt,label = "judd atm P")
# PyPlot.legend(loc = "best")
# PyPlot.ylabel("Altitude (km)")
# PyPlot.xlabel("P ")
#
# PyPlot.figure()
# PyPlot.plot(T,alt,label = "myatm T")
# PyPlot.plot(T2,alt,label = "judd atm T")
# PyPlot.legend(loc = "best")
# PyPlot.ylabel("Altitude (km)")
# PyPlot.xlabel("T")

#HW1.end
b=61 #m
Sref = 325.16 #m2
mass=200000 #kg
weight = mass*9.81
M=0.82
V = M*340.0 #m/s
Cd0cl=0.01
e=0.7

AR = b^2/Sref

function lod(rho,M,weight,V,Sref,e,AR,Cd0cl)
    Cl = 2*weight/(rho*V^2*Sref)#/sqrt(1-M^2)
    Cd = Cl^2/(pi*e*AR)+Cd0cl
    return Cl,Cd
end

PyPlot.figure()
N = 500
alt = linspace(1,25,N) #meters
NN = 3
VL0D = zeros(length(alt),NN)
V = zeros(length(alt),NN)
Lsave= zeros(length(alt),NN)
Dsave = zeros(length(alt),NN)
i = 0
j = 0

weightsweep = linspace(weight,0.6*weight,NN)
weightsweep = [weight,weight-9.81*10000*2,weight-9.81*10000*4]
masssweep = round(Int,weightsweep/9.81)
for j = 1:NN
    for i = 1:N
        rho, P,T,a = stdatm(alt[i])
        # M = V/a # assume constant mach
        V[i,j] = M*a
        Cl,Cd = lod(rho,M,weightsweep[j],V[i,j],Sref,e,AR,Cd0cl)
        L = Cl/2*rho*V[i,j]^2*Sref
        Lsave[i,j] = Cl
        D = Cd/2*rho*V[i,j]^2*Sref
        Dsave[i,j] = Cd
        VL0D[i,j] = V[i,j]*L/D
    end

    PyPlot.plot(VL0D[:,j],alt,label = "Mass (kg) = $(masssweep[j])")

end

PyPlot.legend(loc = "best")
PyPlot.ylabel("Altitude (km)")
PyPlot.xlabel("V*L/D")

#
# PyPlot.figure()
# for j = 1:NN
# PyPlot.plot(V[:,j],alt,label = "Mass (kg) = $(masssweep[j])")
# end
# PyPlot.legend(loc = "best")
# PyPlot.ylabel("Altitude (km)")
# PyPlot.xlabel("V (m/s)")
#
# PyPlot.figure()
# for j = 1:NN
# PyPlot.plot(Lsave[:,j],alt,label = "Mass (kg) = $(masssweep[j])")
# end
# PyPlot.legend(loc = "best")
# PyPlot.ylabel("Altitude (km)")
# PyPlot.xlabel("CL ")
#
# PyPlot.figure()
# for j = 1:NN
# PyPlot.plot(Dsave[:,j],alt,label = "Mass (kg) = $(masssweep[j])")
# end
# PyPlot.legend(loc = "best")
# PyPlot.ylabel("Altitude (km)")
# PyPlot.xlabel("CD")
#
# PyPlot.figure()
# for j = 1:NN
# PyPlot.plot(Lsave[:,j]./Dsave[:,j],alt,label = "Mass (kg) = $(masssweep[j])")
# end
# PyPlot.legend(loc = "best")
# PyPlot.ylabel("Altitude (km)")
# PyPlot.xlabel("L/D")
#
# PyPlot.figure()
# for j = 1:NN
# PyPlot.plot(Lsave[:,j],alt,".-",label = "CL, Mass (kg) = $(masssweep[j])")
# PyPlot.plot(Dsave[:,j],alt,"*-",label = "CD, Mass (kg) = $(masssweep[j])")
# # PyPlot.plot(V[:,j]/1000,alt,label = "V(km/s) (Mass (kg) = $(masssweep[j])")
# end
# PyPlot.legend(loc = "best")
# PyPlot.ylabel("Altitude (km)")
# PyPlot.xlabel("Coefficient")
