using PyPlot
include("atm.jl")

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
