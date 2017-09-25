#KRM


# h in km
# T in K
# a in m/s
# rho in kg/m3
# P in pascals
# mu in N-s/m2


function stdatm(h)
    TSL = 288.15 #sea level temp
    PSL = 1.01325e5 #sea level pressure
    mu_sea = 1.79e-5 #sea level viscosity
    g = 9.80665
    R = 287.053
    Sc = 110.4  #Sutherlands formula constant
    y = 1.4 #air ratio of specific heat
    RoM = 8.314/0.0289644 #287.058


    T = TSL - 71.5 + 2.0 *log(1 + exp(35.75 - 3.25*h) + exp(-3.0 + 0.0003*h^3)) #fit
    P = PSL*exp(-0.118*h- 0.0015*h^2 /(1-0.018*h+0.0011*h^2)) #fit
    rho = P/(R*T) #gas constant
    mu = mu_sea*(T/TSL)^(3.0/2)*(TSL+Sc)/(T+Sc) # southerlands

    a=sqrt.(y*RoM*T)
    return rho, P,T,a,mu
end


# fileLoc = splitdir(@__FILE__)
# push!(LOAD_PATH,"$(fileLoc[1])/../hale-module-miscellany/")
# using StandardAtmosphere
# using PyPlot
# # Verification
# N = 1000
# alt = linspace(0,47,N)
# rho = zeros(alt)
# P = zeros(alt)
# T = zeros(alt)
# a= zeros(alt)
# rho2 = zeros(alt)
# P2 = zeros(alt)
# T2 = zeros(alt)
# a2= zeros(alt)
# mu= zeros(alt)
# mu2= zeros(alt)
# for i = 1:N
#     rho[i], P[i],T[i],a[i],mu[i] = stdatm(alt[i])
#     rho2[i], mu2[i],q,a2[i] = atmosphere(alt[i]*1000,"Fit")
# end
#
#
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
# # PyPlot.figure()
# # PyPlot.plot(P,alt,label = "myatm P")
# # PyPlot.plot(P2*1000,alt,label = "judd atm P")
# # PyPlot.legend(loc = "best")
# # PyPlot.ylabel("Altitude (km)")
# # PyPlot.xlabel("P ")
# #
# # PyPlot.figure()
# # PyPlot.plot(T,alt,label = "myatm T")
# # PyPlot.plot(T2,alt,label = "judd atm T")
# # PyPlot.legend(loc = "best")
# # PyPlot.ylabel("Altitude (km)")
# # PyPlot.xlabel("T")
#
# PyPlot.figure()
# PyPlot.plot(mu,alt,label = "myatm T")
# PyPlot.plot(mu2,alt,label = "judd atm T")
# PyPlot.legend(loc = "best")
# PyPlot.ylabel("Altitude (km)")
# PyPlot.xlabel("Mu")
