using Roots
#Function to calculate entire wing parasitic drag

#MAC mean aerodynamic chord
include("atm.jl")

function cf_calc(xcrit,Vinf,mu,rho)


    Re=rho*Vinf*xcrit/mu
    Re_lam = Re*xcrit
    cf_lam = 1.328/sqrt.(Re_lam) # Laminar
    cf_turb = 0.074/(Re^0.2) # Turb

    cf = xcrit*cf_lam+(1-xcrit)*cf_turb

    return cf,Re
end

function wingparasitic(S,Sref,Vinf,sweep_d,MAC,tc,h=1.2,Recrit=500000) #h = 1.2km or 4000 feet


    rho, P,T,a,mu = stdatm(h)

    sweep = sweep_d*pi/180
    function zeroxcrit(x)
        Re = rho*Vinf*x*MAC/mu
        zero = Re-Recrit
        return zero
    end
    xcrit = 1
    try
        xcrit = fzero(0,1,zeroxcrit)

    catch
        xcrit = 1
    end

    cf,Re = cf_calc(xcrit,Vinf,mu,rho)

    k = 1+2*cos(sweep)*tc+100*(tc)^4

    Swet = 2*(1+0.2*tc)*S

    CDp = k*cf*Swet/Sref

    return CDp,Re
end

function bluntbodydrag(Vinf,S_xc,Sref,l,h=1.2,Recrit=500000)

    rho, P,T,a,mu = stdatm(h)

    function zeroxcrit_blunt(x)
        Re = rho*Vinf*x/mu
        zero = Re-Recrit
        return zero
    end
    xcrit = 1
    try
        xcrit = fzero(0,l,zeroxcrit_blunt)
        xcrit = xcrit/l
    catch
        xcrit = 1
    end

    cf,Re = cf_calc(xcrit,Vinf,mu,rho)

    #Noncircular cross section
    deff = sqrt.(4*S_xc/pi)

    #Ellipsoid area
    a = l/2
    b = deff/2
    c=b
    Swet = 4*pi*(((a*b)^1.6+(a*c)^1.6+(b*c)^1.6)/3)^(1/1.6)

    #Form factor body of revolution
    fr = l/deff
    if fr<15
        k = 1.675-0.09*fr+0.003*fr^2
    else
        k = 1
    end

    #Fuselage
    CDpb = k*cf*Swet/Sref

    return CDpb
end

function liftinduced(CL,e,AR,fuselagewidth=0,span=1)
    #Function to calculate entire wing induced drag
    # Fuselage e effect: e = e*(1-2*(fuselagewidth/span)*2)
    if fuselagewidth!=0
        e = e*(1-2*(fuselagewidth/span)^2)
    end

    CDi = CL^2/(pi*e*AR)

    return CDi
end

function viscousliftinduced(CDp,CL,K = 0.38, einv=0.98, AR = 10)
    # Viscous induced
    CDiv = K*CDp*CL^2
    oswald = 1/(1/einv+K*CDp*pi*AR)
    return CDiv, oswald
end



# #Test
# using PyPlot
# span = 1.0 #meter
# AR = 8
# MAC =span/AR
# S = span*MAC
# Sref = S
# sweep_d = 5.0
# fuselagewidth = .15
# tc = .1
# l = 5.0
# S_xc = fuselagewidth^2
# Vinf = 5.0
# N = 100
# CL_a = linspace(-1.5,1.5,N)
#
# CDp_a = zeros(CL_a)
# CDpb_a = zeros(CL_a)
# CDi_a = zeros(CL_a)
# CDiv_a = zeros(CL_a)
# Re = 0
# for i = 1:N
#     CL = CL_a[i]
#     CDp,Re = wingparasitic(S,Sref,Vinf,sweep_d,MAC,tc,1.2,500000)
#     CDpb = bluntbodydrag(Vinf,S_xc,Sref,l,1.2,500000)
#     CDi = liftinduced(CL,e,AR,fuselagewidth,span)
#     CDiv = viscousliftinduced(CDp,CL)
#
#     CDp_a[i] = CDp
#     CDpb_a[i] = CDpb
#     CDi_a[i] = CDi
#     CDiv_a[i] = CDiv
# end
#
# CD = CDp_a + CDi_a + CDiv_a + + CDpb_a
# LoD = CL_a./CD
# PyPlot.figure()
# PyPlot.plot(CL_a,CD,"k-",label = "Total")
# PyPlot.plot(CL_a,CDi_a+CDiv_a,"r-",label = "Induced")
# PyPlot.plot(CL_a,CDp_a,"b-",label = "Wing Viscous")
# PyPlot.plot(CL_a,CDpb_a,"g-",label = "Fuselage Drag")
# PyPlot.xlabel("CL")
# PyPlot.ylabel("CD")
# PyPlot.legend()
#
# PyPlot.figure()
# PyPlot.plot(CL_a,LoD,"k-")
# PyPlot.xlabel("CL")
# PyPlot.ylabel("L/D")
# # PyPlot.legend()
