#Test
using PyPlot
include("Lift.jl")
include("Drag.jl")
include("atm.jl")
span = 1.0 #meter
AR = 12
MAC =span/AR
S = span*MAC
Sref = S
sweep_d = 5.0
fuselagewidth = 0 #0.15
tc = .1
l = .75
S_xc = fuselagewidth^2
N = 500
weight = 0.56*9.81 #Newtons
Vinf = linspace(5,30,N)
einv = 0.95

CL_a = zeros(N)
CDp_a = zeros(CL_a)
CDpb_a = zeros(CL_a)
CDi_a = zeros(CL_a)
CDiv_a = zeros(CL_a)
Re = zeros(CL_a)
q = zeros(CL_a)
oswald = zeros(CL_a)
for i = 1:N
    CL,q[i] = lift_coefficient(weight,Sref,Vinf[i],1.2)
    CDp,Re[i] = wingparasitic(S,Sref,Vinf[i],sweep_d,MAC,tc,1.2,500000)
    CDpb = bluntbodydrag(Vinf[i],S_xc,Sref,l,1.2,500000)
    CDi = liftinduced(CL,einv,AR,fuselagewidth,span)
    CDiv,oswald[i] = viscousliftinduced(CDp,CL, 0.38, einv,  AR)

    CDp_a[i] = CDp
    CDpb_a[i] = CDpb
    CDi_a[i] = CDi
    CDiv_a[i] = CDiv
    CL_a[i] = CL
end

CD = CDp_a + CDi_a + CDiv_a + CDpb_a
CDwing = CDp_a + CDi_a + CDiv_a
LoD = CL_a./CD
LoDwing = CL_a./CDwing
LoDmax,LoDidxMax = findmax(LoD)
LoDwingmax,LoDwingidxMax = findmax(LoDwing)
idx2 = 0
oswald_av = (oswald[LoDwingidxMax+idx2])
CDp_av = (CDp_a[LoDwingidxMax+idx2])
maxLoDanalytical = sqrt.(CDp_av*pi*AR.*oswald_av)-CL_a
LoDanalytical = CL_a.*pi.*oswald_av*AR./(CL_a.^2+CDp_av.*pi.*oswald_av*AR)
LoDanalmax,LoDanalidxMax = findmax(LoDanalytical)
tol = 1e-3
i=0
for i = 1:length(maxLoDanalytical)
    if maxLoDanalytical[i]>0-tol && maxLoDanalytical[i]<0+tol
        break
    end
end
PyPlot.figure()
PyPlot.plot(Vinf,maxLoDanalytical,"k-",label = "Total CD")
PyPlot.xlabel("Vinf")
PyPlot.ylabel("sqrt(CDp_a*pi*AR.*oswald)-CL_a")
PyPlot.legend()

PyPlot.figure()
PyPlot.plot(Vinf,CD,"k-",label = "Total CD")
PyPlot.plot(Vinf,CDi_a+CDiv_a,"r-",label = "Induced")
PyPlot.plot(Vinf,CDp_a,"b-",label = "Wing Viscous")
PyPlot.plot(Vinf,CDpb_a,"g-",label = "Fuselage Viscous")
PyPlot.xlabel("Vinf")
PyPlot.ylabel("CD")
PyPlot.legend()

PyPlot.figure()
PyPlot.plot(Vinf,(CD).*q.*Sref,"k-",label = "Total Drag")
# PyPlot.plot(Vinf,(CDwing).*q.*Sref,"k--",label = "Total Wing Only Drag")
PyPlot.plot(Vinf,(CDi_a+CDiv_a).*q.*Sref,"r-",label = "Induced")
PyPlot.plot(Vinf,(CDp_a).*q.*Sref,"b-",label = "Wing Viscous")
PyPlot.plot(Vinf,(CDpb_a).*q.*Sref,"g-",label = "Fuselage Viscous")
PyPlot.xlabel("Vinf")
PyPlot.ylabel("Drag (Newtons)")
PyPlot.legend()

PyPlot.figure()
# PyPlot.plot(Vinf,(CD).*q.*Sref,"k-",label = "Total Drag")
PyPlot.plot(Vinf,(CDwing).*q.*Sref,"k--",label = "Total Wing Only Drag")
PyPlot.plot(Vinf,(CDi_a+CDiv_a).*q.*Sref,"r-",label = "Induced")
PyPlot.plot(Vinf,(CDp_a).*q.*Sref,"b-",label = "Wing Viscous")
# PyPlot.plot(Vinf,(CDpb_a).*q.*Sref,"g-",label = "Fuselage Viscous")
PyPlot.xlabel("Vinf")
PyPlot.ylabel("Drag (Newtons)")
# PyPlot.xlim([5,25])
# PyPlot.ylim([0,0.5])
PyPlot.legend()

PyPlot.figure()
PyPlot.plot(Vinf,LoD,"k-",label = "Fuselage Included")
PyPlot.plot(Vinf[LoDidxMax],LoDmax,"kx",label = "Max L/D: $(round(LoDmax,2)) at $(round(Vinf[LoDidxMax],2)) m/s
Re = $(round(Re[LoDidxMax],-3))
CL = $(round(CL_a[LoDidxMax],1))")
# PyPlot.xlabel("Vinf")
# PyPlot.ylabel("L/D")
# PyPlot.legend()
#
# PyPlot.figure()
PyPlot.plot(Vinf,LoDwing,"b-", label = "Wing Only")
PyPlot.plot(Vinf[LoDwingidxMax],LoDwingmax,"bx",label = "Max L/D: $(round(LoDwingmax,2)) at $(round(Vinf[LoDwingidxMax],2)) m/s
Re = $(round(Re[LoDwingidxMax],-3))
CL = $(round(CL_a[LoDwingidxMax],1))")
PyPlot.xlabel("Vinf")
PyPlot.ylabel("L/D")
PyPlot.legend(loc = 4,ncol = 2)

PyPlot.figure()
# PyPlot.plot(Vinf,LoD,"k-",label = "Fuselage Included")
# PyPlot.plot(Vinf[LoDidxMax],LoDmax,"kx",label = "Max L/D: $(round(LoDmax,2)) at $(round(Vinf[LoDidxMax],2)) m/s
# Re = $(round(Re[LoDidxMax],-3))
# CL = $(round(CL_a[LoDidxMax],1))")
PyPlot.plot(Vinf,LoDwing,"b-", label = "Wing Only")
PyPlot.plot(Vinf[LoDwingidxMax],LoDwingmax,"bx",label = "Max L/D: $(round(LoDwingmax,2)) at $(round(Vinf[LoDwingidxMax],2)) m/s
Re = $(round(Re[LoDwingidxMax],-3))
CL = $(round(CL_a[LoDwingidxMax],2))")
PyPlot.plot(Vinf,LoDanalytical,"r-", label = "Analytical Wing Only L/D")
LoDanalmax,LoDanalidxMax
PyPlot.plot(Vinf[i],LoDanalytical[i],"rx", label = "Analytical Max L/D: $(round(LoDanalmax,2)) at $(round(Vinf[i],2)) m/s
Re = $(round(Re[i],-3))
CL = $(round(CL_a[i],2))")
PyPlot.xlabel("Vinf")
PyPlot.ylabel("L/D")
PyPlot.legend(loc = 4,ncol = 2)
