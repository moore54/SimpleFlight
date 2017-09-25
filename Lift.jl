#KRM

include("atm.jl")
function lift_coefficient(weight,Sref,Vinf,h=1.2)
    rho, P,T,a,mu = stdatm(h)
    CL = 2*weight/(rho*Vinf^2*Sref)#/sqrt.(1-M^2)
    q = 0.5*rho*Vinf^2
    return CL,q
end
