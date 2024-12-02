"""
OBJECTIVE: VHF upper cockpit antenna mod in-flight response estimation

Steps:
1.- Characterization of the input, with VHF premod and postmod changes. 
    We have: 
            Premod:
            In-flight natural freq. w1 = 64Hz. 2nd mode at 246Hz and 3rd at 417Hz NOTE: Clean Air!
            Mass = 1.23 Kg
            Antenna tip response: measured in FT up to 10g's at highest speed or beta +-4º at lower speed.
            K = m*w1^2 = 1.23 * (2*π*64) ^2 = 198895 N/m 
            
            Postmod:
            In-flight natural freq. ? Hz.
            Mass = 2.538 Kg (2 times more than premod)  A delta of 1.308Kg is due to the whole dielectric cover.
            Antenna tip response: ? This is the Objective
            K = m*w2^2 = 2.538 * (2*π* ?) ^2 = ? N/m 

2.- Observations:
        1.- The in-flight VHF premod structural response gives a separated 1st mode on all flight conditions tested. 
            Therefore as hypothesis we can consider 1dof system as 1st bending 100% contribution to the response being fully separated freqs,
            the contribution of the 2nd and 3rd modes will be minor to the response.
        2.- VHF is located on the front of the A/C:
                No A/C geometry obstacles upwash, no wakes
                Boundary Layer effets are negligible
                Excitation will come from 4 sources:
                    1.-At beta zero and with beta non-zero (i.e.4º): Flow separation of the antenna itself --> unknown wideband spectrum but expected to be in the range of 0-250Hz. Tip vortex freq?
                        Whatever the spectrum is we know that the tip response acceleration is 10g so this will be the reference point to account unknown external loading with known resopnse.
                    2.-A/C mechanical vibration --> assumed to be negligible in comparison with 1
                    3.-A/C Boundary Layer interaction is expeted to be low. TBInvestigated.
                    4.-External aerodynamic turbulence: Here is the main contributor on "rough air" conditions. Needs to be estimated! 
                    But first to estimate the response in "clean air" conditions, zero excitation of this.
       
3.-Resolution:
    1.-Assume harmonic external excitation
    2.-Get the response of the premod as funtion of the generic amplitude of harmonic excitation signal
    3.-Determine the amplitude of the harmonic excitation that produce a response of 10gs
    4.-Assume that the amplitude of the harmonic excitation is the same for the postmod.
    5.-Get the level of response in g on the antenna tip....should be compared to the 10gs. 
        If we are in range there is no risk.
        If we are significant higher then there is some risk
        Compare against qualification test report.

METHODS:
1dof model based on:
Option 1:
by Harmonic oscilator considering damping solutions
Option 2:
by https://github.com/ONSAS/StructuralDynamicsODESolvers.jl/tree/master
Option 3:
to explore solving actual equations thru:
https://docs.sciml.ai/DiffEqDocs/dev/types/dynamical_types/#dynamical_prob
"""
#Option 1:
# Harmonic oscillator with damping is:
# ü + 2*ξ*ωn*ú + ωn^2*u = P0/m*sin(ωt)  w is the excitation freqency close to wn being ωn the natural frequency
# ξ = c/cr is the damping coefficient as the ratio of the dmaping c over the critical damping cr
# cr= 2*SQRT(K*m)=2*m*ωn
# r = ω/ωn

# Resolve 1dof, get u(t) = f (c,k,m,p(t)), get F0 that on the peak ü(t) produces 10g's measured in flight for pre-mod.
# u(t)postmod = f (c2,k2,m2,F0), get ü(t)max under same p(t).

#Premod key parameters:
f1 = 64.0       #Hz
ω1 = 2*pi*f1    #rad/sec
m1 = 1.23       #Kg
k1 = m1 * ω1^2  #N/m
c1 = 3          #Ns/m Assumed!
cr1= 2*m1*ω1    #Ns/m Assumed!
ξ1 = c1/cr1     #Assumed!
#Postmod key parameters:
# f2 =
# ω2 =
m2 = 2.538      #Kg
# k2 = m2 * ω2^2  #N/m
# c2 =
#External appllied aerodynamic force assumed to be equal in pre and post mods.
p0 = 0.0          #N 
ωn = ω1
ω = 0.95*ωn #frequcencia de excitation close to the natural freq.

#Initial conditions: null
# u(0) = u0
# ú(0) = ú0
#general sol:
r = ω/ωn
ξ = ξ1
k = k1
a = p0 / k * ((1-r^2)/((1-r^2)^2+(2*r*ξ)^2))
b = p0 / k * ((-2*r*ξ)/((1-r^2)^2+(2*r*ξ)^2))
#Null Initial conditions:
A = -a
B = -ωn/ω*(b*r + a*ξ) 
# Displacement:
u(t) = exp.(-ξ*ωn .*t) .*(A .* cos.(ω .*t)+ B .* sin.(ω .*t)) .+ a .* sin.(ω .* t) .+ b .* cos.(ω .*t)
# acceleration:
ü(t) = -ξ*ωn*exp.(-ξ*ωn .*t) .*(cos.(ω.*t).*(B*ω -ξ*ωn*A) .- sin.(ω .*t) .*(ξ*ωn*B + ω*A)) .+ exp.(-ξ*ωn .*t) .*(ω .*sin.(ω.*t).*(ξ*ωn*A-B*ω) .- ω.*cos.(ω.*t).*(ξ*ωn*B + A*ω))  .+ p0/k .* (ω^2/((1-r^2)^2 + (2*r*ξ)^2)) .* ((2*r*ξ .*cos.(ω .*t))-((1-r^2) .*sin.(ω .*t)))
#acceleration, only stationary part:
üs(t) = p0/k .* (ω^2/((1-r^2)^2 + (2*r*ξ)^2)) .* ((2*r*ξ .*cos.(ω .*t))-((1-r^2) .*sin.(ω .*t)))

t = 0:0.1:10

using GLMakie
fig = Figure(; size = (600, 400))
ax1 = Axis(fig[1, 1], yticklabelcolor =:black,     xlabel = "Time [sec]")
ax2 = Axis(fig[1, 1], yticklabelcolor =:dodgerblue, yaxisposition = :right)
hidespines!(ax2)
# hidexdecorations!(ax2)
lines!(ax1,t,u(t).*1000.0, color = :black, label = "Displacmeent in mm")
# lines!(t,üs(t)./9.81, color = :dodgerblue, label = "Stationary accel [g]")
lines!(ax2,t,ü(t)./9.81, color = :dodgerblue, label = "Accel [g]")
# legend(ax2)
colgap!(fig.layout, 5)
fig