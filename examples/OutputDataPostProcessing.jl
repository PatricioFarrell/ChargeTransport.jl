#Postprocessing

using CSV
using DataFrames
using Plots
using GLM
using Polynomials

biasValues        = zeros(0)
staticCapacitance = zeros(0)
biasValues        = CSV.read("biasValues.csv"       , DataFrame; header=false)
#I assume that the values in "biasValues" are equally spaced and that vector start witk 0 (biasValues.Column1[0]=0):
bias = biasValues.Column1.+biasValues.Column1[2]/2
bias = bias[1:100]
staticCapacitance  = CSV.read("staticCapacitance.csv", DataFrame; header=false)
plot(bias,abs.(staticCapacitance.Column1))

y    = abs.(staticCapacitance.Column1).^-2
bias = abs.(bias)

plot(bias,y)

data = DataFrame(X=bias,Y=y)
ols = lm(@formula(Y ~ X), data)

coef(ols)
model(x) = coef(ols)[1] + coef(ols)[2] * x
#1.4775483873006742e7
build_in_potential = coef(ols)[1]/coef(ols)[2]
scatter(data.X, data.Y)
plot(data.X, model.(data.X), legend=false)

eps0       = 8.85*10^-14
epsRelativ = 13.6
S          = 0.5*0.5
q          = 1.602176565*10^-19
N_derived = 2/(coef(ols)[2]*q*epsRelativ*eps0*S^2)