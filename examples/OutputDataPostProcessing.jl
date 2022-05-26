#Postprocessing

module OutputDataPostProcessing

using CSV
using DataFrames
using Plots
using GLM
using Polynomials
using Unitful


function main()
    biasValues        = zeros(0)*u"V"
    staticCapacitance = zeros(0)*u"F"
    biasValues        = CSV.read("biasValues.csv"       , DataFrame; header=false)
    #I assume that the values in "biasValues" are equally spaced and that vector start with 0 (biasValues.Column1[0]=0):
    bias = biasValues.Column1.+biasValues.Column1[2]/2
    bias = bias[1:95]
    staticCapacitance  = CSV.read("staticCapacitance.csv", DataFrame; header=false)
    display(plot(bias,abs.(staticCapacitance.Column1[1:95])))


    #One can get from capacitance-voltage (CV) relationship the values of "build in potential" and "doping density"
    #For details see eq. 3.10 in https://in.ncu.edu.tw/ncume_ee/SchottkyDiode.htm  
    #We will try to calculate "build in potential" and "doping density" from linear fit after plotting (1/C^2) vs. V
    y    = abs.(staticCapacitance.Column1[1:95]).^-2
    bias = abs.(bias)

    #display(plot(bias,y))

    data = DataFrame(X=bias,Y=y)
    ols = lm(@formula(Y ~ X), data)

    coef(ols)
    model(x) = coef(ols)[1] + coef(ols)[2] * x
    build_in_potential = coef(ols)[1]/coef(ols)[2]
    print(build_in_potential)
    println("\n")
    scatter(data.X, data.Y)
    #display(plot(data.X, model.(data.X), legend=false))

    eps0       = 8.85*10^-14 *u"F/cm" 
    #relative permittivity for CIGS
    epsRelativ = 13.6
    #with and depth of device: be carefull! Check if this parameters are the same as in the script you use to get "biasValues" and "staticCapacitance"! 
    w_device = 1.0 *u"cm"  # width of device
    z_device = 1.0 *u"cm"  # depth of device
    S        = w_device*z_device #1cm2
    q        = 1.602176565*10^-19*u"C"
    #Na = 10^15cm^-3
    N_derived = 2/(coef(ols)[2]*q*epsRelativ*eps0*S^2)
    print(N_derived)
    println("\n")
end

end