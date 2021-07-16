##
# Christoph Aurnhammer <aurnhammer@coli.uni-saarland.de>
#
# Apply lmerERP method to ASBCD data
##

## SESSION PREPARATION
pwd()

# load packages
using CSV, DataFrames, DataFramesMeta, MixedModels, PooledArrays, Distributed
using StatsBase: mean, zscore

# load lmerERP functions
include("lmerERP.jl")

## READ DATA
dt = read_data("data/CAPExp.csv");

dt.srp = dt.srp .* -1;
dt.rcnoun = dt.rcnoun .* -1;

#dt[:arcctr] = Array(recode(dt[:Condition], ["A", "C"] => 0.5, ["B", "D"] => -0.5));
#dt[:clctr] = Array(recode(dt[:Condition], ["A", "B"] => 0.5, ["C", "D"] => -0.5));
#dt[:interctr] = dt[:arcctr] .* dt[:clctr];

# describe(dt)

## MODEL PREPARATION
contr = Dict(:Subject => Grouping(), :Item => Grouping())
# f0 = @formula(EEG ~ 1 + (1 | Item) + (1 | Subject))
f1 = @formula(EEG ~ 1 + cl + (1 + cl | Item) + (1 + cl | Subject))
f2 = @formula(EEG ~ 1 + srp + (1 + srp | Item) + (1 + srp | Subject))
f3 = @formula(EEG ~ 1 + rcnoun + (1 + rcnoun | Item) + (1 + rcnoun | Subject))
f4 = @formula(EEG ~ 1 + rcverb + (1 + rcverb | Item) + (1 + rcverb | Subject))
# f4b = @formula(EEG ~ 1 + rcverb + rcnoun + (1 + rcverb + rcnoun | Item) + (1 + rcverb + rcnoun | Subject))
# f5 = @formula(EEG ~ 1 + rcnoun + cl + (1 + rcnoun + cl | Item) + (1 + rcnoun + cl | Subject))
f6 = @formula(EEG ~ 1 + rcnoun + srp + (1 + rcnoun + srp | Item) + (1 + rcnoun + srp | Subject))
# f6 = @formula(EEG ~ 1 + rcnoun + rcverb + srp + (1 + rcnoun + rcverb + srp | Item) + (1 + rcnoun + rcverb + srp | Subject))

f = f2;

## MODEL FITTING
lmerERP = fit_models(dt, f, contr, parallel=true)

############ READ DATA
numpred = length(f.rhs) - 2
coefs = [Symbol(x) for x in string.(f.rhs[1:numpred])];

lmerERP = combine_datasets(dt, lmerERP, coefs[2:length(coefs)]); 
lmerERP = compute_waveforms(lmerERP, coefs);
estimates = names(lmerERP)[occursin.("est_", names(lmerERP))]
lmerERP = compute_residuals(lmerERP, estimates); 
lmerERP = compute_coefs(lmerERP, coefs); 

## Write to disk
CSV.write("data/lmerERP_A.csv", lmerERP)
