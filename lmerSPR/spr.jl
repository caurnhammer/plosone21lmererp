## SESSION PREPARATION
pwd()

using CSV, DataFrames, DataFramesMeta, MixedModels, PooledArrays
using StatsBase: mean, std, zscore

include("spr_fns.jl")

## DATA PREPARATION
# Load data and scale continuous covariates
dt = read_spr_data("data/CAPSPR.csv");

dt.srp = dt.srp .* -1
dt.rcnoun = dt.rcnoun .* -1

# exclude data
dt = exclude_trial(dt[((dt.Region .!== "critical -2") .& (dt.Duplicated .!== "multi")),:], 50, 2500, 50, 6000)

#### TODO: zscore rcnoun / rcverb. Update formula. Run models anew.

## MODEL PREPARATION
contr = Dict(:Subject => Grouping(), :Item => Grouping())

f1 = @formula(logRT ~ 1 + Cloze + (1 + Cloze | Item) + (1 + Cloze | Subject))
f2 = @formula(logRT ~ 1 + srp + (1 + srp | Item) + (1 + srp | Subject))
f3 = @formula(logRT ~ 1 + rcnoun + (1 + rcnoun | Item) + (1 + rcnoun | Subject))
f3a = @formula(logRT ~ 1 + logrcnoun + (1 + logrcnoun | Item) + (1 + logrcnoun | Subject))
f4 = @formula(logRT ~ 1 + rcverb + (1 + rcverb | Item) + (1 + rcverb | Subject))
f5 = @formula(logRT ~ 1 + rcnoun + rcverb + (1 + rcnoun + rcverb | Item) + (1 + rcnoun + rcverb | Subject))
f6 = @formula(logRT ~ 1 + rcnoun + srp + (1 + rcnoun + srp | Item) + (1 + rcnoun + srp | Subject))
f7 = @formula(logRT ~ 1 + srp + rcnoun + inter + (1 + srp + rcnoun + inter | Item) + (1 + srp + rcnoun + inter | Subject))
f8 = @formula(logRT ~ 1 + Cloze + rcnoun + (1 + Cloze + rcnoun | Item) + (1 + Cloze + rcnoun | Subject))

f = f2

lmerSPR = fit_models(dt, f, contr);

numpred = length(f.rhs) - 2
coefs = [Symbol(x) for x in string.(f.rhs[1:numpred])];

lmerSPR = combine_datasets(dt, lmerSPR, coefs[2:length(coefs)])
lmerSPR = compute_RTs(lmerSPR, coefs);
estimates = names(lmerSPR)[occursin.("est_", names(lmerSPR))]
lmerSPR = compute_residuals(lmerSPR, estimates); 
lmerSPR = compute_coefs(lmerSPR, coefs); 

CSV.write("data/lmerSPR_Asrp.csv", lmerSPR);
