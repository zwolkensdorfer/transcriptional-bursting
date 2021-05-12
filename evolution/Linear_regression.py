import pandas as pd
import numpy as np
from statsmodels import api as sm
from sklearn import metrics

def Fit_robust_linear_trend(df, x_var, y_var):
    x = sm.add_constant(df[x_var])
    y = df[y_var]
    model = sm.RLM(y, x, M=sm.robust.norms.HuberT(t=1.345))
    rlm_results = model.fit()
    results = {
        "intercept": rlm_results.params[0],
        "slope": rlm_results.params[1],
        "intercept_error": rlm_results.bse[0],
        "slope_error": rlm_results.bse[1],
        "intercept_pval": rlm_results.pvalues[0],
        "slope_pval": rlm_results.pvalues[1],
        "r2_unweighted": metrics.r2_score(y, rlm_results.fittedvalues),
        "r2_weighted": metrics.r2_score(y, rlm_results.fittedvalues, sample_weight=rlm_results.weights),
    }
    return pd.Series(results)
