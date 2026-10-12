"""Dev-only cross-check of Little's test on airquality (needs numpy/scipy).

Maximises the observed-data normal log-likelihood directly with scipy, then
computes Little's d2 from those ML estimates and compares with our EM and with
R naniar (whose EM, norm::em.norm, stops at a looser convergence criterion).
Run from librestats/: ../.venv-dev/bin/python tools/check_little_ml.py
"""
import sys, os
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, "tests"); sys.path.insert(0, "src/pythonpath")
from test_missing import airquality
from librestats import missing as mi

rows = airquality()
X = np.array([[np.nan if v is None else v for v in r] for r in rows])
n, p = X.shape
mu0, sd0 = np.nanmean(X, 0), np.nanstd(X, 0)
Z = (X - mu0) / sd0                     # optimise on standardized scale
pats = {}
for i, r in enumerate(~np.isnan(Z)):
    pats.setdefault(tuple(r), []).append(i)
tri = np.tril_indices(p)

def unpack(theta):
    mu = theta[:p]; L = np.zeros((p, p)); L[tri] = theta[p:]
    L[np.diag_indices(p)] = np.exp(L[np.diag_indices(p)])
    return mu, L @ L.T

def nll(theta):
    mu, S = unpack(theta); tot = 0.0
    for obs, idx in pats.items():
        o = np.array(obs); Y = Z[np.ix_(idx, o)] - mu[o]; So = S[np.ix_(o, o)]
        sign, logdet = np.linalg.slogdet(So)
        tot += 0.5 * (len(idx) * logdet + np.einsum("ij,ij->", Y @ np.linalg.inv(So), Y))
    return tot

theta0 = np.concatenate([np.zeros(p), np.eye(p)[tri]]); theta0[p:][tri[0] == tri[1]] = 0.0
best = minimize(nll, theta0, method="BFGS", options={"gtol": 1e-10, "maxiter": 20000})
mu, S = unpack(best.x)
mu = mu * sd0 + mu0; S = S * np.outer(sd0, sd0)

d2 = 0.0
for obs, idx in pats.items():
    o = np.array(obs); diff = X[np.ix_(idx, o)].mean(0) - mu[o]
    d2 += len(idx) * diff @ np.linalg.solve(S[np.ix_(o, o)], diff)
em = mi.em_mvnorm(rows)
print("scipy converged:", best.success, "| max |mean diff| vs our EM:", np.max(np.abs(mu - em.mean)))
print(f"Little d2 from scipy ML : {d2:.10f}")
print(f"Little d2 from our EM   : {mi.little_mcar(rows).statistic:.10f}")
print(f"naniar (norm::em.norm)  : 35.1061288690")
