# import core packages
import warnings
warnings.filterwarnings("ignore")
from itertools import combinations

# import semi-core packages
import numpy as np
import pandas as pd

# import open2c libraries
import bioframe

import cooler
import cooltools
import cooltools.expected

from .genome import get_arms
# cvd == contacts-vs-distance
def get_cvd(clr, arms):
    cvd = cooltools.expected.diagsum_symm(
        clr=clr,
        view_df=arms,
        weight_name = None,
        ignore_diags=2
    )
    return cvd
def agg_cvd(cvd,clr):
    # Aggregate diagonals from different genomic regions together.
    # Since all three calcuated statistics are additive, they can be aggregated
    # easily via sum() functiond.
    cvd_agg = (
        cvd
        .groupby("diag")
        .agg(
            {"n_valid":"sum",
            "count.sum":"sum"
            })
        .reset_index()
    )
    # Convert indices of diagonals into genomic separation, expressed in basepairs.
    cvd_agg["s_bp"] = (
        cvd_agg["diag"]
        * clr.binsize
    )
    # Now we can calculate the average raw interaction counts and normalized contact frequencies.
    cvd_agg["count.avg"] = (
        cvd_agg["count.sum"]
        / cvd_agg["n_valid"]
    )
    #since pandas groupby replaces nan with zero
    cvd_agg.loc[0:1,:] = np.nan
    return cvd_agg
def logbin_smooth(cvd,clr):
    lb_cvd, lb_slopes, lb_distbins = cooltools.expected.logbin_expected(
        cvd,
        summary_name = "count.sum"
    )
    lb_cvd_agg, lb_slopes_agg = cooltools.expected.combine_binned_expected(
        lb_cvd,
        Pc_name='count.avg',
        binned_exp_slope=lb_slopes
    )
    lb_cvd_agg['s_bp'] = lb_cvd_agg['diag.avg'] * clr.binsize
    lb_slopes_agg['s_bp'] = lb_slopes_agg['diag.avg'] * clr.binsize
    return lb_cvd_agg
def get_ps_curve(filep, ref_name):
    clr = cooler.Cooler(filep)
    arms = get_arms(ref_name, clr)
    cvd = get_cvd(clr, arms)
    cvd_agg = agg_cvd(cvd,clr)
    return cvd_agg
def get_lb_ps_curve(filep, ref_name):
    clr = cooler.Cooler(filep)
    arms = get_arms(ref_name, clr)
    cvd = get_cvd(clr, arms)
    lb_cvd_agg = logbin_smooth(cvd, clr)
    return lb_cvd_agg