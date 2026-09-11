import numpy as np
import pandas as pd
import scipy
import statsmodels.formula.api as smf
from statsmodels.nonparametric.smoothers_lowess import lowess


def gc_correction(df):
    df["rdr"] = df["reads"]

    # Filtering and sorting
    df_regression = df[df["valid"]].copy()

    df_regression.sort_values(by="gc", inplace=True)

    df_regression = modal_quantile_regression(df_regression, lowess_frac=0.2)

    df["gc_correction"] = np.nan

    df["rdr_cor"] = np.nan

    df.loc[df_regression.index, "gc_correction"] = df_regression["modal_curve"]

    df["rdr_cor"] = df["rdr"] / df["gc_correction"]

    return df


def modal_quantile_regression(df_regression, lowess_frac=0.2, degree=2, knots=(0.38,)):
    """
    Fits a B-spline polynomial curve through the "modal" quantile of the data:
    * Runs quantile regression to fit a B-spline curve for each percentile 10-90
    * Estimates the modal quantile as the quantile where difference in AUC is minimized
    * Uses the curve fit to this modal quantile for normalization

    Parameters:
        df_regression: pandas.DataFrame with at least columns [chr, start, end, rdr, gc]
        lowess_frac: float, fraction of data used to estimate each y-value in Lowess smoothing of AUC curve
        degree: int, degree of polynomial to fit to each section of the B-spline curve
        knots: list of floats, GC values where B-spline polynomial is allowed to change

    Returns:
        pandas.DataFrame with additional columns
            modal_curve: modal curve's predicted # rdr for GC value in this row
            modal_quantile: quantile selected as the mode (should be the same for all bins)
            modal_corrected: corrected read count (i.e., rdr / modal_curve)
    """

    q_range = range(10, 91, 1)
    quantiles = np.array(q_range) / 100
    quantile_names = [str(x) for x in q_range]

    # need at least 3 values to compute the quantiles
    if len(df_regression) < 10 or sum(df_regression["rdr"]) < 100:
        df_regression["modal_quantile"] = None
        df_regression["modal_curve"] = None
        df_regression["modal_corrected"] = None
        return df_regression

    poly_quantile_model = smf.quantreg(
        f"rdr ~ bs(gc, degree={degree}, knots={knots}, include_intercept = True)",
        data=df_regression,
    )

    poly_quantile_fit = [
        poly_quantile_model.fit(q=q, max_iter=10000) for q in quantiles
    ]

    poly_quantile_predict = [
        poly_quantile_fit[i].predict(df_regression) for i in range(len(quantiles))
    ]

    poly_quantile_params = pd.DataFrame()

    for i in range(len(quantiles)):
        df_regression[quantile_names[i]] = poly_quantile_predict[i]
        poly_quantile_params[quantile_names[i]] = poly_quantile_fit[i].params

    # integration and mode selection
    gc_min = df_regression["gc"].quantile(q=0.10)
    gc_max = df_regression["gc"].quantile(q=0.90)

    true_min = df_regression["gc"].min()
    true_max = df_regression["gc"].max()

    poly_quantile_integration = np.zeros(len(quantiles) + 1)

    # form (k+1)-regular knot vector
    repeats = degree + 1

    my_t = np.r_[[true_min] * repeats, knots, [true_max] * repeats]

    for i in range(len(quantiles)):
        # compose params into piecewise polynomial
        params = poly_quantile_params[quantile_names[i]].to_numpy()
        pp = scipy.interpolate.PPoly.from_spline((my_t, params[1:] + params[0], degree))

        # compute integral
        poly_quantile_integration[i + 1] = pp.integrate(gc_min, gc_max)

    # find the modal quantile
    distances = poly_quantile_integration[1:] - poly_quantile_integration[:-1]

    df_dist = pd.DataFrame(
        {
            "quantiles": quantiles,
            "quantile_names": quantile_names,
            "distances": distances,
        }
    )
    dist_max = df_dist["distances"].quantile(q=0.95)
    df_dist_filter = df_dist[df_dist["distances"] < dist_max].copy()
    df_dist_filter["lowess"] = lowess(
        df_dist_filter["distances"],
        df_dist_filter["quantiles"],
        frac=lowess_frac,
        return_sorted=False,
    )

    modal_quantile = df_dist_filter.set_index("quantile_names")["lowess"].idxmin()

    # add values to table
    df_regression["modal_quantile"] = modal_quantile

    df_regression["modal_curve"] = df_regression[modal_quantile]

    return df_regression
