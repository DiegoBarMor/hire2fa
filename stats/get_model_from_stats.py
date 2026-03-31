import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

NBINS = dict(
    dist  = 200,
    angle = 360,
    dihed = 720,
)
RANGES = dict(
    dist  = (0.5, 6.0),
    angle = (0, np.pi),
    dihed = (-2*np.pi, 2*np.pi),
)

THRESHOLD_IMSHOW = 0.1

STATS_NAMES = [
    "op1", "op2", "c5", "c3", "c2", "o2", "o4",
    "buC6", "buN1", "buC2", "buO2", "buN3", "buC4", "buO4", "buC5",
    "bcC6", "bcN1", "bcC2", "bcO2", "bcN3", "bcC4", "bcN4", "bcC5",
    "baN9", "baC4", "baN3", "baC2", "baN1", "baC6", "baN6", "baC5", "baN7", "baC8",
    "bgN9", "bgC4", "bgN3", "bgC2", "bgN2", "bgN1", "bgC6", "bgO6", "bgC5", "bgN7", "bgC8",
]

THRESHOLD_PEAK_HEIGHT = {
    "op1": 0.2, "op2": 0.2, "c5": 0.05, "c3": 0.1, "c2": 0.06, "o2": 0.04, "o4": 0.08,
    "buC6": 0.2, "buN1": 0.2, "buC2": 0.2, "buO2": 0.2, "buN3": 0.2, "buC4": 0.2, "buO4": 0.2, "buC5": 0.2,
    "bcC6": 0.2, "bcN1": 0.2, "bcC2": 0.2, "bcO2": 0.2, "bcN3": 0.2, "bcC4": 0.2, "bcN4": 0.2, "bcC5": 0.2,
    "baN9": 0.2, "baC4": 0.2, "baN3": 0.2, "baC2": 0.2, "baN1": 0.2, "baC6": 0.2, "baN6": 0.2, "baC5": 0.2, "baN7": 0.2, "baC8": 0.2,
    "bgN9": 0.2, "bgC4": 0.2, "bgN3": 0.2, "bgC2": 0.2, "bgN2": 0.2, "bgN1": 0.2, "bgC6": 0.2, "bgO6": 0.2, "bgC5": 0.2, "bgN7": 0.2, "bgC8": 0.2,
}
PEAK_WIDTH_SIGMAS = {
    "op1": 8, "op2": 8, "c5": 8, "c3": 8, "c2": 8, "o2": 10, "o4": 10,
    "buC6": 8, "buN1": 8, "buC2": 8, "buO2": 8, "buN3": 8, "buC4": 8, "buO4": 8, "buC5": 8,
    "bcC6": 8, "bcN1": 8, "bcC2": 8, "bcO2": 8, "bcN3": 8, "bcC4": 8, "bcN4": 8, "bcC5": 8,
    "baN9": 8, "baC4": 8, "baN3": 8, "baC2": 8, "baN1": 8, "baC6": 8, "baN6": 8, "baC5": 8, "baN7": 8, "baC8": 8,
    "bgN9": 8, "bgC4": 8, "bgN3": 8, "bgC2": 8, "bgN2": 8, "bgN1": 8, "bgC6": 8, "bgO6": 8, "bgC5": 8, "bgN7": 8, "bgC8": 8,
}

NAMES_DETAILED = ["c5", "c3", "c2", "o2", "o4"]
REF_DETAILED = "c5"

NAMES_FIX_PERIODICITY = [ # dihedral peaks that get cut at the -pi|pi boundary
    "buC2", "buC4", "buC5", "buN3", "buO2", "buO4",
    "bcC2", "bcC4", "bcC5", "bcN3", "bcN4", "bcO2",
    "baC5", "baC6", "baC8", "baN1", "baN6", "baN7",
    "bgC5", "bgC6", "bgC8", "bgN1", "bgN7", "bgO6",
]

BUFFER_OUT = []

################### HISTOGRAMS 1D NUMERICS ###################
# ------------------------------------------------------------------------------
def _get_histograms_1d(df: pd.DataFrame) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    return dict(
        dist  = np.histogram(df["dist"],  bins = NBINS["dist"] , range = RANGES["dist"] ),
        angle = np.histogram(df["angle"], bins = NBINS["angle"], range = RANGES["angle"]),
        dihed = np.histogram(df["dihed"], bins = NBINS["dihed"], range = RANGES["dihed"]),
    )


# ------------------------------------------------------------------------------
def _fix_periodicity(atomname: str, geometry: str, points: np.ndarray) -> None:
    if geometry != "dihed": return
    if atomname not in NAMES_FIX_PERIODICITY: return

    quarter = len(points) // 4 # number of points between (n)pi and (n+1)pi
    points[:quarter] = points[2*quarter:3*quarter]
    points[2*quarter:3*quarter] = 0


# ------------------------------------------------------------------------------
def _calc_gaussian(x, normalization, mu, sigma):
    return normalization * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


# ------------------------------------------------------------------------------
def _fit_gaussians(
    hist1d, resname: str, atomname: str, geometry: str, report_fitting: bool
) -> tuple[float, float, float]:
    points,xedges = hist1d
    _fix_periodicity(atomname, geometry, points)

    x = np.linspace(*RANGES[geometry], NBINS[geometry])
    pos_peak = xedges[np.argmax(points)]

    g = lambda params: _calc_gaussian(x, *params) - points
    initial_guess_params = [np.max(points), pos_peak, 0.1]
    result = least_squares(g, initial_guess_params)
    normalization, mu, sigma = result.x
    sigma = abs(sigma)

    if report_fitting:
        print(f"...>>> Fitted Gaussian ({resname}, {atomname}, {geometry}): normalization={normalization:.2f}, mu={mu:.2f}, sigma={sigma:.2f}")
        ### 0 is a placeholder for the variant index
        BUFFER_OUT.append((atomname, geometry, 0, normalization, mu, sigma))

    return normalization, mu, sigma



################### HISTOGRAMS 2D NUMERICS ###################
# ------------------------------------------------------------------------------
def _get_histograms_2d(df: pd.DataFrame) -> tuple[tuple[np.ndarray, np.ndarray, np.ndarray],...]:
    hist_dist_angle = np.histogram2d(df["dist"], df["angle"],
        bins = (NBINS["dist"], NBINS["angle"]), range = (RANGES["dist"], RANGES["angle"])
    )
    hist_dist_dihed = np.histogram2d(df["dist"], df["dihed"],
        bins = (NBINS["dist"], NBINS["dihed"]), range = (RANGES["dist"], RANGES["dihed"])
    )
    hist_angle_dihed = np.histogram2d(df["angle"], df["dihed"],
        bins = (NBINS["angle"], NBINS["dihed"]), range = (RANGES["angle"], RANGES["dihed"])
    )
    return hist_dist_angle, hist_dist_dihed, hist_angle_dihed


# ------------------------------------------------------------------------------
def _normalize_2dhist(
    hist: tuple[np.ndarray, np.ndarray, np.ndarray], atomname: str
) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    points, xedges, yedges = hist
    points = np.log(points + 1)
    points /= np.max(points)
    points[points < THRESHOLD_IMSHOW] = 0

    relevant_cols = np.any(points > 0, axis = 0)
    relevant_rows = np.any(points > 0, axis = 1)

    relevant_cols_boundaries = np.logical_xor(relevant_cols, np.roll(relevant_cols, 1))
    relevant_rows_boundaries = np.logical_xor(relevant_rows, np.roll(relevant_rows, 1))

    idx_cols_start = 0 # default full range
    idx_cols_end   = len(relevant_cols_boundaries)
    idx_rows_start = 0
    idx_rows_end   = len(relevant_rows_boundaries)

    atomnames_exclude_trimming = [] # place here any atoms that behave weird
    if atomname not in atomnames_exclude_trimming:
        idx_cols_start += np.argmax(relevant_cols_boundaries)
        idx_cols_end   -= np.argmax(relevant_cols_boundaries[::-1])
        idx_rows_start += np.argmax(relevant_rows_boundaries)
        idx_rows_end   -= np.argmax(relevant_rows_boundaries[::-1])

    trimmed_points = points[idx_rows_start:idx_rows_end, idx_cols_start:idx_cols_end]
    extent = (
        xedges[idx_rows_start], xedges[idx_rows_end],
        yedges[idx_cols_start], yedges[idx_cols_end],
    )
    return trimmed_points.T, extent



################### HISTOGRAM 1D PLOTS ###################
# ------------------------------------------------------------------------------
def _plotsave_histogram1d(
    hist1d, resname: str, atomname: str, geometry: str,
    report_fitting: bool
) -> None:
    points,xedges = hist1d
    points_max = points.max()
    points = points.astype(float) / points_max
    _fix_periodicity(atomname, geometry, points)

    plt.figure()
    plt.plot(xedges[:-1], points, label = "data")

    max_iters = 2
    for _ in range(max_iters):
        if np.all(hist1d[0] < THRESHOLD_PEAK_HEIGHT[atomname] * points_max): break

        norm, mu, sigma = _fit_gaussians(
            hist1d, resname, atomname, geometry, report_fitting
        )
        x = np.linspace(*RANGES[geometry], NBINS[geometry])
        g = _calc_gaussian(x, norm, mu, sigma)
        g *= norm / (points_max * g.max())

        peak_start = mu - sigma*PEAK_WIDTH_SIGMAS[atomname]/2
        peak_end   = mu + sigma*PEAK_WIDTH_SIGMAS[atomname]/2
        mask = (peak_start <= x) & (x <= peak_end)
        hist1d[0][mask] = 0
        plt.plot(x, g, linestyle = "--", linewidth = 1.5, label = "gaussian fit")

    plt.title(f"{resname}: {atomname} ({geometry})")
    plt.legend()

    folder_name = "merged" if resname == "all" else "individual"
    plt.savefig(FOLDER_PLOTS / f"{folder_name}_{geometry}" / f"{atomname}.{resname}.png", dpi = 300)
    plt.close()



################### HISTOGRAM 2D PLOTS ###################
# ------------------------------------------------------------------------------
def _plotsave_histogram2d_hmap(hist2d, atomname: str, geo_0: str, geo_1: str) -> None:
    points , extent = _normalize_2dhist(hist2d, atomname)
    plt.figure()
    plt.imshow(points, extent = extent, aspect = "auto", origin = "lower")
    plt.xlabel(geo_0)
    plt.ylabel(geo_1)
    plt.title(f"{atomname}: {geo_0} vs {geo_1}")
    plt.savefig(FOLDER_PLOTS / f"hmaps_{geo_0}_{geo_1}" / f"{atomname}.png", dpi = 300)
    plt.close()



################### MAIN FUNCTIONS ###################
# ------------------------------------------------------------------------------
def plotsave_geometries_correlated(df: pd.DataFrame) -> None:
    print(">>> Plotting histograms...")
    for atomname in STATS_NAMES:
        df_atomname = df[df["atomname"] == atomname]
        if len(df_atomname) == 0: continue

        ##### 1D HISTOGRAMS #####
        histograms = _get_histograms_1d(df_atomname)
        for geometry in ["dist", "angle", "dihed"]:
            _plotsave_histogram1d(
                histograms[geometry], "all", atomname, geometry,
                report_fitting = True
            )

        ##### 2D HISTOGRAMS #####
        hist_dist_angle, hist_dist_dihed, hist_angle_dihed = _get_histograms_2d(df_atomname)

        _plotsave_histogram2d_hmap(hist_dist_angle,  atomname,  "dist", "angle")
        _plotsave_histogram2d_hmap(hist_dist_dihed,  atomname,  "dist", "dihed")
        _plotsave_histogram2d_hmap(hist_angle_dihed, atomname, "angle", "dihed")


# ------------------------------------------------------------------------------
def plotshow_detailed(df: pd.DataFrame) -> None:
    df_ref = df[df["atomname"] == REF_DETAILED]

    for atomname in NAMES_DETAILED:
        if atomname == REF_DETAILED: continue

        arr_ref = df_ref["dihed"].values
        df_atomname = df[df["atomname"] == atomname]
        arr_atomname = df_atomname["dihed"].values

        if len(arr_atomname) > len(arr_ref):
            arr_atomname = arr_atomname[:len(arr_ref)]
        elif len(arr_atomname) < len(arr_ref):
            arr_ref = arr_ref[:len(arr_atomname)]

        hist = np.histogram2d(
            arr_ref, arr_atomname,
            bins = (NBINS["dihed"], NBINS["dihed"]),
            range = (RANGES["dihed"], RANGES["dihed"]),
        )
        points , extent = _normalize_2dhist(hist, atomname)
        plt.figure()
        plt.imshow(points, extent = extent, aspect = "auto", origin = "lower")
        plt.xlabel(REF_DETAILED)
        plt.ylabel(atomname)
        plt.title(f"Dihed: {REF_DETAILED} vs {atomname}")
    plt.show()


# ------------------------------------------------------------------------------
def export_fitted_gaussians_to_csv() -> None:
    """Must be called after `plotsave_geometries_correlated` (relies on `BUFFER_OUT`)"""
    header = "atomname,geometry,variant,normalization,mu,sigma"
    buffer = [list(buf) for buf in BUFFER_OUT]

    for i,buf_this in enumerate(buffer):
        buf_this[3] = f"{buf_this[3]:.2f}" # normalization
        buf_this[4] = f"{buf_this[4]:.2f}" # mu
        buf_this[5] = f"{buf_this[5]:.2f}" # sigma

        if not i: continue
        buf_prev = buffer[i-1]
        if buf_this[:2] != buf_prev[:2]: continue # atomname, geometry
        buf_this[2] = buf_prev[2] + 1 # variant

    PATH_CSV_MODEL.write_text(
        '\n'.join((header,
            *(','.join(str(x) for x in buf) for buf in buffer)
        ))
    )


# ------------------------------------------------------------------------------
def pack_model_to_json() -> dict:
    def get_geometry_weight_pairs(df: pd.DataFrame, geometry: str) -> list[tuple[float, float]]:
        return [
            (round(row["mu"], 3), row["normalization"])
            for _, row in df[df["geometry"] == geometry].iterrows()
        ]

    model = {}
    df = pd.read_csv(PATH_CSV_MODEL)
    for atomname in df["atomname"].unique():
        df_atomname = df[df["atomname"] == atomname]
        pairs_dist  = get_geometry_weight_pairs(df_atomname, "dist")
        pairs_angle = get_geometry_weight_pairs(df_atomname, "angle")
        pairs_dihed = get_geometry_weight_pairs(df_atomname, "dihed")

        dist_m,  _ = pairs_dist[0]
        angle_m, _ = pairs_angle[0]
        sum_dws = sum(dihed_w for _,dihed_w in pairs_dihed)
        model[atomname] = [
            ### hardcoded assumption all entries consider only one unique "dist" peak and one unique "angle" peak
            ### so the entry's probability is only controlled by the dihedral variants
            [dihed_w/sum_dws, dist_m, angle_m, dihed_mu]
            for dihed_mu,dihed_w in pairs_dihed
        ]

    with open(PATH_JSON_MODEL, 'w') as file:
        json.dump(model, file)


# ------------------------------------------------------------------------------
def main():
    FOLDER_PLOTS.mkdir(exist_ok = True)
    for geometry in ["dist", "angle", "dihed"]:
        (FOLDER_PLOTS / f"merged_{geometry}").mkdir(exist_ok = True)

    (FOLDER_PLOTS / f"hmaps_dist_angle").mkdir(exist_ok = True)
    (FOLDER_PLOTS / f"hmaps_dist_dihed").mkdir(exist_ok = True)
    (FOLDER_PLOTS / f"hmaps_angle_dihed").mkdir(exist_ok = True)

    df_stats = pd.read_csv(PATH_CSV_STATS)

    # plotshow_detailed(df_stats)

    plotsave_geometries_correlated(df_stats)
    export_fitted_gaussians_to_csv() # must be called after plotsave_geometries_correlated
    pack_model_to_json()


################################################################################
if __name__ == "__main__":
    PATH_CSV_STATS = Path(sys.argv[1])
    PATH_CSV_MODEL = Path(sys.argv[2])
    PATH_JSON_MODEL = PATH_CSV_MODEL.with_suffix(".json")
    FOLDER_PLOTS = PATH_CSV_STATS.parent / "plots"
    FOLDER_PLOTS.mkdir(parents = True, exist_ok = True)
    PATH_CSV_MODEL.parent.mkdir(parents = True, exist_ok = True)
    main()


################################################################################
# python3 stats/get_model_from_stats.py testdata/stats.csv hire2fa/_data/model.csv
