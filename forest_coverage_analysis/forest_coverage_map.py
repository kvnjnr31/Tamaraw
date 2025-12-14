"""
forest_coverage_map.py

Estimate basic forest coverage statistics from a WorldCover raster and
compute a box-counting fractal dimension for the forest mask.

This code produces a one-by-two PNG:

  - LEFT panel:  Forest + mangroves map tile (~100 m resolution)
  - RIGHT panel: Dimension-versus-"embedding" curve and a randomized
                 forest mask for comparison. 
"""

import rasterio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from typing import Tuple

WORLD_PATH = "Malawan_WorldCover_100m_Map.tif"
MAX_DE = 200        # max embedding index to show on the x-axis
N_PLOT_POINTS = 25  # max number of points to show per curve in right panel

# ----------------------------------------------------------------------
# Box-counting utilities
# ----------------------------------------------------------------------
def boxcount(mask: np.ndarray, k: int) -> int:
    """
    Count the number of non-empty k×k boxes needed to cover the 1-valued
    pixels in 'mask', using a simple explicit for-loop implementation.
    """
    nrows = (mask.shape[0] // k) * k
    ncols = (mask.shape[1] // k) * k

    if nrows == 0 or ncols == 0:
        return 0

    Zc = mask[:nrows, :ncols]

    S = []
    for i in range(0, nrows, k):
        row_blocks = []
        for j in range(0, ncols, k):
            block = Zc[i:i+k, j:j+k]
            row_blocks.append(block.sum())
        S.append(row_blocks)

    S = np.array(S)
    return np.count_nonzero(S)


def compute_boxcount_curve(mask: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the box-counting curve (k, N(k)) using all integer box sizes.
    """
    rows, cols = mask.shape
    min_side = min(rows, cols)

    eps = np.arange(1, min_side + 1, dtype=int)
    N_eps = np.array([boxcount(mask, int(k)) for k in eps], dtype=float)

    valid = (N_eps > 0) & (N_eps < (rows * cols))
    return eps[valid], N_eps[valid]


def fit_fractal_dimension(eps: np.ndarray, N_eps: np.ndarray) -> Tuple[float, float]:
    """
    Linear least-squares fit to log N(k) versus log(1/k):

        log N(k) ≈ D · log(1/k) + c
    """
    if eps.size < 2:
        raise ValueError("Need at least two scales (k) to fit a line.")

    x = np.log(1.0 / eps.astype(float))
    y = np.log(N_eps.astype(float))

    A = np.vstack([x, np.ones_like(x)]).T
    m, b = np.linalg.lstsq(A, y, rcond=None)[0]
    return m, b


def incremental_dimension(eps: np.ndarray, N_eps: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Incremental estimate of the fractal dimension as a function of
    an "embedding index" k (number of scales used in the fit).
    """
    m = len(eps)
    if m < 2:
        raise ValueError("Need at least two scales (k) to compute incremental dimension.")

    x = np.log(1.0 / eps.astype(float))
    y = np.log(N_eps.astype(float))

    De_list = []
    D_list = []

    for k in range(2, m + 1):
        A = np.vstack([x[:k], np.ones(k)]).T
        slope_k, _ = np.linalg.lstsq(A, y[:k], rcond=None)[0]
        De_list.append(k)
        D_list.append(slope_k)

    return np.asarray(De_list, dtype=float), np.asarray(D_list, dtype=float)


def truncate_De(De: np.ndarray, D_seq: np.ndarray, max_de: int = MAX_DE) -> Tuple[np.ndarray, np.ndarray]:
    """Truncate (De, D_seq) to De <= max_de."""
    if De.size == 0:
        return De, D_seq
    mask = De <= max_de
    return De[mask], D_seq[mask]


def downsample_for_plot(De: np.ndarray, D_seq: np.ndarray, n_points: int = N_PLOT_POINTS) -> Tuple[np.ndarray, np.ndarray]:
    """
    Downsample (De, D_seq) to at most n_points for visualization,
    using approximately evenly spaced indices.
    """
    if De.size == 0 or De.size <= n_points:
        return De, D_seq

    idx = np.linspace(0, De.size - 1, n_points, dtype=int)
    idx = np.unique(idx)
    return De[idx], D_seq[idx]


# ----------------------------------------------------------------------
# Main script: read raster, build mask, stats, fractal dimension
# ----------------------------------------------------------------------

# ---- READ WORLD COVER MAP AT NATIVE RESOLUTION -----------------------
with rasterio.open(WORLD_PATH) as src:
    data = src.read(1)
    bounds = src.bounds
    orig_h = src.height
    orig_w = src.width

print("Raster shape (rows, cols):", data.shape)

# -------- FOREST MASK (10=tree, 95=mangrove) --------
forest_mask = np.isin(data, [10, 95]).astype(np.uint8)

rows, cols = forest_mask.shape
total_cells = forest_mask.size
forest_cells = int(forest_mask.sum())
forest_frac = float(forest_cells) / float(total_cells)

print("Total cells:", total_cells)
print("Forest cells:", forest_cells)
print("Forest fraction:", forest_frac, " (", forest_frac * 100.0, "% )")

# -------- Coverage Statistics (km^2) --------
dlat = (bounds.top - bounds.bottom) / float(orig_h)
dlon = (bounds.right - bounds.left) / float(orig_w)
lat_c = 0.5 * (bounds.top + bounds.bottom)

deg_to_km_lat = 111.32
deg_to_km_lon = 111.32 * np.cos(np.deg2rad(lat_c))

pixel_area_km2 = dlat * deg_to_km_lat * dlon * deg_to_km_lon

total_area_km2 = total_cells * pixel_area_km2
forest_area_km2 = forest_cells * pixel_area_km2

print("Approx total area (km^2):", total_area_km2)
print("Approx forest area (km^2):", forest_area_km2)

# -------- BOX-COUNTING FRACTAL DIMENSION: FOREST MASK --------
eps_forest, N_eps_forest = compute_boxcount_curve(forest_mask)

if eps_forest.size >= 2:
    D_forest_global, intercept_forest = fit_fractal_dimension(eps_forest, N_eps_forest)
    print("Global box-counting fractal dimension D (forest mask, LSQ fit):", D_forest_global)
    De_forest, D_seq_forest = incremental_dimension(eps_forest, N_eps_forest)
else:
    D_forest_global = np.nan
    De_forest = np.array([])
    D_seq_forest = np.array([])
    print("Not enough scales to estimate fractal dimension for forest mask.")

# -------- BOX-COUNTING FRACTAL DIMENSION: RANDOMIZED MASK (NULL MODEL) --------
rng = np.random.default_rng(seed=42)
random_mask = (rng.random((rows, cols)) < forest_frac).astype(np.uint8)

eps_rand, N_eps_rand = compute_boxcount_curve(random_mask)

if eps_rand.size >= 2:
    D_rand_global, intercept_rand = fit_fractal_dimension(eps_rand, N_eps_rand)
    print("Global box-counting fractal dimension D (random mask, LSQ fit):", D_rand_global)
    De_rand, D_seq_rand = incremental_dimension(eps_rand, N_eps_rand)
else:
    D_rand_global = np.nan
    De_rand = np.array([])
    D_seq_rand = np.array([])
    print("Not enough scales to estimate fractal dimension for random mask.")

# ---- Truncate sequences to De <= MAX_DE for asymptotics and plotting ---
De_forest_trunc, D_seq_forest_trunc = truncate_De(De_forest, D_seq_forest, MAX_DE)
De_rand_trunc,   D_seq_rand_trunc   = truncate_De(De_rand,   D_seq_rand,   MAX_DE)

# Asymptotic forest D: mean of last 20% of truncated points
if D_seq_forest_trunc.size > 0:
    n_tail = max(1, int(0.2 * D_seq_forest_trunc.size))
    D_forest_tail_mean = float(D_seq_forest_trunc[-n_tail:].mean())
    print("Asymptotic D (forest, mean of last 20% of D_k, De<=MAX_DE):", D_forest_tail_mean)
else:
    D_forest_tail_mean = np.nan

# ---- Downsample for visualization (max ~25 points per curve) ----------
De_forest_plot, D_seq_forest_plot = downsample_for_plot(De_forest_trunc, D_seq_forest_trunc, N_PLOT_POINTS)
De_rand_plot,   D_seq_rand_plot   = downsample_for_plot(De_rand_trunc,   D_seq_rand_trunc,   N_PLOT_POINTS)

# ----------------------------------------------------------------------
# one-by-two FIGURE
# ----------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

# ===== LEFT PANEL: FOREST MAP =====
ax_map = axes[0]

cmap = mcolors.ListedColormap(["#eaf5e1", "#004b23"])
bounds_norm = [0, 0.5, 1]
norm = mcolors.BoundaryNorm(bounds_norm, cmap.N)

ax_map.imshow(forest_mask, cmap=cmap, norm=norm, interpolation="nearest")
ax_map.set_title("Malawan tile: forest + mangroves (~100 m)")
ax_map.set_xlabel("Column index")
ax_map.set_ylabel("Row index")

forest_patch = mpatches.Patch(color="#004b23", label="Forest (1)")
nonforest_patch = mpatches.Patch(color="#eaf5e1", label="Non-forest (0)")
ax_map.legend(
    handles=[forest_patch, nonforest_patch],
    loc="lower left",
    frameon=True,
    facecolor="white",
    edgecolor="black",
)

# ===== RIGHT PANEL: DIMENSION VS "EMBEDDING DIMENSION" =====
ax_corr = axes[1]

if De_rand_plot.size > 0:
    ax_corr.plot(
        De_rand_plot,
        D_seq_rand_plot,
        "o-",
        color="black",
        label="Random mask (white-noise example)",
    )

if De_forest_plot.size > 0:
    ax_corr.plot(
        De_forest_plot,
        D_seq_forest_plot,
        "o:",
        color="black",
        label="Forest mask (incremental fit)",
    )

# Bold dashed asymptotic lines
ax_corr.axhline(
    2.0,
    linestyle="--",
    color="black",
    linewidth=2.0,
    label="White-noise limit D = 2",
)

if not np.isnan(D_forest_tail_mean):
    ax_corr.axhline(
        D_forest_tail_mean,
        linestyle="--",
        color="black",
        linewidth=2.0,
        alpha=0.8,
        label=f"Forest asymptotic D ≈ {D_forest_tail_mean:.2f}",
    )
    ax_corr.text(
        MAX_DE * 0.98,
        D_forest_tail_mean,
        f"{D_forest_tail_mean:.2f}",
        fontsize=10,
        ha="right",
        va="center",
    )

ax_corr.set_xlabel("Embedding Dimension: De (fit index)")
ax_corr.set_ylabel("Dimension D (incremental slope)")
ax_corr.set_title("Dimension vs Embedding Dimension (box-counting)")

if (De_forest_trunc.size > 0) or (De_rand_trunc.size > 0):
    xmin = 1.5
    xmax = MAX_DE + 5
    ymin = 0.5
    ymax = max(
        D_seq_forest_trunc.max() if De_forest_trunc.size > 0 else 0,
        D_seq_rand_trunc.max() if De_rand_trunc.size > 0 else 0,
        2.0,
        D_forest_tail_mean if not np.isnan(D_forest_tail_mean) else 0,
    ) + 0.3
    ax_corr.set_xlim(xmin, xmax)
    ax_corr.set_ylim(ymin, ymax)

# Legend moved to lower right to reduce clutter
ax_corr.legend(loc="lower right", frameon=False)

fig.tight_layout()

out_png = "malawan_forest_corrdim_1x2.png"
fig.savefig(out_png, dpi=200)
print("Saved 1×2 forest + dimension-vs-embedding figure to", out_png)
