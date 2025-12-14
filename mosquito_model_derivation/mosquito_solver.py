import numpy as np
import matplotlib.pyplot as plt

try:
    import rasterio
except ImportError as e:
    raise SystemExit(
        "This script requires rasterio to load the GeoTIFF.\n"
        "Install it with:  pip install rasterio"
    ) from e


# ===============================================================
# Global settings
# ===============================================================

TIF_PATH = "Malawan_WorldCover_100m_Map.tif"  # Path to the GeoTIFF file

DX = 1.0
DT = 1.0
NT = 20000

# Random-walk hop probabilities
P_E = 0.25
P_W = 0.15
P_N = 0.10
P_S = 0.10

# Seed location in FULL (uncropped) grid indices
SEED_ROW = 790
SEED_COL = 590

# Cropping window (around the seed) to reduce domain size
USE_CROP = True
CROP_HALF_HEIGHT = 250   # pixels in row direction
CROP_HALF_WIDTH = 200    # pixels in column direction


# ===============================================================
# 1. Load TIFF and build mask (forest + mangroves only)
# ===============================================================

def load_mask_from_tiff(path, valid_values=(10, 95)):
    """
    Load the WorldCover GeoTIFF and build a binary mask for the
    forest + mangrove classes.
    """
    with rasterio.open(path) as src:
        data = src.read(1)

    if valid_values is None:
        mask = (data != 0)
    else:
        mask = np.isin(data, np.asarray(valid_values))

    return mask.astype(float)


# ===============================================================
# 2. Single initial bump inside the (possibly cropped) mask
# ===============================================================

def make_initial_bump(mask, cy, cx, radius=5):
    """
    Place a single Gaussian bump at a specified location (cy, cx)
    in the CURRENT mask coordinates (which may be cropped).
    """
    ny, nx = mask.shape
    if not (0 <= cy < ny and 0 <= cx < nx):
        raise ValueError("Seed indices fall outside the current mask domain.")

    yy, xx = np.indices(mask.shape)

    # Gaussian bump
    r2 = (xx - cx) ** 2 + (yy - cy) ** 2
    bump = np.exp(-r2 / (2.0 * radius ** 2))    # exp(-r^2 / (2 sigma^2))

    # Restrict to land and renormalize
    bump *= mask
    total = bump.sum()
    if total == 0.0:
        raise ValueError("Chosen seed location lies completely outside the mask.")
    bump /= total

    print(f"Seed center (row, col) in local (cropped) coords: ({cy}, {cx})")
    return bump


# ===============================================================
# 3. Random walk with absorbing domain boundary
# ===============================================================

def step_random_walk(u, mask, pE, pW, pN, pS):
    """
    One time step of the discrete 2D random walk on the forest mask.
    Absorbing at the rectangular boundary, then clipped by the mask.
    """
    r = 1.0 - (pE + pW + pN + pS)
    if r < 0:
        raise ValueError("Probabilities must sum to ≤ 1.")

    u_new = np.zeros_like(u)

    i = slice(1, -1)
    j = slice(1, -1)

    uC = u[i, j]
    uE = u[i, 2:]
    uW = u[i, :-2]
    uN = u[:-2, j]
    uS = u[2:, j]

    u_new[i, j] = (
        r * uC +
        pE * uW +
        pW * uE +
        pN * uS +
        pS * uN
    )

    u_new *= mask
    return u_new


# ===============================================================
# 4. PDE with masked Laplacian (absorbing, consistent with walk)
# ===============================================================

def masked_laplacian(m, mask, dx):
    """
    5-point Laplacian with the same characteristic-function masking
    used in the random walk.
    """
    mC = m
    mE = np.roll(m, -1, axis=1)
    mW = np.roll(m, 1, axis=1)
    mN = np.roll(m, 1, axis=0)
    mS = np.roll(m, -1, axis=0)

    mE *= mask
    mW *= mask
    mN *= mask
    mS *= mask

    lap = (mE + mW + mN + mS - 4.0 * mC) / (dx * dx)
    lap *= mask
    return lap


def step_pde(m, mask, vx, vy, D, dx, dt):
    """
    Explicit time step for:

        ∂_t m = -vx ∂_x m - vy ∂_y m + D (∂_xx m + ∂_yy m)

    Using upwind advection and masked Laplacian, with absorbing
    boundaries consistent with the random walk.
    """
    i = slice(1, -1)
    j = slice(1, -1)

    mC = m[i, j]
    mE = m[i, 2:]
    mW = m[i, :-2]
    mN = m[:-2, j]
    mS = m[2:, j]

    if vx >= 0:
        dmx = (mC - mW) / dx
    else:
        dmx = (mE - mC) / dx

    if vy >= 0:
        dmy = (mC - mN) / dx
    else:
        dmy = (mS - mC) / dx

    lap = masked_laplacian(m, mask, dx)[i, j]

    m_new = m.copy()
    m_new[i, j] = mC + dt * (-vx * dmx - vy * dmy + D * lap)
    m_new *= mask
    return m_new


# ===============================================================
# 5. Main driver
# ===============================================================

def main():
    # -----------------------------------------------------------
    # Load FULL mask
    # -----------------------------------------------------------
    mask_full = load_mask_from_tiff(TIF_PATH)
    nrows_full, ncols_full = mask_full.shape
    print("Full mask size:", nrows_full, "x", ncols_full)

    # -----------------------------------------------------------
    # Optionally crop around the seed BEFORE simulation
    # -----------------------------------------------------------
    if USE_CROP:
        r0 = max(0, SEED_ROW - CROP_HALF_HEIGHT)
        r1 = min(nrows_full, SEED_ROW + CROP_HALF_HEIGHT)
        c0 = max(0, SEED_COL - CROP_HALF_WIDTH)
        c1 = min(ncols_full, SEED_COL + CROP_HALF_WIDTH)

        mask = mask_full[r0:r1, c0:c1]
        seed_row_local = SEED_ROW - r0
        seed_col_local = SEED_COL - c0

        print(f"Using cropped domain: rows [{r0}, {r1}), cols [{c0}, {c1})")
        print("Cropped mask size:", mask.shape)
    else:
        mask = mask_full
        seed_row_local = SEED_ROW
        seed_col_local = SEED_COL

    dx = DX
    dt = DT
    Nt = NT

    pE = P_E
    pW = P_W
    pN = P_N
    pS = P_S

    r = 1.0 - (pE + pW + pN + pS)
    print("Stay-put probability =", r)

    vx = (pE - pW) * dx / dt
    vy = (pS - pN) * dx / dt
    Dx = (pE + pW) * dx ** 2 / (2.0 * dt)
    Dy = (pN + pS) * dx ** 2 / (2.0 * dt)
    D = 0.5 * (Dx + Dy)

    print("Derived vx, vy, D =", vx, vy, D)

    # Single-seed initial condition at the desired Palawan location
    u0 = make_initial_bump(mask, seed_row_local, seed_col_local, radius=5)

    u_walk = u0.copy()
    m_pde = u0.copy()

    # -----------------------------------------------------------
    # Time stepping on the CROPPED domain only
    # -----------------------------------------------------------
    for _ in range(Nt):
        u_walk = step_random_walk(u_walk, mask, pE, pW, pN, pS)
        m_pde = step_pde(m_pde, mask, vx, vy, D, dx, dt)

    diff = u_walk - m_pde
    abs_L2 = np.sqrt(np.sum(diff ** 2))
    rel_L2 = abs_L2 / np.sqrt(np.sum(m_pde ** 2))

    mass_walk = u_walk.sum()
    mass_pde = m_pde.sum()

    print("\nAfter", Nt, "steps:")
    print("Absolute L2 error:", abs_L2)
    print("Relative L2 error:", rel_L2)
    print("Mass walk:", mass_walk)
    print("Mass PDE :", mass_pde)

    # ==========================================================
    # 6. FIGURE: flip mask for outline and label seed
    #    (no extra cropping needed – domain is already cropped)
    # ==========================================================

    fig, ax = plt.subplots(1, 3, figsize=(18, 6))

    # Flip vertically so outline matches the forest_coverage_map.py figure
    mask_plot = np.flipud(mask)
    u_walk_plot = np.flipud(u_walk)
    m_pde_plot = np.flipud(m_pde)
    diff_plot = np.flipud(diff)

    # Seed location in the flipped (plot) coordinates
    nrows, ncols = mask_plot.shape
    seed_row_plot = nrows - 1 - seed_row_local   # row index flips vertically
    seed_col_plot = seed_col_local

    # Shared color scale for walk + PDE panels
    vmax_main = max(u_walk_plot.max(), m_pde_plot.max())
    vmin_main = 0.0

    # Panel 1: Random walk
    im0 = ax[0].imshow(
        u_walk_plot,
        origin="upper",
        vmin=vmin_main,
        vmax=vmax_main,
        cmap="viridis",
    )
    ax[0].contour(mask_plot, levels=[0.5], colors="white", linewidths=0.8)
    ax[0].plot(
        seed_col_plot,
        seed_row_plot,
        marker="o",
        markersize=6,
        markeredgecolor="white",
        markerfacecolor="none",
    )
    ax[0].text(
        seed_col_plot + 40,
        seed_row_plot + 20,
        "Mosquito\nseed",
        color="white",
        fontsize=14,
        fontweight="bold",
        ha="left",
        va="bottom",
    )
    ax[0].set_aspect("equal")
    ax[0].set_title("Random walk density")
    # ax[0].set_xticks([])
    # ax[0].set_yticks([])
    plt.colorbar(im0, ax=ax[0])

    # Panel 2: PDE solution
    im1 = ax[1].imshow(
        m_pde_plot,
        origin="upper",
        vmin=vmin_main,
        vmax=vmax_main,
        cmap="viridis",
    )
    ax[1].contour(mask_plot, levels=[0.5], colors="white", linewidths=0.8)
    ax[1].plot(
        seed_col_plot,
        seed_row_plot,
        marker="o",
        markersize=6,
        markeredgecolor="white",
        markerfacecolor="none",
    )
    ax[1].text(
        seed_col_plot + 40,
        seed_row_plot + 20,
        "Mosquito\nseed",
        color="white",
        fontsize=14,
        fontweight="bold",
        ha="left",
        va="bottom",
    )
    ax[1].set_aspect("equal")
    ax[1].set_title("PDE density")
    # ax[1].set_xticks([])
    # ax[1].set_yticks([])
    plt.colorbar(im1, ax=ax[1])

    # Panel 3: Difference (walk - PDE)
    im2 = ax[2].imshow(
        diff_plot,
        origin="upper",
        cmap="viridis",
    )
    ax[2].contour(mask_plot, levels=[0.5], colors="white", linewidths=0.8)
    ax[2].plot(
        seed_col_plot,
        seed_row_plot,
        marker="o",
        markersize=6,
        markeredgecolor="white",
        markerfacecolor="none",
    )
    ax[2].text(
        seed_col_plot + 40,
        seed_row_plot + 30,
        "Mosquito\nseed",
        color="white",
        fontsize=14,
        fontweight="bold",
        ha="left",
        va="bottom",
    )
    ax[2].set_aspect("equal")
    ax[2].set_title("Difference")
    # ax[2].set_xticks([])
    # ax[2].set_yticks([])
    plt.colorbar(im2, ax=ax[2])

    fig.tight_layout()

    out_png = "malawan_mosquito_walk_vs_pde_1x3.png"
    fig.savefig(out_png, dpi=300)
    print("Saved figure to", out_png)


if __name__ == "__main__":
    main()
