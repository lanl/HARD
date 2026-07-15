from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import imageio.v2 as imageio
from tqdm import tqdm

# ----------------------------
# Input / output paths
# ----------------------------
input_dir = Path(".").expanduser()
output_dir = Path(".").expanduser()
output_dir.mkdir(parents=True, exist_ok=True)

# Create temporary directory for frames
frames_dir = output_dir / "gif_frames"
frames_dir.mkdir(exist_ok=True)

# ----------------------------
# Configuration
# ----------------------------
start_step = 0
end_step = 3001
step_increment = 100

# Generate list of timesteps
timesteps = list(range(start_step, end_step + 1, step_increment))

print(f"Creating GIF for timesteps: {start_step} to {end_step} (every {step_increment})")
print(f"Total frames: {len(timesteps)}")

# ----------------------------
# Analytical solution functions
# ----------------------------

def wave_function(p_star, rK, pK, cK, ggK):
    """Wave function f_K(p_star)."""
    if p_star <= pK:
        return (
            2.0*cK/(ggK - 1.0)
            * ((p_star/pK)**((ggK - 1.0)/(2.0*ggK)) - 1.0)
        )
    else:
        AK = 2.0/((ggK + 1.0)*rK)
        BK = (ggK - 1.0)/(ggK + 1.0)*pK
        return (p_star - pK)*np.sqrt(AK/(p_star + BK))


def pressure_equation(p_star, rL, pL, uL, cL, ggL, rR, pR, uR, cR, ggR):
    """Nonlinear equation for p_star."""
    fL = wave_function(p_star, rL, pL, cL, ggL)
    fR = wave_function(p_star, rR, pR, cR, ggR)
    return fL + fR + uR - uL


def compute_p_star(rL, pL, uL, cL, ggL, rR, pR, uR, cR, ggR, tol=1e-12):
    """Solve for the star-region pressure p_star."""
    p_low = 1e-12*max(pL, pR, 1.0)
    p_high = max(pL, pR, 1.0)

    def F(p):
        return pressure_equation(p, rL, pL, uL, cL, ggL, rR, pR, uR, cR, ggR)

    while F(p_high) < 0.0:
        p_high *= 2.0

    p_star = brentq(F, p_low, p_high, xtol=tol, rtol=tol)
    return p_star


def SodShockAnalytic2Gamma(rL, pL, uL, rR, uR, pR, xs, x0, T, ggL, ggR):
    """Analytical solution for a two-gamma Sod-like Riemann problem."""
    Nx = len(xs)
    x_interface = xs[x0]
    v_analytic = np.zeros((3, Nx), dtype="float64")

    if T == 0.0:
        mask_L = xs < x_interface
        mask_R = xs >= x_interface
        v_analytic[0, mask_L] = rL
        v_analytic[1, mask_L] = pL
        v_analytic[2, mask_L] = uL
        v_analytic[0, mask_R] = rR
        v_analytic[1, mask_R] = pR
        v_analytic[2, mask_R] = uR
        return v_analytic

    cL = np.sqrt(ggL*pL/rL)
    cR = np.sqrt(ggR*pR/rR)

    p_star = compute_p_star(rL, pL, uL, cL, ggL, rR, pR, uR, cR, ggR)

    fL = wave_function(p_star, rL, pL, cL, ggL)
    fR = wave_function(p_star, rR, pR, cR, ggR)

    u_star = 0.5*(uL + uR) + 0.5*(fR - fL)

    if not (p_star < pL and p_star > pR):
        raise ValueError(
            "This implementation assumes left rarefaction and right shock. "
            f"Got p_star={p_star}, pL={pL}, pR={pR}."
        )

    # Region R: right undisturbed state
    PR = p_star/pR
    S_shock = uR + cR*np.sqrt((ggR - 1.0 + PR*(ggR + 1.0))/(2.0*ggR))
    x_shock = x_interface + S_shock*T
    mask_R = xs >= x_shock
    v_analytic[0, mask_R] = rR
    v_analytic[1, mask_R] = pR
    v_analytic[2, mask_R] = uR

    # Region 2: right star state
    S_contact = u_star
    x_contact = x_interface + S_contact*T
    alphaR = (ggR + 1.0)/(ggR - 1.0)
    rho_star_R = rR*(1.0 + alphaR*PR)/(alphaR + PR)
    mask_2 = (xs >= x_contact) & (xs < x_shock)
    v_analytic[0, mask_2] = rho_star_R
    v_analytic[1, mask_2] = p_star
    v_analytic[2, mask_2] = u_star

    # Region 3: left star state
    rho_star_L = rL*(p_star/pL)**(1.0/ggL)
    c_star_L = np.sqrt(ggL*p_star/rho_star_L)
    S_fanright = u_star - c_star_L
    x_fanright = x_interface + S_fanright*T
    mask_3 = (xs >= x_fanright) & (xs < x_contact)
    v_analytic[0, mask_3] = rho_star_L
    v_analytic[1, mask_3] = p_star
    v_analytic[2, mask_3] = u_star

    # Region 4: inside left rarefaction fan
    S_fanleft = uL - cL
    x_fanleft = x_interface + S_fanleft*T
    mask_4 = (xs >= x_fanleft) & (xs < x_fanright)
    xi = (xs[mask_4] - x_interface)/T
    u4 = 2.0/(ggL + 1.0)*(cL + 0.5*(ggL - 1.0)*uL + xi)
    c4 = cL + 0.5*(ggL - 1.0)*(uL - u4)
    v_analytic[0, mask_4] = rL*(c4/cL)**(2.0/(ggL - 1.0))
    v_analytic[1, mask_4] = pL*(c4/cL)**(2.0*ggL/(ggL - 1.0))
    v_analytic[2, mask_4] = u4

    # Region L: left undisturbed state
    mask_L = xs < x_fanleft
    v_analytic[0, mask_L] = rL
    v_analytic[1, mask_L] = pL
    v_analytic[2, mask_L] = uL

    return v_analytic


# ----------------------------
# Initial conditions for analytical solution
# ----------------------------
gg = 1.4
rL, uL, pL = 1.0, 0.0, 1.0
rR, uR, pR = 0.125, 0.0, 0.1

Nx = 1000
X = 1.0
xs = np.linspace(0, X, Nx)
x0 = Nx // 2

# ----------------------------
# Plot settings
# ----------------------------
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 18
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["axes.labelsize"] = 18
plt.rcParams["axes.linewidth"] = 1.0
plt.rcParams["xtick.labelsize"] = 16
plt.rcParams["ytick.labelsize"] = 16

# Field configuration
# Column indices based on your data structure:
# 0:time, 1:cellid_x, 2:coord_x, 3:density, 4:pressure, 5:sound_speed,
# 6:specific_internal_energy, 7:total_energy_density, 8:psi, 9:material_id,
# 10:velocity_x, 11:momentum_density_x

fields = ["density", "pressure", "velocity_x"]
fields_j = [3, 4, 10]  # ✅ CORRECTED: velocity_x is at column 10, not 8

units = {
    "density": r" [$kg \ m^{-3}$]",
    "pressure": r" $\times 10^5$ [$Pa$]",
    "velocity_x": r" [$m/s$]",
}

ylims = {
    "density": (0.0, 1.5),
    "pressure": (0.0, 1.5),
    "velocity_x": (0.0, 1.2),
}


# ----------------------------
# Generate frames for each field
# ----------------------------
for field_idx, field in enumerate(fields):
    print(f"\nGenerating frames for {field}...")
    
    frame_files = []
    j = fields_j[field_idx]
    
    for timestep in tqdm(timesteps):
        # Format timestep with leading zeros
        timestep_str = f"{timestep:05d}"
        
        # Construct input filename
        input_file = input_dir / f"output-1D-0-{timestep_str}.csv"
        
        # Check if file exists
        if not input_file.exists():
            print(f"Warning: {input_file} not found, skipping...")
            continue
        
        # Read data
        try:
            df = pd.read_csv(input_file)
            header = df.columns.values
            mm_ray = np.array(df)
        except Exception as e:
            print(f"Error reading {input_file}: {e}")
            continue
        
        # Extract current time from the first column
        current_time = mm_ray[0, 0]
        
        # Extract spatial data - SIMPLIFIED, NO OFFSET
        i = 2  # coord_x column index
        x_mm = mm_ray[:, i]
        y_mm = mm_ray[:, j]
        x_label = header[i]
        y_label = header[j]
        
        # Compute analytical solution at current time
        v_analytic_current = SodShockAnalytic2Gamma(
            rL, pL, uL, rR, uR, pR, xs, x0, current_time, ggL=1.4, ggR=2.0
        )
        
        # Create plot
        fig, ax = plt.subplots(figsize=(8, 6))
        
        ax.plot(x_mm, y_mm, 'o', label="Sod (MM)", color='C0', markersize=4)
        ax.plot(xs, v_analytic_current[field_idx], linestyle="-", 
                linewidth=3, label="Analytical", color='C1')
        
        ax.set_xlim(0, 1)
        ax.set_ylim(*ylims[field])
        
        ax.set_xlabel(x_label + r" [$m$]")
        ax.set_ylabel(y_label + units[field])
        
        # Add time to title
        ax.set_title(f"Sod Test: {field} - Step {timestep} (t = {current_time:.4f} s)")
        ax.legend(loc='best')
        
        # Add time annotation in corner
        ax.text(0.02, 0.98, f"t = {current_time:.4f} s", 
                transform=ax.transAxes, fontsize=16,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        fig.tight_layout()
        
        # Save frame
        frame_file = frames_dir / f"{field}_frame_{timestep:05d}.png"
        fig.savefig(frame_file, dpi=150)
        plt.close(fig)
        
        frame_files.append(frame_file)
    
    # ----------------------------
    # Create GIF
    # ----------------------------
    if frame_files:
        print(f"\nCreating GIF for {field}...")
        
        images = []
        for frame_file in frame_files:
            images.append(imageio.imread(frame_file))
        
        # Save GIF
        gif_file = output_dir / f"SOD_MM_{field}_animation.gif"
        imageio.mimsave(gif_file, images, duration=0.5, loop=0)
        
        print(f"Saved: {gif_file}")
        
        # Optionally clean up frame files
        # Uncomment to delete frames after creating GIF:
        # for frame_file in frame_files:
        #     frame_file.unlink()
    else:
        print(f"No frames generated for {field}")

print("\n✅ All GIFs created successfully!")
print(f"Output directory: {output_dir}")
print(f"Frames saved in: {frames_dir}")
