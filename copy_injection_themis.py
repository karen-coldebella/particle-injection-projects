# ======================================================================
# PURPOSE:
# Read and process THEMIS FGM data (CDF) to identify particle injections
# using magnetic field signatures (dipolarization).
# ======================================================================

# %% ---------------------------- IMPORTS ----------------------------
import datetime
import os
import cdflib
from cdflib import cdfepoch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata


# %% ---------------------------- 1. TIME INTERVAL -------------------
start_time = datetime.datetime(2007, 3, 23, 0, 0)
end_time   = datetime.datetime(2007, 3, 24, 0, 0)

start_plot_dt = datetime.datetime(2007, 3, 23, 11, 0)
end_plot_dt   = datetime.datetime(2007, 3, 23, 12, 0)


# %% ---------------------------- 2. DIRECTORY -----------------------
# Base directory (project root)
BASE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..')
)

# FGM instrument
directory_FGM_thd = os.path.join(BASE_DIR, 'Data', 'FGM')
file_FGM_thd = 'thd_l2_fgm_20070323_v01.cdf'
path_FGM_thd = os.path.join(directory_FGM_thd, file_FGM_thd)

# SST instrument
directory_SST_thd = os.path.join(BASE_DIR, 'Data', 'SST')
file_SST_thd = 'thd_l2_sst_20070323_v01.cdf'
path_SST_thd = os.path.join(directory_SST_thd, file_SST_thd)


# %% ---------------------------- 3. OPEN CDF ------------------------
#FGM instrument
cdf_FGM_thd = cdflib.CDF(path_FGM_thd)
print("Variables in CDF:", cdf_FGM_thd.cdf_info().zVariables)

#SST instrument
cdf_SST_thd = cdflib.CDF(path_SST_thd)
print("Variables in CDF:", cdf_SST_thd.cdf_info().zVariables)


# %% ---------------------------- 4. READ VARIABLES FGM ------------------
#FGM instrument
epoch_raw_FGM = cdf_FGM_thd.varget('thd_fgs_time')  # em segundos desde 1970
Btotal    = cdf_FGM_thd.varget('thd_fgs_btotal')
Bgse      = cdf_FGM_thd.varget('thd_fgs_gse')

# sanity check
#print("Epoch sample:", epoch_raw_FGM[:5])
#print("Btotal sample:", Btotal[:5])
#print("Bgse sample:\n", Bgse[:5, :])

# Covert epoch to datetime
# Unix epoch (seconds since 1970-01-01)
epoch_dt_FGM = pd.to_datetime(epoch_raw_FGM, unit='s', origin='unix')

# create DataFrame
df_thed_FGM = pd.DataFrame(
    {
        'Btotal': Btotal,
        'Bx': Bgse[:, 0],
        'By': Bgse[:, 1],
        'Bz': Bgse[:, 2],
    },
    index=pd.DatetimeIndex(epoch_dt_FGM)
)
df_thed_FGM.index.name = 'time'

# check
print(df_thed_FGM.index.min())
print(df_thed_FGM.index.max())




# %% ---------------------------- 5. READ VARIABLES SST ------------------
# SST instrument (Plasma Science Ion Flux – PSIF)

# Central time of each SST data sample
# UTC, in seconds since 1970-01-01 00:00:00 (Unix time)
# Defines the main temporal axis of the instrument
epoch0 = cdf_SST_thd.varget('thd_psif_time')

# Duration of each data sample (integration time)
# Unit: seconds
# Typically ~3 s in Full Mode
# Used to estimate the temporal window of the measurement
delta_time = cdf_SST_thd.varget('thd_psif_delta_time')

# Ion differential energy flux (energy spectrogram)
# Dimension: (time, energy) → (N, 16)
# Unit: eV / (cm^2 · s · sr · eV)
# Logarithmic scale
# Main variable for identifying energetic particle injections
ion_flux = cdf_SST_thd.varget('thd_psif_en_eflux')

# Central energy of each spectrogram channel
# Unit: eV
# Dimension: (time, energy), but the 16 energy channels are fixed in time
# Used as the Y-axis of the spectrogram
energy = cdf_SST_thd.varget('thd_psif_en_eflux_yaxis')

# Integrated (non-spectral) ion particle flux
# Flux vector with three components
# Unit: # / (cm^2 · s)
# Linear scale
# Suitable for time series analysis and comparison with FGM measurements
vector_ion_flux = cdf_SST_thd.varget('thd_psif_flux')

# Labels of the ion flux vector components
# Define the physical meaning of each column of thd_psif_flux
# e.g., instrumental directions or DSL components
vector_flux_labls = cdf_SST_thd.varget('thd_psif_flux_labl')


# -------------------- INSPECTION --------------------
var_info = cdf_SST_thd.varattsget('thd_psif_time')
print("Informações da variável 'thd_psif_time':", var_info)


variables = {
    'epoch0': epoch0,
    'delta_time': delta_time,
    'ion_flux': ion_flux,
    'energy': energy,
    'vector_ion_flux': vector_ion_flux,
    'vector_flux_labls': vector_flux_labls,
}

for name, var in variables.items():
    print(f"\n{name.upper()}:")
    print(f"  Tipo: {type(var)}")
    try:
        print(f"  Shape: {np.shape(var)}")
    except:
        print("  Shape: não aplicável")
    
    # imprimir os primeiros elementos
    if isinstance(var, (np.ndarray, list)):
        print("  Primeiros valores:", np.array(var).flatten()[:5])
    else:
        print("  Valor:", var)


# -------------------- CONVERT EPOCH TO DATETIME --------------------
# eixo temporal completo do SST
epoch0_dt = pd.to_datetime(epoch0, unit='s')
delta_td = pd.to_timedelta(delta_time, unit='s')

time_sst = epoch0_dt + delta_td

print(time_sst[0], time_sst[-1])

print()
print(type(epoch0))
print(epoch0.shape)
print()

#%% -------------------- DATA FRAME time series --------------------

df_thed_SST_series = pd.DataFrame(
    vector_ion_flux,
    index=pd.DatetimeIndex(time_sst),
    columns=vector_flux_labls
)

df_thed_SST_series.index.name = 'time'


# time_sst é o eixo temporal
# energy_bins = energia central de cada canal
# canais de energia (fixos no tempo)
energy_bins = energy[0, :]  # (16,)

df_spec_series = pd.DataFrame(
    ion_flux,
    index=pd.DatetimeIndex(time_sst),
    columns=energy_bins
)

df_spec_series.index.name = 'time'

energy_bins = energy[0, :]  # (16,)

print("SST energy channels (index : energy [keV])")
for i, e in enumerate(energy_bins):
    if np.isfinite(e):
        print(f"Channel {i:02d}: {e/1e3:.2f} keV")
    else:
        print(f"Channel {i:02d}: NaN (undefined)")



#%%# -------------------- ENERGY SPEC --------------------

df_spec = pd.DataFrame(ion_flux, index=time_sst, columns=energy_bins)
df_spec_interp = df_spec.interpolate(method='time')

Z_plot = df_spec_interp.T.values
Z_plot = np.ma.masked_invalid(Z_plot)
Z_plot = np.ma.masked_less_equal(Z_plot, 0)

time_plot = df_spec_interp.index.to_numpy()

energy_plot = np.array(energy_bins, dtype=float)
valid_energy = np.isfinite(energy_plot) & (energy_plot > 0)

energy_plot = energy_plot[valid_energy]
Z_plot = Z_plot[valid_energy, :]

Z_pos = Z_plot.compressed()

if Z_pos.size == 0:
    raise ValueError("Z_plot não contém valores positivos.")

vmin = np.percentile(Z_pos, 1)
vmax = np.percentile(Z_pos, 99)

# checar
print(np.all(np.diff(time_plot.astype('datetime64[ns]').astype(float)) > 0))
print("NaN em time_plot:", np.any(~np.isfinite(time_plot)))
print("NaN em energy_plot:", np.any(~np.isfinite(energy_plot)))
print("Shapes:", Z_plot.shape, time_plot.shape, energy_plot.shape)


# %% ---------------------------- 6. FILTER TIME INTERVAL -----------
df_plot = df_thed_FGM[
    (df_thed_FGM.index >= start_plot_dt) &
    (df_thed_FGM.index <= end_plot_dt)
]
#epoch = df_plot.index

df_plot = df_plot.reset_index()

df_plot["epoch"] = pd.to_datetime(df_plot["time"])
df_plot.set_index("epoch", inplace=True)

print()
print(type(df_plot))
print(df_plot.shape)
print(df_plot.head(2))
print()
print(df_plot.tail(2))

print()
print(df_plot.info())



# %% ---------------------------- 7. PLOTS ---------------------------

fig = plt.figure(figsize=(12, 15))

gs = GridSpec(
    nrows=5,
    ncols=2,
    figure=fig,
    width_ratios=[20, 1],      # coluna do colorbar
    height_ratios=[1, 1, 1, 1, 1],
    hspace=0.15,
    wspace=0.05
)

# Eixos principais
ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[1, 0], sharex=ax_a)
ax_c = fig.add_subplot(gs[2, 0], sharex=ax_a)
ax_d = fig.add_subplot(gs[3, 0], sharex=ax_a)
ax_e = fig.add_subplot(gs[4, 0], sharex=ax_a)

# Eixo dedicado ao colorbar do painel (c)
cax = fig.add_subplot(gs[2, 1])


# -------------------- (a) |B| --------------------
ax_a.plot(df_plot.index, df_plot['Btotal'], color='black')
ax_a.set_ylabel('|B| (nT)', fontsize=12)
ax_a.grid(True)
ax_a.tick_params(axis='both', labelsize=12)
ax_a.text(0.02, 0.85, '(a)', transform=ax_a.transAxes,
          fontsize=12, fontweight='bold')
ax_a.set_xlim(df_plot.index.min(), df_plot.index.max())


# -------------------- (b) B components (GSE) --------------------
ax_b.plot(df_plot.index, df_plot['Bx'], color='tab:red',   label=r'$B_x$ (GSE)')
ax_b.plot(df_plot.index, df_plot['By'], color='tab:blue',  label=r'$B_y$ (GSE)')
ax_b.plot(df_plot.index, df_plot['Bz'], color='tab:green', label=r'$B_z$ (GSE)')
ax_b.set_ylabel('B (nT)', fontsize=12)
ax_b.grid(True)
ax_b.legend(
    fontsize=12,
    loc='upper left',
    bbox_to_anchor=(1.02, 1),
    borderaxespad=0
)
ax_b.tick_params(axis='both', labelsize=12)
ax_b.text(0.02, 0.85, '(b)', transform=ax_b.transAxes,
          fontsize=12, fontweight='bold')


# -------------------- (c) Ion energy spectrogram --------------------

nan_fraction_time = np.mean(np.all(~np.isfinite(Z_plot), axis=0))

im = ax_c.pcolormesh(
    time_plot,
    energy_plot,
    Z_plot,
    shading='nearest',
    cmap='plasma',
    norm=LogNorm(vmin=vmin, vmax=vmax)
)


ax_c.set_yscale('log')
ax_c.set_ylabel('Ion energy (keV)', fontsize=12)
ax_c.tick_params(axis='both', labelsize=12)
ax_c.text(0.02, 0.85, '(c)', transform=ax_c.transAxes,
          fontsize=12, fontweight='bold')

# Colorbar do painel (c)
cbar = fig.colorbar(im, cax=cax)
cbar.set_label(
    r'(eV cm$^{-2}$ s$^{-1}$ sr$^{-1}$ eV$^{-1}$)',
    fontsize=12
)
cbar.ax.tick_params(labelsize=12)


# -------------------- (d) Ion flux vector (integrated) --------------------
flux_labels = ['Fx', 'Fy', 'Fz']
flux_colors = ['red', 'blue', 'green']

for i, (lab, col) in enumerate(zip(flux_labels, flux_colors)):
    ax_d.plot(time_sst, vector_ion_flux[:, i],
              color=col, label=lab)

#ax_d.set_ylabel('Ion flux vector (#/cm$^{2}$·s)', fontsize=12)
ax_d.set_ylabel('Ion flux vector', fontsize=12)
ax_d.grid(True)
ax_d.legend(
    fontsize=12,
    loc='upper left',
    bbox_to_anchor=(1.02, 1),
    borderaxespad=0
)
ax_d.set_ylim(4e4, 4e6)   # exemplo

ax_d.tick_params(axis='both', labelsize=12)
ax_d.text(0.02, 0.85, '(d)', transform=ax_d.transAxes,
          fontsize=12, fontweight='bold')



# -------------------- (e) Ion differential energy flux (series) --------------------
channels_to_plot = [0, 1, 2, 3, 4, 5]
colors = ['blue', 'green', 'orange', 'red', 'purple', 'brown']

for ch, c in zip(channels_to_plot, colors):
    ax_e.plot(
        df_spec_series.index,
        df_spec_series.iloc[:, ch],
        color=c,
        label=f'{energy_bins[ch]/1e3:.1f} keV'
    )

ax_e.set_yscale('log')
#ax_e.set_ylabel('Ion flux (eV cm$^{-2}$ s$^{-1}$ sr$^{-1}$ eV$^{-1}$)', fontsize=12)
ax_e.set_ylabel('Ion flux', fontsize=12)
ax_e.set_xlabel('Time (UT)', fontsize=12)
ax_e.grid(True, which='both', ls='--', alpha=0.5)
ax_e.legend(
    fontsize=12,
    loc='upper left',
    bbox_to_anchor=(1.02, 1),
    borderaxespad=0
)
ax_e.tick_params(axis='both', labelsize=11)
ax_e.text(0.02, 0.85, '(e)', transform=ax_e.transAxes,
          fontsize=12, fontweight='bold')


# -------------------- Ajustes finais --------------------
x_label_pos = -0.05

for ax in [ax_a, ax_b, ax_c, ax_d, ax_e]:
    ax.yaxis.set_label_coords(x_label_pos, 0.5)


for ax in [ax_a, ax_b, ax_c, ax_d]:
    plt.setp(ax.get_xticklabels(), visible=False)

fig.autofmt_xdate()
plt.show()




# %% SPEC GRID TEST
from scipy.interpolate import griddata
from matplotlib.colors import LogNorm
import numpy as np
import matplotlib.pyplot as plt

# --- 1. Converter tempo para float (dias desde 1970, sem timezone) ---
t0 = np.datetime64('1970-01-01T00:00:00')
time_float = (np.array(time_sst) - t0) / np.timedelta64(1, 'D')

# --- 2. Selecionar canais de energia válidos ---
energy_valid = energy_bins[np.isfinite(energy_bins)]

# --- 3. Selecionar valores válidos de fluxo ---
valid = np.isfinite(ion_flux)

# Flatten para griddata
time_flat = np.repeat(time_float, ion_flux.shape[1])[valid.flatten()]
energy_flat = np.tile(energy_bins, ion_flux.shape[0])[valid.flatten()]
flux_flat = ion_flux.flatten()[valid.flatten()]

# --- 4. Grid alvo ---
# Tempo regular (0.01 h)
dt_hours = 0.01
dt_days = dt_hours / 24
time_grid = np.arange(time_float.min(), time_float.max(), dt_days)

# Energia regular (log-spaced, fisicamente mais adequado)
energy_grid = np.logspace(
    np.log10(energy_valid.min()),
    np.log10(energy_valid.max()),
    64
)

# Meshgrid 2D
T_grid, E_grid = np.meshgrid(time_grid, energy_grid)

# --- 5. Interpolação linear ---
flux_grid = griddata(
    points=np.column_stack([time_flat, energy_flat]),
    values=flux_flat,
    xi=(T_grid, E_grid),
    method='linear',
    fill_value=np.nan
)

# --- 6. Plot ---
fig, ax = plt.subplots(figsize=(12, 4))

im = ax.pcolormesh(
    T_grid,
    E_grid / 1e3,   # energia em keV
    flux_grid,
    shading='nearest',
    cmap='plasma',
    norm=LogNorm()
)

ax.set_yscale('log')
ax.set_ylabel('Ion energy (keV)')
ax.set_xlabel('Time (days since 1970)')

cbar = fig.colorbar(im, ax=ax)
cbar.set_label('Ion flux (eV cm$^{-2}$ s$^{-1}$ sr$^{-1}$ eV$^{-1}$)')

plt.tight_layout()
plt.show()

# %%


#%%
# %% -------------------- PANEL (c): 2D INTERPOLATION (time × energy) --------------------

# ---------- 1. Preparar eixos ----------
# tempo em float (segundos desde epoch) para o griddata
time_float = time_plot.astype('datetime64[s]').astype(float)

energy_plot = np.array(energy_plot, dtype=float)

# ---------- 2. Grade original ----------
T, E = np.meshgrid(time_float, energy_plot)

Z = Z_plot.filled(np.nan)

# ---------- 3. Pontos válidos ----------
mask = np.isfinite(Z)

points = np.column_stack((T[mask], E[mask]))
values = Z[mask]

# ---------- 4. Nova grade (mesma resolução) ----------
Ti, Ei = np.meshgrid(time_float, energy_plot)

# ---------- 5. Interpolação 2D controlada ----------
Z_interp = griddata(
    points,
    values,
    (Ti, Ei),
    method='linear'
)

# mascarar valores inválidos ou não físicos
Z_interp = np.ma.masked_invalid(Z_interp)
Z_interp = np.ma.masked_less_equal(Z_interp, 0)

# ---------- 6. Recalcular limites do LogNorm ----------
Z_pos = Z_interp.compressed()
vmin = np.percentile(Z_pos, 1)
vmax = np.percentile(Z_pos, 99)

# ---------- 7. Plot ----------
# ---------- 7. Plot (figura separada) ----------
fig, ax = plt.subplots(figsize=(12, 4))

im = ax.pcolormesh(
    time_plot,
    energy_plot,
    Z_interp,
    shading='auto',
    cmap='plasma',
    norm=LogNorm(vmin=vmin, vmax=vmax)
)

ax.set_yscale('log')
ax.set_ylabel('Ion energy (keV)')
ax.set_xlabel('Time (UT)')

cbar = fig.colorbar(im, ax=ax)
cbar.set_label(
    r'Ion flux (eV cm$^{-2}$ s$^{-1}$ sr$^{-1}$ eV$^{-1}$)'
)

ax.tick_params(axis='both', labelsize=11)

plt.tight_layout()
plt.show()

# %%
