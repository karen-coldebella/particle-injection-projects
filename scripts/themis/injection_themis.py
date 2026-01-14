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
spec_ion_flux = cdf_SST_thd.varget('thd_psif_en_eflux')

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
series_ion_flux = cdf_SST_thd.varget('thd_psif_flux')

# Labels of the ion flux vector components
# Define the physical meaning of each column of thd_psif_flux
# e.g., instrumental directions or DSL components
flux_labls = cdf_SST_thd.varget('thd_psif_flux_labl')


# -------------------- INSPECTION --------------------
var_info = cdf_SST_thd.varattsget('thd_psif_time')
print("Informações da variável 'thd_psif_time':", var_info)


variables = {
    'epoch0': epoch0,
    'delta_time': delta_time,
    'spec_ion_flux': spec_ion_flux,
    'energy': energy,
    'series_ion_flux': series_ion_flux,
    'flux_labls': flux_labls,
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
    series_ion_flux,
    index=pd.DatetimeIndex(time_sst),
    columns=flux_labls
)

df_thed_SST_series.index.name = 'time'


# time_sst já é o eixo temporal
# energy_bins = energia central de cada canal
energy_bins = np.ravel(energy)

df_spec_series = pd.DataFrame(
    spec_ion_flux, 
    index=time_sst, 
    columns=energy_bins  # 16 canais de energia
)

df_spec_series.index.name = 'time'


# -------------------- ENERGY SPEC --------------------
#energy_bins = energy[0, :]      # (16,)
#time = pd.to_datetime(epoch0, unit='s')
#Z = spec_ion_flux.T             # (16, N)

# Cria DataFrame temporário com tempo como índice
df_spec = pd.DataFrame(spec_ion_flux, index=time_sst, columns=energy_bins)

# Interpolação linear ao longo do tempo
df_spec_interp = df_spec.interpolate(method='time')

# Prepara para pcolormesh
Z_plot = df_spec_interp.T.values  # shape (16, N)
time_plot = df_spec_interp.index
energy_plot = df_spec_interp.columns.values

time_plot = np.array(time_plot)
energy_plot = np.array(energy_bins)

energy_plot = pd.Series(energy_plot).interpolate().to_numpy()

# checar
print(np.any(np.isnan(time_plot)), np.any(np.isnan(energy_plot)))


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
fig, axs = plt.subplots(
    5, 1,                      # 4 painéis originais + 1 painel para séries de energia
    figsize=(12, 10), 
    sharex=True,
    gridspec_kw={'height_ratios': [1, 1, 1, 1, 1]}  # todos com mesma altura
)


# -------------------- (a) |B| --------------------
axs[0].plot(df_plot.index, df_plot['Btotal'], color='black')
axs[0].set_ylabel('|B| (nT)', fontsize=13)
axs[0].grid(True)
axs[0].tick_params(axis='both', labelsize=11)
axs[0].text(0.02, 0.85, '(a)', transform=axs[0].transAxes,
            fontsize=13, fontweight='bold')
axs[0].set_xlim(np.min(df_plot.index), np.max(df_plot.index))


# -------------------- (b) B components (GSE) --------------------
axs[1].plot(df_plot.index, df_plot['Bx'], label='Bx (GSE)')
axs[1].plot(df_plot.index, df_plot['By'], label='By (GSE)')
axs[1].plot(df_plot.index, df_plot['Bz'], label='Bz (GSE)')
axs[1].set_ylabel('B (nT)', fontsize=13)
axs[1].grid(True)
axs[1].legend(
    fontsize=11,
    loc='upper left',
    bbox_to_anchor=(1.02, 1),
    borderaxespad=0
)
axs[1].tick_params(axis='both', labelsize=11)
axs[1].text(0.02, 0.85, '(b)', transform=axs[1].transAxes,
            fontsize=13, fontweight='bold')
axs[1].set_xlim(np.min(df_plot.index), np.max(df_plot.index))


# -------------------- (c) Ion energy spectrogram --------------------
im = axs[2].pcolormesh(
    time_plot,                 
    energy_plot / 1e3,         
    Z_plot,                    
    shading='auto',
    cmap='plasma',
    norm=plt.matplotlib.colors.LogNorm(
        vmin=np.nanmin(Z_plot[Z_plot>0]),   
        vmax=np.nanmax(Z_plot)
    )
)
axs[2].set_ylabel('Energy (keV)', fontsize=13)
axs[2].set_yscale('log')
axs[2].tick_params(axis='both', labelsize=11)
axs[2].text(0.02, 0.85, '(c)', transform=axs[2].transAxes,
            fontsize=13, fontweight='bold')
fig.colorbar(im, ax=axs[2], label='Flux (eV/cm²·s·sr·eV)')


# -------------------- (d) Ion flux series (integrated) --------------------
for i, label in enumerate(flux_labls):
    axs[3].plot(time_sst, series_ion_flux[:, i], label=label)
axs[3].set_ylabel('Ion flux (#/cm²·s)', fontsize=13)
axs[3].grid(True)
axs[3].legend(
    fontsize=11,
    loc='upper left',
    bbox_to_anchor=(1.02, 1),
    borderaxespad=0
)
axs[3].tick_params(axis='both', labelsize=11)
axs[3].text(0.02, 0.85, '(d)', transform=axs[3].transAxes,
            fontsize=13, fontweight='bold')


# -------------------- (e) Ion differential energy flux (séries) --------------------
channels_to_plot = [0, 4, 8, 12, 15]  # índices dos canais
colors = ['blue', 'green', 'orange', 'red', 'purple']

for ch, c in zip(channels_to_plot, colors):
    axs[4].plot(df_spec_series.index, df_spec_series.iloc[:, ch], 
                label=f'{energy_bins[0][ch]/1e3:.1f} keV', color=c)

axs[4].set_yscale('log')
axs[4].set_ylabel('Ion flux (eV/cm²·s·sr·eV)', fontsize=13)
axs[4].set_xlabel('Time (UT)', fontsize=13)
axs[4].grid(True, which='both', ls='--', alpha=0.5)
axs[4].legend(
    fontsize=11,
    loc='upper left',
    bbox_to_anchor=(1.02, 1),
    borderaxespad=0
)
axs[4].tick_params(axis='both', labelsize=11)
axs[4].text(0.02, 0.85, '(e)', transform=axs[4].transAxes,
            fontsize=13, fontweight='bold')


# -------------------- Ajuste de formato de data --------------------
plt.gcf().autofmt_xdate()
plt.tight_layout()
plt.show()


# %%
