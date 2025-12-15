# Data Sources for Tamaraw Project

## 1. Land Cover Data (ESA WorldCover 2021)

**Dataset:** ESA WorldCover 10-meter 2021 Land Cover  
**Provider:** European Space Agency (ESA) / Copernicus  
**Data Access:** https://esa-worldcover.org/en/data-access  
**Format:** GeoTIFF  
**Resolution:** 10 meters  

### Acquisition Notes
- Accessed via the ESA WorldCover Data Access Portal.  
- Located and downloaded the tile covering Palawan, Philippines (e.g., `N09E117`).  
- Used in the project to generate forest and mangrove masks and to downsample to ~100 m resolution for simulation.

---

## 2. Climate Data (ERA5 Reanalysis)

**Dataset:** ERA5 Hourly Data on Single Levels (1940–present)  
**Provider:** Copernicus Climate Data Store (ECMWF)  
**Data Access:** https://cds.climate.copernicus.eu  
**Format:** NetCDF  
**Resolution:** ~0.25° (~27 km)  
**Temporal Resolution:** Hourly  

### Acquisition Notes
- Searched for “ERA5 single levels” on the CDS portal.  
- Selected the dataset *ERA5 hourly data on single levels from 1940 to present*.  
- Configured a geographic subset for Palawan using the bounding box:  
  - North: 12.0  
  - West: 117.0  
  - South: 8.0  
  - East: 120.0  
- Downloaded key climate variables including temperature, dewpoint, precipitation, wind components, and surface pressure.

---

_Last updated: 2025-12-14_
