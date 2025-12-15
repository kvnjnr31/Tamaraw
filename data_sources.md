# Data Sources for Tamaraw Project

## 1. Land Cover Data (ESA WorldCover 2021)

**Dataset:** ESA WorldCover 10-meter 2021 Land Cover  
**Provider:** European Space Agency (ESA) / Copernicus  
**Data Access:** https://esa-worldcover.org/en/data-access  
**Format:** GeoTIFF  
**Resolution:** 10 meters  

### Acquisition Notes
- Went to the ESA WorldCover Data Access Portal.  
- Used the interactive map to locate the Palawan tile (`N09E117`).  
- Downloaded the 10‑m GeoTIFF and later downsampled it to ~100 m for simulation.  

---

## 2. Climate Data (ERA5 Reanalysis)

**Dataset:** ERA5 Hourly Data on Single Levels (1940–present)  
**Provider:** Copernicus Climate Data Store (ECMWF)  
**Data Access:** https://cds.climate.copernicus.eu  
**Format:** NetCDF  
**Resolution:** ~0.25° (~27 km)  
**Temporal Resolution:** Hourly  

### Step‑by‑Step Acquisition Breakdown

#### **Step 1 — Navigate to the ERA5 Single‑Levels Page**
Screenshot:  
![ERA5 Overview](/mnt/data/80322c83-2601-4c7b-8b61-8e009122cd48.png)

- From the Climate Data Store homepage, searched for **“ERA5 single levels”**.
- Selected **ERA5 hourly data on single levels from 1940 to present**.

---

#### **Step 2 — Open the Download Tab and Select Variables**
Screenshot:  
![ERA5 Variable Selection](/mnt/data/299e52d6-3c34-43a7-b18d-65a38c325f28.png)

Variables selected:

- **10m u-component of wind**  
- **10m v-component of wind**  
- **2m dewpoint temperature**  
- **2m temperature**  
- **Surface pressure**  
- **Total precipitation**

These are the core climate drivers for mosquito diffusion and habitat suitability modeling.

---

#### **Step 3 — Select Year, Months, Days, and Hours**

- Year: *e.g., 2025*  
- All months  
- All days  
- All 24 hourly timesteps  

This provides a complete annual climate dataset at hourly resolution.

---

#### **Step 4 — Choose the Geographic Sub‑Region**

Bounding box used:

```
North: 12
West: 117
South: 8
East: 120
```

This frame isolates the Palawan region.

---

#### **Step 5 — Choose Output Format and Submit**

- Output format: **NetCDF (.nc)**  
- Pressed **Submit Form**  
- File was generated under *Your Requests* and downloaded as a .zip containing the NetCDF file.

---

### Summary
The ERA5 dataset was filtered to include only climate variables relevant to vector‑borne disease dynamics in Palawan. Hourly data were extracted for the entire year, clipped to the Palawan bounding box, and saved as a NetCDF file for downstream processing.

---

_Last updated: 2025‑12‑14_
