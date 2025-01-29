# Gamma Shielding Material Analysis: Lead, HY80, and 30MnB5

This project aims to compare the gamma shielding properties of three materials—Lead, HY80, and 30MnB5—for potential applications in military ships. The analysis involves gamma spectroscopy to evaluate key parameters such as peak intensity, full width at half maximum (FWHM), and Compton edge stability across varying material thicknesses.

---

## **Project Objectives**
1. Evaluate the gamma attenuation efficiency of Lead, HY80, and 30MnB5.
2. Analyze the scattering profiles using FWHM and Gaussian fits.
3. Investigate the stability of Compton edges to understand photon energy interactions.
4. Identify the most suitable material or hybrid shielding design for submarine applications.

---

## **Project Workflow**
1. **Data Preparation**:
   - ROOT files were generated for each material at thicknesses ranging from 0 cm to 25 cm.
   - Gamma spectra were extracted using ROOT’s TTree and TCanvas structures.

2. **Analysis**:
   - Key parameters such as peak intensity, FWHM, and Compton edges were calculated for each material and thickness.
   - Trends were plotted to highlight material performance.

3. **Visualization**:
   - Plots were created to compare peak intensity, FWHM, and Compton edges across materials and thicknesses.
   - Gaussian fits were applied to analyze scattering profiles.

4. **Recommendations**:
   - Proposed hybrid designs combining HY80 and 30MnB5 for optimal performance.
   - Emphasized weight savings and structural integrity for submarine applications.

---

## **Key Findings**
1. **Lead**:
   - Best gamma attenuation per unit thickness.
   - Narrow FWHM indicates minimal scattering.
   - Stable Compton edges across thicknesses.
   - Limitation: High density makes it impractical for weight-sensitive applications.

2. **HY80**:
   - Moderate attenuation efficiency.
   - Broader FWHM compared to Lead, indicating more scattering.
   - Slight shifts in Compton edges with increasing thickness.
   - Lightweight and structurally strong, making it viable for hybrid designs.

3. **30MnB5**:
   - Similar performance to HY80 with slightly better attenuation.
   - Increased scattering compared to Lead.
   - Potential for neutron shielding due to boron content.

---

## **Key Visualizations**
1. **Peak Intensity vs. Thickness**:
   - Demonstrates attenuation efficiency.
   - Lead shows the steepest decline in intensity.

2. **FWHM vs. Thickness**:
   - Highlights scattering contributions for each material.

3. **Compton Edge Stability**:
   - Indicates photon energy redistribution due to material interactions.

4. **Weight vs. Performance**:
   - Compares weight savings of HY80 and 30MnB5 relative to Lead.

---

## **How to Run**
1. Open Terminal

2. Clone the repository:
    ```bash
    git clone https://github.com/ammanielg/Final_Project.git
    cd Final_Project
    conda create --name root_env python=3.11 -y
    conda activate root_env
    pip install -r requirements.txt
    conda install -c conda-forge root -y
    code .

3. Navigate to notebooks/main.ipynb and open it.

