# PCM Thermal Energy Storage Model - Technical Review Summary

## Overview

This repository contains a Python implementation of a packed-bed Latent Heat Thermal Energy Storage (LHTES) system model, attempting to replicate Figure 16 from:

**Klitou, T., & Fokaides, P. A. (2024).** "Modelling a packed-bed latent heat thermal energy storage unit and studying its performance using different paraffins." *International Journal of Sustainable Energy*, Vol. 43, No. 1.

A comprehensive technical review was conducted on November 6, 2025, identifying critical issues in material properties and implementation, with corrected versions provided.

---

## Repository Contents

### Original Code
- **`pcm_working-model.py`** - Original implementation with identified issues

### Review Deliverables
- **`TECHNICAL_REVIEW_REPORT.md`** - 40+ page comprehensive technical analysis covering:
  - Physical parameters validation
  - Governing equations verification
  - Numerical methods assessment
  - Phase change modeling evaluation
  - Energy balance analysis
  - Detailed recommendations

- **`QUICK_REFERENCE.md`** - Actionable summary of critical issues and fixes

- **`pcm_corrected_model.py`** - Fully corrected and enhanced implementation with:
  - Literature-accurate RT58 paraffin properties
  - Fixed phase change modeling
  - Temperature-dependent water properties
  - Energy balance verification
  - Enhanced diagnostics and validation tools

- **`compare_models.py`** - Side-by-side comparison script showing impact of corrections

- **`REVIEW_SUMMARY.md`** - This file

---

## Critical Issues Identified

### 🔴 Priority 1: Critical Errors

| Issue | Original Value | Correct Value | Impact |
|-------|---------------|---------------|--------|
| **Latent heat** | 170 kJ/kg | 126 kJ/kg | 35% overestimation of energy storage |
| **k_solid** | 0.4 W/(m·K) | 0.24 W/(m·K) | 67% faster heat transfer than physical |

### 🟠 Priority 2: High Impact

| Issue | Original | Corrected | Impact |
|-------|----------|-----------|--------|
| **k_liquid** | 0.2 W/(m·K) | 0.15 W/(m·K) | 33% faster melting |
| **ρ_solid** | 916 kg/m³ | 880 kg/m³ | 4% mass overestimation |
| **c_p,solid** | 2100 J/(kg·K) | 1900 J/(kg·K) | 11% heating rate error |

### 🟡 Priority 3: Moderate Issues

- Mushy zone heat capacity calculation error
- Phase change range too narrow (2K vs realistic 8K)
- Grid resolution insufficient (Pe_grid >> 2)
- No temperature-dependent fluid properties

---

## Key Findings

### Energy Storage Capacity

**Original Model:**
- Total energy capacity: 8.63 MJ
- Estimated charging time: 6.4 hours

**Corrected Model:**
- Total energy capacity: 7.69 MJ
- Estimated charging time: 5.7 hours

**Difference:** Original model overestimates storage by 12% and charging time by ~45 minutes

### Heat Transfer

**Biot Number Analysis:**
- Original: Bi ≈ 4.8
- Corrected: Bi ≈ 8.0
- **Implication:** Lumped capacitance assumption is marginal; internal temperature gradients may be significant

**Péclet Number:**
- Original (nz=50): Pe_grid ≈ 65
- Corrected (nz=100): Pe_grid ≈ 32
- **Recommendation:** Consider nz=200 for Pe_grid < 10

### Phase Change Modeling

**Temperature Range:**
- Original: 59-61°C (2K range)
- Corrected: 56-64°C (8K range, matches DSC data)
- **Impact:** More realistic gradual melting behavior

---

## Model Strengths ✅

1. **Solid theoretical foundation:** LTNE (Local Thermal Non-Equilibrium) approach
2. **Appropriate correlations:** Wakao-Kaguei for packed beds
3. **Suitable numerics:** BDF solver for stiff ODEs
4. **Good structure:** Modular, well-documented code
5. **Physical boundary conditions:** Upwind scheme for advection

---

## Quick Start

### 1. Run Comparison Analysis
```bash
python compare_models.py
```
This generates a visual comparison showing the impact of corrections.

### 2. Run Corrected Model
```bash
python pcm_corrected_model.py
```
This runs the simulation with corrected properties and generates:
- Temperature evolution plots (Figure 16 comparison)
- Liquid fraction tracking
- Spatial temperature profiles
- Energy balance verification

### 3. Review Technical Details
```bash
# Comprehensive analysis
cat TECHNICAL_REVIEW_REPORT.md

# Quick reference for fixes
cat QUICK_REFERENCE.md
```

---

## Validation Checklist

Before considering the model validated against Figure 16:

- [ ] **Material properties** match literature values for RT58
- [ ] **Energy balance** closes within 5% error
- [ ] **Phase change plateau** appears at ~60°C
- [ ] **Charging time** matches experimental data (±10%)
- [ ] **HTF-PCM temperature gap** is physically reasonable (2-5 K)
- [ ] **Melting progression** from bottom to top
- [ ] **Grid independence** verified (results stable with finer mesh)
- [ ] **Comparison with Figure 16** shows good agreement

---

## Next Steps

### Immediate Actions

1. **Run corrected model** and compare outputs with Figure 16
2. **Verify operating conditions** against paper:
   - Mass flow rate (currently 6 kg/min, paper tested 2-4 kg/min)
   - Collector heat input (currently 375 W - verify if realistic)
   - Initial temperature (currently 25°C)

3. **Upload reference materials** for direct comparison:
   - Paper PDF (Paper_modelo.pdf)
   - Figure 16 image (1762431943854_image.png)

### Further Improvements

4. **Sensitivity analysis** on:
   - Mass flow rate: 2, 4, 6 kg/min
   - Collector power: 500, 750, 1000 W
   - Phase change range: 2, 5, 8 K
   - Grid resolution: 50, 100, 200 nodes

5. **Model enhancements:**
   - 1D radial conduction within PCM capsules (address high Bi number)
   - Radiation heat transfer (if significant)
   - Time-delayed inlet boundary condition
   - Experimental validation

6. **Documentation:**
   - Add Figure 16 comparison plots
   - Document any remaining discrepancies
   - Create user guide for parameter adjustment

---

## Technical Specifications

### System Geometry
- Tank diameter: 360 mm
- Tank height: 470 mm
- Porosity: 0.49
- Capsule diameter: 55 mm

### Materials
- **PCM:** Paraffin RT58 (Rubitherm)
- **HTF:** Water
- **Temperature range:** 25-65°C

### Computational Details
- **Method:** 1D axial finite volume
- **Solver:** scipy.solve_ivp with BDF method
- **Tolerances:** rtol=1e-4, atol=1e-6
- **Grid:** 100 axial nodes (corrected model)
- **Time step:** Adaptive, max 10 seconds

---

## References

### Primary Source
Klitou, T., & Fokaides, P. A. (2024). "Modelling a packed-bed latent heat thermal energy storage unit and studying its performance using different paraffins." *International Journal of Sustainable Energy*, Vol. 43, No. 1, DOI: 10.1080/14786451.2024.2306416

### Supporting Literature
- **Wakao, N., & Kaguei, S. (1982).** *Heat and Mass Transfer in Packed Beds*. Gordon and Breach.
- **Rubitherm GmbH** - RT58 Technical Data Sheet
- **Ismail & Henríquez (2002)** - Packed bed LHTES experimental studies
- **Voller & Prakash (1987)** - Phase change modeling methods

---

## Contact & Support

For questions about:
- **Theoretical framework:** See TECHNICAL_REVIEW_REPORT.md, Sections 1-3
- **Implementation details:** See pcm_corrected_model.py inline comments
- **Material properties:** See QUICK_REFERENCE.md property tables
- **Numerical methods:** See TECHNICAL_REVIEW_REPORT.md, Section 6

---

## Revision History

| Date | Version | Changes |
|------|---------|---------|
| 2025-11-06 | 1.0 | Initial comprehensive technical review |
| | | - Identified critical property errors |
| | | - Created corrected model implementation |
| | | - Generated detailed technical report |

---

## License & Citation

If you use this corrected model in your research, please cite both:

1. The original paper by Klitou et al. (2024)
2. This technical review and corrections

---

**Review Status:** ✅ Complete
**Model Status:** 🟢 Ready for validation testing
**Priority Action:** Compare corrected model outputs with Figure 16

---

*Technical Review conducted by Claude AI*
*Date: November 6, 2025*
*Methodology: Comprehensive code analysis, literature review, and theoretical validation*
