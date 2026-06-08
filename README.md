# Reflection fraction

Calculation of the reflection fraction for thick discs with a lamp-post corona

## Some relevant reading

### General papers

- [A light bending model for the X-ray temporal and spectral properties of accreting black holes](https://ui.adsabs.harvard.edu/abs/2004MNRAS.349.1435M/abstract) (2004) by Miniutti and Fabian. This paper that shows how the continuum and reflection flux (and hence reflection fraction) change as the corona height changes.
- [Exploring the Effects of Disk Thickness on the Black Hole Reflection Spectrum](https://ui.adsabs.harvard.edu/abs/2018ApJ...855..120T/abstract) (2018) by Taylor and Reynolds. This paper shows how the reflection spectrum changes with Eddington fraction / disc thickness.
- [Normalizing a relativistic model of X-ray reflection. Definition of the reflection fraction and its implementation in relxill](https://ui.adsabs.harvard.edu/abs/2016A%26A...590A..76D/abstract) (2016) by Dauser, Garcia, and Walton. This paper discusses the technicalities of defining the "reflection fraction".
- [Gradus.jl: spacetime-agnostic general relativistic ray-tracing for X-ray spectral modelling](https://ui.adsabs.harvard.edu/abs/2026MNRAS.545f1770B/abstract) (2026) by Baker and Young. This paper describes the relativistic ray tracing model [Gradus.jl](https://codeberg.org/astro-group/Gradus.jl) that we will use to perform the calculations we need.
- [General relativistic effects in black hole X-ray coronal models](https://cosroe.com/thesis.html) (2025) by Baker. This PhD thesis is an excellent pedagogical guide and includes the code used to produce the figures.

### Reflection modelling

- [Towards modelling X-ray reverberation in AGN: piecing together the extended corona](https://ui.adsabs.harvard.edu/abs/2016MNRAS.458..200W/abstract) by Wilkins et al. (2016). This paper shows 2D reverberation transfer functions convolved with XILLVER reflection spectra for a range of extended corona geometries.
- [Venturing beyond the ISCO: detecting X-ray emission from the plunging regions around black holes](https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.5532W/abstract) by Wilkins et al. (2020). This uses a similar technique to the previous paper but includes emission from the plunging region (we're not including this here but perhaps we could in the future). This paper has nice figures showing the calculated 2D reverberation transfer functions and those that are convolved with XILLVER reflection spectra.
- In the [Gradus paper](https://ui.adsabs.harvard.edu/abs/2026MNRAS.545f1770B/abstract) the routine [`reverberation.lag-energy.jl`](https://github.com/fjebaker/gradus-paper/blob/main/code/reverberation.lag-energy.jl) is used to calculate the lag-energy spectrum shown in Figure 13. This makes use of [`reverberation.jl`](https://codeberg.org/astro-group/Gradus.jl/src/branch/main/src/reverberation.jl) functions in Gradus.
