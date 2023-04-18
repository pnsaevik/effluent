---
title: 'Effluent: A Python package for modelling effluent discharge'
tags:
  - Python
  - wastewater
  - turbulence
  - galactic dynamics
  - milky way
authors:
  - name: Pål Næverlid Sævik
    orcid: 0000-0002-7301-2008
    equal-contrib: true
    affiliation: 1

affiliations:
 - name: Institute of Marine Research, Norway
   index: 1

date: 18 April 2023
bibliography: paper.bib

---

# Summary

Effluent dispersion modelling is an essential tool in the management of water
resources. Both domestic and industrial wastewater carry substances that may
be harmful to the environment unless sufficiently diluted. Also, 
wastewater containing nutrients may lead to eutrophication of the receiving
water body if the outfall rises to the surface. Computer modelling of the
dilution process may help discover this type of problems before they appear,
and can guide the design of an outfall system to minimize the impact to the
environment.

# Statement of need

`effluent` is an open-source python package for simulating the dispersion of
effluent discharges from wastewater pipes. The underlying model is based on
`@lee:2003`. The implementation also contains couplings to the popular
open-source ocean model ROMS `[@shchepetkin:2005]`.

There already exists multiple closed-source applications for modelling outfall
dispersion, such as CORMIX `[@jirka:2004]`, VISUAL PLUMES `[@frick:2004]`
or VISJET `[@lee:2003]`. These are well-suited for governmental and industrial
applications. In research, however, it is desirable to have an open-source,
transparent implementation which can also be adapted and tailored to the
specific needs of a reserach project. The package `effluent` provides this,
along with a detailed online documentation with usage examples.

# Acknowledgements

This work is financed by the Institute of Marine Sciences, Norway.

# References
