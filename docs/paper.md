---
title: 'Effluent: A Python package for modelling effluent discharge'
tags:
  - Python
  - wastewater
  - turbulence
  - outfalls
  - plumes

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
and can guide the design of an outfall system to minimize environmental
impact.

# Statement of need

`effluent` is an open-source python package for simulating the dispersion of
effluent discharges from wastewater pipes. The underlying model is based on
[@Lee2003]. The implementation also contains couplings to the popular
open-source ocean model ROMS [@Shchepetkin2005].

There already exists multiple closed-source applications for modelling outfall
dispersion, such as CORMIX [@Jirka2004], VISUAL PLUMES [@Frick2004]
or VISJET [@Lee2003]. These are well-suited for governmental and industrial
applications. In research, however, it is desirable to have an open-source,
transparent implementation which can also be adapted and tailored to the
specific needs of a research project. The package `effluent` provides this,
along with a detailed online documentation with usage examples.

# Acknowledgements

This work is financed by the Institute of Marine Sciences, Norway.

# References

---
references:
- id: Frick2004
  author:
    - family: Frick
      given: Walter E.
  citation-key: Frick2004
  accessed:
    - year: 2021
      month: 2
      day: 19
  collection-title: Facilitating Constructive Government-Industry Interations
  container-title: Environmental Modelling & Software
  container-title-short: Environmental Modelling & Software
  DOI: 10.1016/j.envsoft.2003.08.018
  ISSN: 1364-8152
  issue: '7'
  issued:
    - year: 2004
      month: 7
      day: 1
  language: en
  page: 645-654
  source: ScienceDirect
  title: Visual Plumes mixing zone modeling software
  type: article-journal
  URL: https://www.sciencedirect.com/science/article/pii/S1364815203001890
  volume: '19'

- id: Jirka2004
  accessed:
    - year: 2021
      month: 2
      day: 22
  author:
    - family: Jirka
      given: Gerhard H.
  citation-key: Jirka2004
  container-title: Environmental Fluid Mechanics
  container-title-short: Environmental Fluid Mechanics
  DOI: 10.1023/A:1025583110842
  ISSN: 1573-1510
  issue: '1'
  issued:
    - year: 2004
      month: 3
      day: 1
  language: en
  page: 1-56
  source: Springer Link
  title: >-
    Integral Model for Turbulent Buoyant Jets in Unbounded Stratified Flows.
    Part I: Single Round Jet
  title-short: >-
    Integral Model for Turbulent Buoyant Jets in Unbounded Stratified Flows.
    Part I
  type: article-journal
  URL: https://doi.org/10.1023/A:1025583110842
  volume: '4'

- id: Lee2003
  accessed:
    - year: 2021
      month: 2
      day: 19
  author:
    - family: Lee
      given: Joseph H. W.
    - family: Chu
      given: Vincent H.
  citation-key: Lee2003
  DOI: 10.1007/978-1-4615-0407-8
  event-place: Boston, MA
  ISBN: 978-1-4613-5061-3 978-1-4615-0407-8
  issued:
    - year: 2003
  language: en
  publisher: Springer US
  publisher-place: Boston, MA
  source: DOI.org (Crossref)
  title: Turbulent Jets and Plumes
  type: book
  URL: http://link.springer.com/10.1007/978-1-4615-0407-8

- id: Shchepetkin2005
  accessed:
    - year: 2019
      month: 4
      day: 8
  author:
    - family: Shchepetkin
      given: Alexander F.
    - family: McWilliams
      given: James C.
  citation-key: Shchepetkin2005
  container-title: Ocean Modelling
  DOI: 10.1016/j.ocemod.2004.08.002
  ISSN: '14635003'
  issue: '4'
  issued:
    - year: 2005
      month: 1
  language: en
  page: 347-404
  source: DOI.org (Crossref)
  title: >-
    The regional oceanic modeling system (ROMS): a split-explicit, free-surface,
    topography-following-coordinate oceanic model
  title-short: The regional oceanic modeling system (ROMS)
  type: article-journal
  URL: https://linkinghub.elsevier.com/retrieve/pii/S1463500304000484
  volume: '9'
...

