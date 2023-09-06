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
or VISJET [@Lee2003]. A thorough comparison between these models has been done
by Palomar et al. [@Palomar2012]. Since `effluent` has the same theoretical
foundation as VISJET, the capabilities of the two models are similar. The main
difference is that `effluent` does not contain any internal visualization
capabilities, and results must be visualized using external packages. On the
other hand, `effluent` is easier to incorporate into a scripting environment
as its input and output formats are standardized and well documented. In
addition, `effluent` allows ambient co- and crossflow currents to vary in both
time and depth. This makes it possible to combine with modelled ambient
current data. Being open source, `effluent` also allows users to fork and
modify the program source code to suit specific needs.

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
  
- id: Palomar2012
  accessed:
    - year: 2023
      month: 9
      day: 6
  author:
    - family: Palomar
      given: P.
    - family: Lara
      given: J. L.
    - family: Losada
      given: I. J.
    - family: Rodrigo
      given: M.
    - family: Alvárez
      given: A.
  citation-key: Palomar2012
  container-title: Desalination
  container-title-short: Desalination
  DOI: 10.1016/j.desal.2011.11.037
  ISSN: 0011-9164
  issued:
    - year: 2012
      month: 3
      day: 30
  page: 14-27
  source: ScienceDirect
  title: 'Near field brine discharge modelling part 1: Analysis of commercial tools'
  title-short: Near field brine discharge modelling part 1
  type: article-journal
  URL: https://www.sciencedirect.com/science/article/pii/S0011916411009702
  volume: '290'
...

