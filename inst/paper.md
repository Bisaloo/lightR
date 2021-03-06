---
title: "lightr: import spectral data and metadata in R"
tags:
  - R
  - colour
  - spectrophotometry
authors: 
  - name: Hugo Gruson
    orcid: "0000-0001-8493-9450"
    affiliation: 1
  - name: Thomas White
    orcid: "0000-0002-3976-1734"
    affiliation: 2
  - name: Rafael Maia
    orcid: "0000-0002-7563-9795"
    affiliation: 3
affiliations:
  - name: CEFE, University Montpellier, CNRS, University Paul Valery Montpellier 3, EPHE, IRD, Montpellier, France
    index: 1
  - name: School of Life and Environmental Sciences, University of Sidney, Sidney, NSW, Australia
    index: 2
  - name: Department of Ecology, Evolution and Environmental Biology, Columbia University, New York, NY
    index: 3
date: "October 24, 2019"
bibliography: lightr.bib
---

# Summary

Living organisms wildly differ in their ability to see colours [@Osorio2008].
For this reason, colour science relies on the use of objective measurements of
reflectance, transmittance, or aborbance spectra rather than human vision
[@Bennett1994;@Cuthill1999;@Eaton2005]. These spectra are then used in vision
models that allow scientists to predict how a given object is seen through the
eyes of a given species (e.g., how a male bird is seen by a potential mate).
This is the basis of all studies in for example the study of the evolution of
colours of animals and plants as communication signals.

Spectrometers record the amount of captured photons at different wavelengths
(usually between 300-700 nm for colour science, as many species are sensitive to
ultraviolet radiations). However, there is no standard file format for
spectrometry data and different scientific instrumentation companies use wildly
different formats to store spectral data. This use of non-standard file formats
jeopardises scientific reproducibility [@Peng2009] as other researchers might
not have the (paid) tools to open these files, and it makes us dependent on a
third-party which might vanish anytime, leaving a trove of scientific data
impossible to access. Vendors' proprietary software sometimes have an option to
convert those formats into human readable files such as `csv` but such software
are often expensive and they discard most metadata in the process. Yet, those
metadata are critical to ensure reproducibility of the measurements, and
ultimately of the scientific findings [@White2015].

In this article, we present `lightr`, an R package that aims at offering a
unified user-friendly interface for users to read reflectance, transmittance,
and absorbance spectra files from various formats in a single line of code.
Additionally, it provides for the first time a fully free and open source
solution to read proprietary spectra file formats on all operating systems.

`lightr` started as a fork from the popular R package `pavo`, which provides a
large suite of colour analysis tools [@Maia2013;@Maia2019].

# Package design

`lightr` has been designed to provide two levels on the complexity / 
customability trade-off:

* Spectral data and metadata for each file format are extracted using 
specialized parsers. Parsers are also aliased with many different names so that
users can often use `lr_parse_$extension()` where `$extension` is the file 
extension of the file to parse. For convenience, we also provide a generic 
fallback, named `lr_parse_generic()` that works for many "simple" formats, often
derived from `csv` or `tsv`. Specialized parsers should usually be preferred to
`lr_parse_generic()` because `lr_parse_generic()` is not able to parse metadata.

* Because spectrometers store each measurement in a separate file, the number of
files for a single study can quickly increase. To ensure easy and efficient
processing of those files, `lightr` also provides three high-levels functions
that can recursively find files and process them with a parallelized loop using
the `parallel` R package: `lr_get_spec()` and `lr_get_metadata()`, which import
respectively spectral data and metadata as `data.frame` in R, as well as
`lr_convert_tocsv()`, which converts all spectra files in a given folder as
`csv`, with the same filename (minus the file extension).

```
library(lightr)
lr_convert_tocsv(where = "yourfolder", ext = "ProcSpec")
```

# Recommended workflow

As mentioned earlier, proprietary spectrometry software can also export 
spectral data into a human-readable format (usually a kind of tabulation
separated values, or `tsv`, with a complex header). `lightr` can read files 
generated by this export step. We however
**do not recommend you use the software's built-in export function**,
because it will apply possibly unwanted transformation to your data
(interpolation and  subsetting) and may discard important metadata.

Instead, we recommend you keep the files in the proprietary format (such as
Avantes `ABS`, `ROH` and `TRM`, or OceanOptics `ProcSpec` and `jdx`) and that
you use `lightr` to convert them into your preferred file format (such as 
`csv`). 

# Usage and future directions

`lightr` can serve as a basis for colour analysis R packages to deal with the
file import step. Most of them can only read a limited variety of file
formats currently. Future versions of `pavo`, for example, will include `lightr`
as a dependency. Below is an illustration of a workflow where `lightr` is used
to import the spectral data, which is then analysed with `pavo`:

```
library(lightr)
specs <- lr_get_spec(where = "yourfolder", ext = "ProcSpec")

library(pavo)
plot(specs, col = spec2rgb(specs))
```

![](specs_joss.pdf)

```
summary(specs, subset = TRUE)
```

|       B2|       S8|  H1|
|--------:|--------:|---:|
|  9.31682| 1.915661| 561|
| 11.26643| 2.156246| 551|
| 12.78053| 2.128401| 557|
| 13.41558| 2.123076| 551|
| 13.44852| 2.118632| 562|
| 12.14774| 2.210931| 557|
| 11.76633| 2.076845| 547|
| 10.62519| 2.204452| 551|
| 10.14280| 2.272771| 547|

The first column indicates the brightness (in % relative to a white reference),
the second is the saturation (also called spectral purity) and the last contains
the hue (in nm).

`lightr` can also prove useful for developers of 
other programming languages, providing a free and open source template that can
easily be translated to such other languages.
We also plan on providing a web application based on shiny 
(<https://github.com/rstudio/shiny>), which uses `lightr` in the background, and
provides users with limited R or technical knowledge with a simple and 
convenient way to convert all their proprietary files to `csv`.

# Acknowledgements

We thank the two rOpenSci reviewers, Jeroen Ooms and Karthik Ram, for their 
helpful feedback that improved this package, as well as JOSS editor Daniel S. 
Katz for his comments on this manuscript.

# References
