
# Hnoss
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

## Introduction

Hnoss is an extension of the Freyja<sup>[1](#references)</sup> SARS-CoV-2 strain deconvolution package<sup>[2](#references)</sup> to allow for easier manipulation and more flexibility. 

## Table of Contents

- [Introduction](#introduction)
- [Quick-Start Guide](#quick-start%guide)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Input](#input)
- [Output](#output)
- [Workflow](#workflow)
- [References](#references)

## Quick-Start Guide

Import into your script:
```
conda activate hnoss-env
import hnoss
```
For details on available arguments, enter:
```
hnoss --help
```

## Dependencies

[Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) is required to build an environment with required workflow dependencies. To create the environment
```
conda create -n hnoss-env
```
then add the following channels
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
and then install hnoss
```
conda install hnoss
```

See REQUIREMENTS.txt for package dependancies.

## Input
This tool runs after Freyja via a flexible import functions. For strain deconvolution, it allows any combination of output files from ```freyja demix``` and ```freyja aggregate```. Please refer to the [Freyja documentation](https://github.com/andersen-lab/Freyja).

## Output
#### <ins>Typical Hnoss structure:</ins>

<table>
<thead>
  <tr>
    <th colspan="3">SampleInfo</th>
    <th colspan="4"  bgcolor="gray"><font color="black">SummarizedLineages</font></th>
    <th colspan="4">RawLineages</th>
  </tr>
</thead>
<tbody>
  <tr>
    <th>file</th>
    <th>resid</th>
    <th>coverage</th>
    <th bgcolor="gray"><font color="black">Alpha</font></th>
    <th bgcolor="gray"><font color="black">Beta</font></th>
    <th bgcolor="gray"><font color="black">XBB* (XBB.X)</font></th>
    <th bgcolor="gray"><font color="black">...</font></th>
    <th>A.1</th>
    <th>BA.2.1</th>
    <th>XBB.1.11</th>
    <th>...</th>
  </tr>
  <tr>
    <td>230101_1_S1_L001.variants.tsv</td>
    <td>8.28</td>
    <td>98.475</td>
    <td bgcolor="gray"><font color="black">0.121</font></td>
    <td bgcolor="gray"><font color="black">0.312</font></td>
    <td bgcolor="gray"><font color="black">0.426</font></td>
    <td bgcolor="gray"><font color="black">...</font></td>
    <td>NA</td>
    <td>0.061</td>
    <td>0.678</td>
    <td>...</td>
  </tr>
  <tr>
    <td>230101_2_S2_L001.variants.tsv</td>
    <td>9.42</td>
    <td>99.234</td>
    <td bgcolor="gray"><font color="black">0.015</font></td>
    <td bgcolor="gray"><font color="black">0.168</font></td>
    <td bgcolor="gray"><font color="black">0.767</font></td>
    <td bgcolor="gray"><font color="black">...</font></td>
    <td>0.002</td>
    <td>0.007</td>
    <td>0.362</td>
    <td>...</td>
  </tr>
  <tr>
    <td>...</td>
    <td>...</td>
    <td>...</td>
    <td bgcolor= "gray"><font color="black">...</font></td>
    <td bgcolor= "gray"><font color="black">...</font></td>
    <td bgcolor= "gray"><font color="black">...</font></td>
    <td bgcolor= "gray"><font color="black">...</font></td>
    <td>...</td>
    <td>...</td>
    <td>...</td>
    <td>...</td>
  </tr>
</tbody>
</table>

## References

1. Karthikeyan, Smruthi, et al. "Wastewater sequencing uncovers early, cryptic SARS-CoV-2 variant transmission (preprint)." (2021).

2. Andersen-Lab Andersen-Lab/Freyja: Depth-weighted De-Mixing https://github.com/andersen-lab/Freyja (accessed Oct 4, 2023). 