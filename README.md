![GitHub top language](https://img.shields.io/github/languages/top/hungtrannam/clustering-pdf-3d )
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/hungtrannam/clustering-pdf-3d )
[![GitHub contributors](https://img.shields.io/github/contributors/hungtrannam/clustering-pdf-3d )](https://github.com/hungtrannam/clustering-pdf-3d /graphs/contributors)

# Clustering 3D probability density functions 


[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=UniprJRC/FSDA&project=FSDA.prj)
[![GNU Octave](https://img.shields.io/badge/Powered_by-GNU_Octave-blue.svg)](https://www.gnu.org/software/octave/)




This package implements a clustering algorithm based on Fuzzy Clustering, k-mean, k-medoids and possibilistic c-mean, specifically designed for 3D gaussian distributions. This method is particularly effective for applications involving uncertain or overlapping data boundaries.

## Installation

To use this package, simply clone the repository and integrate the MATLAB code into your project. Ensure that you have the necessary MATLAB toolboxes (e.g., Statistics and Machine Learning Toolbox) installed to run the Gaussian probability functions. This code is compatible with MATLAB 2024b but can also be used with Octave.
```
git clone https://github.com/hungtrannam/clustering-pdf-3d.git
```
## Example Results
Data Visualization Before Clustering

<img src="https://github.com/hungtrannam/clustering-pdf-3d/blob/main/Image_Repo/FCM_dat2.jpg" alt="FCM Data" width="60%">

In this image, we visualize the initial dataset, which consists of multiple probability density functions. Each curve represents a Gaussian distribution corresponding to different data points in a 2D grid. The data points are not clustered yet, and the distributions are plotted without any clear groupings.
Clustering Results After Fuzzy Clustering

<img src="https://github.com/hungtrannam/clustering-pdf-3d/blob/main/Image_Repo/FCM_clust.jpg" alt="FCM Clustering results" width="60%">

After applying the fuzzy clustering algorithm, the data points are grouped into clusters based on their probability density functions. Each color represents a different cluster, and the boundaries are soft, meaning some points have partial membership in multiple clusters. The algorithm successfully identified the underlying cluster structure by iteratively optimizing the cluster centers and membership degrees.

<img src="https://github.com/hungtrannam/clustering-pdf-3d/blob/main/Image_Repo/FCM_heat.jpg" alt="FCM Membership" width="60%">


The clusters are visualized with both solid and dashed contours, showing the fitted cluster centers (solid lines) and the initial data points (dashed black lines). The soft boundaries between clusters allow the algorithm to account for data points that span multiple clusters.
