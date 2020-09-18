# ETBRASIL - geeSEBAL
<img src="https://github.com/et-brasil/EESEBAL/blob/master/Images/geeSEBAL_logo_update_cut.png?raw=true" width="200">


geeSEBAL is a open-source implementation of Surface Energy Balance Algorithm for Land (SEBAL) using Google Earth Engine (GEE). geeSEBAL is available in both Javascript and Python API.\
\
A web application is also available (https://etbrasil.org/geesebal).

### Python functions: Image(), Collection() and TimeSeries().
### Javascript functions: image().

![geeSEBAL Fluxogram](../master/Images/ilustration.png?raw=true )


## How to use Google Earth Engine?

You need an account in GEE (https://earthengine.google.com/).
 
### Python API

Using pip to install earthengine-api

```bash
pip install earthengine-api
```
Authenticate Earth Engine library
```bash
import ee; ee.Authenticate()
```
### JavaScript API

Editor code available (https://code.earthengine.google.com/)

## Examples
### Image
```python
from etbrasil.geesebal import Image
Image_ID=ee.Image('LANDSAT/LC08/C01/T1_SR/LC08_222081_20160118')
geeSEBAL_Image=Image(Image_ID)

```
### Collection
```python
from etbrasil.geesebal import Collection

#inputs= init Year, init Month, init dat, end Year, end Month, end day, Cloud Cover
geeSEBAL_Collection=Collection(2000,1,1,2010,5,6,15)
```
### TimeSeries
```python
from etbrasil.eesebal import TimeSeries

#inputs= init Year, init Month, init dat, end Year, end Month, end day, Cloud Cover,ee.Geometry.Point
point=ee.Geometry.Point([-50.161317, -9.824870])

geeSEBAL_Collection=TimeSeries(2000,1,1,2010,5,6,15,point)
```

## What is SEBAL?

Surface Energy Balance Algorithm for Land (SEBAL) was developed and validated by Bastiaanssen (Bastiaanssen, 1995; Bastiaanssen et al., 1998a, 1998b) to 
estimate evapotranspiration (ET) from energy balance equation (Rn – G = LE + H), where LE, Rn, G and H are Latent Heat Flux, Net Radiation, Soil Heat Flux and Sensible Heat Flux, respectively.

Working....

## How geeSEBAL works?
![geeSEBAL Fluxogram](../master/Images/Fluxogram.png?raw=true)

## More about SEBAL:
#### [Bastiaanssen, W.G.M., 1995. Regionalization of surface flux densities and moisture indicators in composite terrain: a remote sensing approach under clear skies in Mediterranean climates. Dr. thesis, Wageningen Agric. Univ. Wageningen Netherlands. SC-DLO, Wageningen. ](https://doi.org/90-5485-465-0)
#### [Bastiaanssen, W.G.M., Menenti, M., Feddes, R.A., Holtslag, A.A.M., 1998a. A remote sensing surface energy balance algorithm for land (SEBAL): 1. Formulation. J. Hydrol. 212–213, 198–212.](https://doi.org/10.1016/S0022-1694(98)00253-4)
#### [Bastiaanssen, W.G.M., Pelgrum, H., Wang, J., Ma, Y., Moreno, J.F., Roerink, G.J., van der Wal, T., 1998c. A remote sensing surface energy balance algorithm for land (SEBAL): 2. Validation. J. Hydrol. 212–213, 213–229.](https://doi.org/10.1016/S0022-1694(98)00254-6)

## Our work:

##### [Ruhoff, A., Paz, A., Collischonn, W., Aragão, L., Rocha, H., S. Malhi, Y., 2012. A MODIS-Based Energy Balance to Estimate Evapotranspiration for Clear-Sky Days in Brazilian Tropical Savannas, Remote Sensing, vol. 4, issue 3, pp. 703-725.](https://doi.org/10.3390/rs4030703)

##### [Ruhoff, A. R. Paz, L. E. O. C. Aragao, Q. Mu, Y. Malhi, W. Collischonn, H. R. Rocha & S. W. Running (2013) Assessment of the MODIS global evapotranspiration algorithm using eddy covariance measurements and hydrological modelling in the Rio Grande basin, Hydrological Sciences Journal, 58:8, 1658-1676](https://DOI:10.1080/02626667.2013.837578)
##### [Laipelt, L.; Ruhoff, A.L.; Fleischmann, A.S.; Kayser, R.H.B.; Kich, E.M.; da Rocha, H.R.; Neale, C.M.U. Assessment of an Automated Calibration of the SEBAL Algorithm to Estimate Dry-Season Surface-Energy Partitioning in a Forest–Savanna Transition in Brazil. Remote Sens. 2020, 12, 1108.](https://doi.org/10.3390/rs12071108)



