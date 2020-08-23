# ETBRASIL - EESEBAL

EESEBAL is a open-source implementation of Surface Energy Balance Algorithm for Land (SEBAL) using Google Earth Engine (GGE). EESEBAL is available in both Javascript and Python API.
A web application is also available (https://eesebal.etbrasil.org)

### Python functions: Image(), Collection() and TimeSeries().
### Javascript functions: image().

![EESEBAL Fluxogram](../master/Images/ilustration.png?raw=true)


## How to use Google Earth Engine?
Create an account in GEE (https://earthengine.google.com/)
 
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
from etbrasil.eesebal import Image
Image_ID=ee.Image('LANDSAT/LC08/C01/T1_SR/LC08_222081_20160118')
EESEBAL_Image=Image(Image_ID)

```
### Collection
```python
from etbrasil.eesebal import Collection

#inputs= init Year, init Month, init dat, end Year, end Month, end day, Cloud Cover
EESEBAL_Collection=Collection(2000,1,1,2010,5,6,15)
```
### TimeSeries
```python
from etbrasil.eesebal import TimeSeries

#inputs= init Year, init Month, init dat, end Year, end Month, end day, Cloud Cover,ee.Geometry.Point
point=ee.Geometry.Point([-50.161317, -9.824870])

EESEBAL_Collection=TimeSeries(2000,1,1,2010,5,6,15,point)
```

## What is SEBAL?

Surface Energy Balance Algorithm for Land (SEBAL) was developed and validated by Bastiaanssen (Bastiaanssen, 1995; Bastiaanssen et al., 1998a, 1998b) to 
estimate evapotranspiration (ET) from energy balance equation (Rn – G = LE + H), where LE, Rn, G and H are Latent Heat Flux, Net Radiation, Soil Heat Flux and Sensible Heat Flux, respectively.
SEBAL estimate LE as a residual of others energy fluxes (LE = Rn - LE - G)....

## How EESEBAL works?
![EESEBAL Fluxogram](../master/Images/Fluxogram.png?raw=true)
