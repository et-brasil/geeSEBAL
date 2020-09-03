//GEESEBAL - GOOGLE EARTH ENGINE APP FOR SURFACE ENERGY BALANCE ALGORITHM FOR LAND (SEBAL)
//CREATE BY: LEONARDO LAIPELT, RAFAEL KAYSER, ANDERSON RUHOFF AND AYAN FLEISCHMANN 
//PROJECT - ET BRASIL https://etbrasil.org/
//LAB - HIDROLOGIA DE GRANDE ESCALA [HGE] website: https://www.ufrgs.br/hge/author/hge/
//UNIVERSITY - UNIVERSIDADE FEDERAL DO RIO GRANDE DO SUL - UFRGS
//BRAZIL, RIO GRANDE DO SUL
//
//DOI
//VERSION 1.0
//CONTATC US: leonardo.laipelt@ufrgs.br

//FOLDERS //[CAUTION] If you create an own repository, change user name.
var MeteorologicalData = require('users/leolaipelt/geeSEBAL:Meteorological_data');
var Spectral_Indices = require('users/leolaipelt/geeSEBAL:spectral_indices');
var tools = require('users/leolaipelt/geeSEBAL:tools');

//IMAGE FUNCTION 
exports.ET_estimation = function(image_eeSEBAL,topNDVI,coldestTs,lowestNDVI,hottestTs,flag_meteorological){

  var image =ee.Image(image_eeSEBAL);
  var index =image.get('system:index'); //LANDSAT ID
  var zenith_angle=image.get("SOLAR_ZENITH_ANGLE"); //SOLAR ZENITH ANGLE FROM LS
  var azimuth_angle=image.get('SOLAR_AZIMUTH_ANGLE'); //SOLAR AZIMUTH ANGLE FROM LS
  var time_start=image.get('system:time_start'); //TIME START FROM LS
  var date=ee.Date(time_start); //GET EE.DATE
  var year=ee.Number(date.get('year')); //YEAR
  var month=ee.Number(date.get('month')); //MONTH
  var day=ee.Number(date.get('day')); //DAY
  var hour=ee.Number(date.get('hour')); //HOUR
  var min = ee.Number(date.get('minutes')); //MINUTES
  var crs = image.projection().crs(); //PROJECTION
  var transform = ee.List(ee.Dictionary(ee.Algorithms.Describe(image.projection())).get('transform'));
  var p_top_NDVI=ee.Number(topNDVI); //TOP NDVI PERCENTILE (FOR COLD PIXEL)
  var p_coldest_Ts=ee.Number(coldestTs); //COLDEST TS (FOR COLD PIXEL)
  var p_lowest_NDVI=ee.Number(lowestNDVI); //LOWEST NDVI (FOR HOT PIXEL)
  var p_hottest_Ts=ee.Number(hottestTs); //HOTTEST TS (FOR HOT PIXEL)
  var sun_elevation=ee.Number(90).subtract(zenith_angle); //SUN ELEVATION 
  var col_GLDAS =MeteorologicalData.get_meteorological_data(image,time_start,flag_meteorological); //METEOROLOGICAL DATA 
  var T_air = col_GLDAS.select('AirT_G'); //AIR TEMPERATURE INSTANTANEOUS [ºC]
  var ux= col_GLDAS.select('ux_G'); //WIND SPEED INSTANTANEOUS [m s-1]
  var UR = col_GLDAS.select('RH_G'); //RELATIVE HUMIDITY INSTANTANEOUS [%]
  var Rn24hobs = col_GLDAS.select('Rn24h_G'); //DAILY NET RADIATION [W m-2] (Bruin, 1987) 
  var SRTM_ELEVATION ='USGS/SRTMGL1_003'; //SRTM PRODUCT
  var srtm = ee.Image(SRTM_ELEVATION).clip(image.geometry().bounds());
  var z_alt = srtm.select('elevation');
  image=Spectral_Indices.spec_ind(image); //NDVI, NDWI, TS, SAVI, EVI, LAI
  image=tools.LST_correction(image,z_alt,T_air,UR,sun_elevation, hour, min) //CORRECT TS (USING DEM, SLOPE AND ASPECT)
  var d_cold_pixel = tools.fexp_cold_pixel(image, p_top_NDVI, p_coldest_Ts); //COLD PIXEL
  image = tools.fexp_radlong_up(image); // OUTGOING LONGWAVE RADIATION (Rl_up) [W m-2]
  image = tools.fexp_radshort_down(image, z_alt, T_air, UR,sun_elevation); // INCOMING SHORTWAVE RADIATION [W m-2]
  var n_Ts_cold = ee.Number(d_cold_pixel.get('temp'));   // COLD PIXEL - NUMBER
  image = tools.fexp_radlong_down(image, n_Ts_cold);  // INCOMING LONGWAVE RADIATION [W m-2]
  image = tools.fexp_radbalance(image); // NET RADIATON BALANCE  [W m-2] 
  image = tools.fexp_soil_heat(image); // SOIL HEAT FLUX  [W m-2] 
  var d_hot_pixel = tools.fexp_hot_pixel(image, p_lowest_NDVI, p_hottest_Ts);  // HOT PIXEL
  image = tools.fexp_sensible_heat_flux(image, ux, UR, Rn24hobs,n_Ts_cold, d_hot_pixel ); // SENSIBLE HEAT FLUX 
  image = tools.fexp_inst_et(image, Rn24hobs); //ET 24h [mm day-1] and ET inst [mm hour-1]
  
  return image;
}