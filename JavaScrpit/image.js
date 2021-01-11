//-------------------------------------------------------------------------------------------//
//----------------------------------------//GEESEBAL//---------------------------------------//
//GEESEBAL - GOOGLE EARTH ENGINE APP FOR SURFACE ENERGY BALANCE ALGORITHM FOR LAND (SEBAL)
//CREATE BY: LEONARDO LAIPELT, RAFAEL KAYSER, ANDERSON RUHOFF AND AYAN FLEISCHMANN 
//PROJECT - ET BRASIL https://etbrasil.org/
//LAB - HIDROLOGIA DE GRANDE ESCALA [HGE] website: https://www.ufrgs.br/hge/author/hge/
//UNIVERSITY - UNIVERSIDADE FEDERAL DO RIO GRANDE DO SUL - UFRGS
//BRAZIL, RIO GRANDE DO SUL
//
//DOI
//VERSION 0.1
//CONTATC US: leonardo.laipelt@ufrgs.br
//----------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------//


//FOLDERS //[CAUTION] If you create an own repository, change user .
var MeteorologicalData = require('users/leolaipelt/geeSEBAL:era5');
var endmembers = require('users/leolaipelt/geeSEBAL:endmembers');
var tools = require('users/leolaipelt/geeSEBAL:tools');
var evapotranspiration = require('users/leolaipelt/geeSEBAL:evapotranspiration');

//IMAGE FUNCTION 
exports.ET_estimation = function(image_geeSEBAL,topNDVI,coldestTs,lowestNDVI,hottestTs){

  var image =ee.Image(image_geeSEBAL);
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
  var sun_elevation=ee.Number(90).subtract(zenith_angle); //SUN ELEVATION 
  
  //ENDMEMBERS
  var p_top_NDVI=ee.Number(topNDVI); //TOP NDVI PERCENTILE (FOR COLD PIXEL)
  var p_coldest_Ts=ee.Number(coldestTs); //COLDEST TS (FOR COLD PIXEL)
  var p_lowest_NDVI=ee.Number(lowestNDVI); //LOWEST NDVI (FOR HOT PIXEL)
  var p_hottest_Ts=ee.Number(hottestTs); //HOTTEST TS (FOR HOT PIXEL)
  
  //METEOROLOGY PARAMETERS - GLDAS 2.1 AND 2.0
  var col_ERA5 =MeteorologicalData.era5(image,time_start);
  
  //AIR TEMPERATURE [C]
  var T_air = col_ERA5.select('AirT_G'); 
  
  //WIND SPEED [M S-1]
  var ux= col_ERA5.select('ux_G'); 
  
  //RELATIVE HUMIDITY [%]
  var UR = col_ERA5.select('RH_G'); 
  
  //NET RADIATION 24H [W M-2]
  var Rn24hobs = col_ERA5.select('Rn24h_G');  
  
  //SRTM DATA ELEVATION
  var SRTM_ELEVATION ='USGS/SRTMGL1_003'; 
  var srtm = ee.Image(SRTM_ELEVATION).clip(image.geometry().bounds());
  var z_alt = srtm.select('elevation');
  
  //SPECTRAL IMAGES (NDVI, EVI, SAVI, LAI, T_LST, e_0, e_NB, long, lat)
  image=tools.spec_ind(image); 
  
  //LAND SURFACE TEMPERATURE 
  image=tools.LST_correction(image,z_alt,T_air,UR,sun_elevation, hour, min);
  
  //COLD PIXEL
  var d_cold_pixel = endmembers.fexp_cold_pixel(image, p_top_NDVI, p_coldest_Ts); 
  
  //COLD PIXEL NUMBER
  var n_Ts_cold = ee.Number(d_cold_pixel.get('temp'));   
  
  //INSTANTANEOUS OUTGOING LONG-WAVE RADIATION [W M-2]
  image = tools.fexp_radlong_up(image); 
  
  //INSTANTANEOUS INCOMING SHORT-WAVE RADIATION [W M-2]
  image = tools.fexp_radshort_down(image, z_alt, T_air, UR,sun_elevation); 
  
  //INSTANTANEOUS INCOMING LONGWAVE RADIATION [W M-2]
  image = tools.fexp_radlong_down(image, n_Ts_cold);  
  
  //INSTANTANEOUS NET RADIATON BALANCE [W M-2]
  image = tools.fexp_radbalance(image); 
  
  //SOIL HEAT FLUX (G) [W M-2]
  image = tools.fexp_soil_heat(image); 
  
  //HOT PIXEL 
  var d_hot_pixel = endmembers.fexp_hot_pixel(image, p_lowest_NDVI, p_hottest_Ts); 
  
  //SENSIBLE HEAT FLUX (H) [W M-2]
  image = tools.fexp_sensible_heat_flux(image, ux, UR, Rn24hobs,n_Ts_cold, d_hot_pixel ); 
  
  //DAILY EVAPOTRANSPIRATION (ET_24H) [MM DAY-1]
  image = evapotranspiration.fexp_inst_et(image, Rn24hobs); 
  
  return image;
};

//=========================================REFERENCES=============================================//								
//================================================================================================//								
//ASCE–EWRI. 2005. “The ASCE standardized reference evapotranspirationequation.”								
//ASCE–EWRI Standardization of Reference Evapotranspiration Task Committe Rep.	 ASCE Reston	 Va.						
       								
//ALLEN RG	 PEREIRA LS	 RAES D	 SMITH M. 1998. Crop evapotranspiration – Guidelines					
//for computing crop water requirements. Irrigation and drainage paper FAO-56.								
//Water Resources	 Development and Management Service. Roma	 Itália. 322 p. 						
//ISBN: 0254-5293. http://www.fao.org/docrep/x0490e/x0490e00.htm.								
								
//ALLEN RG	 TASUMI M	 TREZZA R	 WATERS R	 BASTIAANSSEN W. 2002.				
//Surface Energy Balance Algorithm for Land (SEBAL) – Advanced training and user’s manual.								
//Kimberly	 USA: University of Idaho. 98 p.							
								
//ALLEN RG	 TASUMI M	 MORSE A	 TREZZA R	 WRIGHT JL	 BASTIAANSSEN WGM	 et al. 2007.		
//Satellite-Based Energy Balance for Mapping Evapotranspiration with Internalized								
//Calibration (METRIC)—Model. Journal of Irrigation and Drainage Engineering	 v. 133 (4)							
//p. 380-394. DOI: 10.1061/(ASCE)0733-9437(2007).								
								
//ALLEN R.G.	 BURNETT B.	 KRAMBER W.	 HUNTIGTON J.	 KJAERSGAARD J.	 KILIC A.	 KELLY C.	 TREZZA R.	2013
//Automated Calibration of the METRIC-Landsat Evapotranspiration Process. JAWRA J. Am. Water Resour. Assoc. 49								
//563–576. https://doi.org/10.1111/jawr.12056.								
								
//BASTIAANSSEN WGM	 MENENTI M	 FEDDES RA	 HOLTSLAG AM. 1998. A remote sensing					
//surface energy balance algorithm for land (SEBAL). 1. Formulation. Journal of Hydrology								
//v. 212-213 (1-4)	 p. 198-212. DOI: 10.1016/S0022-1694(98)00253-4.							
								
//BASTIAANSSEN WGM. 2000. SEBAL-based sensible and latent heat fluxes in the irrigated Gediz Basin								
//Turkey. Journal of Hydrology	 v. 229	 p. 87-100. DOI: 10.1016/s0022-1694(99)00202-4.						
								
//BASTIAANSSEN WGM	 1995. Regionalization of surface flux densities and							
//moisture indicators in composite terrain: a remote sensing approach under								
//clear skies in Mediterranean climates. Dr. thesis	 Wageningen Agric. Univ.							
//Wageningen Netherlands. SC-DLO	 Wageningen. https://doi.org/90-5485-465-0.							
								
//BISHT G	 VENTURINIA V	 ISLAM S	 JIANGB L. 2005. Estimation of the net radiation using MODIS					
//(Moderate Resolution Imaging Spectroradiometer) data for clear sky days. Remote Sensing of Environment								
//v. 97	 p. 52-67. DOI: 10.1016/j.rse.2005.03.014.							
								
//BRUIN	 H.A.. de	 1987. From penman to makkink	 J.C. Hooghart (Ed.)	 Proceedings and information:				
//TNO Committee on Hydrological Research. Gravenhage	 The Netherlands.							
								
//BRUTSAERT W. (1982) Evaporation into the Atmosphere: Theory	 History	 and Applications.						
//Springer	 Dordrecht	 299. http://dx.doi.org/10.1007/978-94-017-1497-6.						
 								
//CRAGO RD. 1996. Conservation and variability of the evaporative fraction during the daytime.								
//Journal of Hydrology	 v. 180 (4)	 p. 173-194. DOI: 10.1016/0022-1694(95)02903-6.						
								
//GARRISON  J.D.	 ADLER G.P.	 1990. Estimation of precipitable water over the United States						
//for application to the division of solar radiation into its direct and diffuse components. Sol. Energy 44								
//225–241. https://doi.org/10.1016/0038-092X(90)90151-2.								
								
//JAAFAR H.H.	 Ahmad	 F.A.	 2019. Time series trends of Landsat-based ET using automated					
//calibration in METRIC and SEBAL: The Bekaa Valley	 Lebanon.							
//Remote Sens. Environ. https://doi.org/https://doi.org/10.1016/j.rse.2018.12.033.								
								
//LAGOURADE JP	 BRUNET Y. 1983. A simple model for estimating the daily upward longwave surface radiation flux							
//from NOAA–AVHRR data. International Journal of Remote Sensing	 v. 14 (5)							
//p. 907-925. DOI: 10.1080/01431169308904386.								
        								
//MODIS UCSB Emissivity Library. 2004. “MODIS University of California Santa Barbara.”								
//http://www.icess.ucsb.edu/modis/EMIS/html/em.html								
								
//MARKHAM B. L.	 and BARKER J. L. 1986. “Landsat MSS and TM postcalibration dynamic ranges							
//exoatmospheric reflectances and atsatellite temperatures.”								
//EOSAT Landsat Technical Notes 1:3-8	 Earth Observation Satellite Company	 Lanham	 Md.					
								
//PAULSON CA. 1970. The mathematical representation of wind speed and temperature in								
//the unstable atmospheric surface layer. Journal of Applied Meteorology	 v. 9	 p. 857-861.						
//DOI: 10.1175/1520-0450(1970)009<c0857:atmrows>e2.0.co;b2								
								
//TASUMI M. 2003. Progress in operational estimation of regional evapotranspiration using satellite imagery.								
//PhD Thesis in Biological and Agricultural Engineering. University of Idaho. Kimberly	 USA. 379 p.							
								
//WEBB EK	 PEARMAN GI	 AND LEUNING R. 1980. Correction of flux measurements for density effects due						
//to heat and water vapor transfer. Quaterly Journal Royal Meteorological Society	 v. 106	 p.						
//85-100. DOI: 10.1002/qj.49710644707.								
