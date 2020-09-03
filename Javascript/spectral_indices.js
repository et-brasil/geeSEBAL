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
exports.spec_ind = function(image_eeSEBAL) {

    var image = image_eeSEBAL;
    var SRTM_ELEVATION ='USGS/SRTMGL1_003';
    var srtm = ee.Image(SRTM_ELEVATION).clip(image.geometry().bounds());
    var z_alt = srtm.select('elevation');
    // NDVI
    var ndvi =  image.normalizedDifference(['NIR', 'R']).rename('NDVI');
    //EVI
    var evi = image.expression('2.5 * ((N - R) / (N + (6 * R) - (7.5 * B) + 1))', {
        'N': image.select('NIR').divide(10000), 'R': image.select('R').divide(10000), 'B': image.select('B').divide(10000),}).rename('EVI');
    // SAVI
    var savi = image.expression(
      '((1 + 0.5)*(B5 - B4)) / (0.5 + (B5 + B4))', {
        'B4': image.select('R').multiply(0.0001),
        'B5': image.select('NIR').multiply(0.0001),
    }).rename('SAVI');
 
    // NDWI
    var ndwi =  image.normalizedDifference(['GR', 'NIR']).rename('NDWI');
  
    //LAI
    var savi1 = savi.where(savi.gt(0.689), 0.689);  
    var lai = image.expression(
      '-(log(( 0.69-SAVI)/0.59 ) / 0.91)', {'SAVI': savi1}).rename('LAI');
 
    //BROAD-BAND SURFACE EMISSIVITY
    var e_0 = image.expression(
      '0.95 + 0.01 * LAI',{
        'LAI': lai});
      e_0 = e_0.where(lai.gt(3), 0.98).rename('e_0');
  
    // Narrow band transmissivity
    var e_NB = image.expression(
      '0.97 + (0.0033 * LAI)',{'LAI': lai});
    e_NB = e_NB.where(lai.gt(3), 0.98).rename('e_NB');
    var log_eNB = e_NB.log();
    
  
    // LAND SURFACE TEMPERATURE 
    var comp_onda = ee.Number(1.122e-05); // ee.Number(1.155e-05) - landsat 5 ee.Number(1.089e-05) - landsat 7
    
    var lst = image.expression(
      'Tb / ( 1+ ( ( comp_onda * Tb / fator) * log_eNB))',{
        'Tb': image.select('BRT').multiply(0.1),
        'comp_onda': comp_onda,
        'log_eNB': log_eNB,
        'fator': ee.Number(1.438e-02),
      }).rename('T_LST');
      
  // GET COORDINATES
  var proj = image.select('B').projection();
  var latlon=image.select('B').addBands(ee.Image.pixelLonLat())
  //var latlon = ee.Image.pixelLonLat().reproject(proj);
  var coords = latlon.select(['longitude', 'latitude']);
 
 // FOR COLD AND HOT PIXEL SELECTION
  var pos_ndvi = ndvi.updateMask(ndvi.gt(0)).rename('pos_NDVI');  
  var ndvi_neg =  pos_ndvi.multiply(-1).rename('NDVI_neg');
 // var lst_neg = lst.multiply(-1).rename('T_LST_neg');
  var int = ee.Image(1).rename('int');
  
  //var no_water = ndwi.updateMask(ndwi.lte(0)).rename('NOWAT');  
 // var lst_nw = lst.updateMask(ndwi.lte(0)).rename('LST_NW');  
  //var lst_neg = lst.multiply(-1).rename('LST_neg');

  image = image.addBands([ndvi,evi,savi,lai,lst,e_0,e_NB,coords,ndvi_neg,pos_ndvi,int,ndwi]);
  return image;

};

