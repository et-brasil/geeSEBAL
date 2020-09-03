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

//LANDSAT 5 and 7
exports.f_cloudMaskL457_SR = function(image) {
  var quality = image.select('pixel_qa');
    var c01 = quality.eq(66); 
    var c02 = quality.eq(68); 
  var mask = c01.or(c02);
  return image.updateMask(mask);
  };
//LANDSAT 8
exports.f_cloudMaskL8_SR = function(image) {
  var quality = image.select('pixel_qa');
    var c01 = quality.eq(322); 
    var c02 = quality.eq(324); 
    var c03 = quality.eq(1346); 
  var mask = c01.or(c02).or(c03);
  return image.updateMask(mask);
  };
  //ALBEDO
//LANDSAT 5 and 7
exports.f_albedoL5L7=function(image) {
    var alfa = image.expression(
      '(0.254*B1) + (0.149*B2) + (0.147*B3) + (0.311*B4) + (0.103*B5) + (0.036*B7)',{
        'B1' : image.select(['B']).divide(10000),
        'B2' : image.select(['GR']).divide(10000),
        'B3' : image.select(['R']).divide(10000),
        'B4' : image.select(['NIR']).divide(10000),
        'B5' : image.select(['SWIR_1']).divide(10000),
        'B7' : image.select(['SWIR_2']).divide(10000)
      }).rename('ALFA');
    return image.addBands(alfa);
  };
//LANDSAT 8
exports.f_albedoL8=function(image) {
    var alfa = image.expression(
      '(0.130*B1) + (0.115*B2) + (0.143*B3) + (0.180*B4) + (0.281*B5) + (0.108*B6) + (0.042*B7)',{ 
        'B1' : image.select(['UB']).divide(10000),
        'B2' : image.select(['B']).divide(10000),
        'B3' : image.select(['GR']).divide(10000),
        'B4' : image.select(['R']).divide(10000),
        'B5' : image.select(['NIR']).divide(10000),
        'B6' : image.select(['SWIR_1']).divide(10000),
        'B7' : image.select(['SWIR_2']).divide(10000)
      }).rename('ALFA');
    return image.addBands(alfa);
  };