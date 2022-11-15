function configureMap(_collection, _vis_params, _label, label_position) {
  var _map = ui.Map();
  _map.addLayer(_collection, _vis_params, _label);
  _map.add(ui.Label(_label, {position:label_position}));
  return _map;
}


// Load the Sentinel-1 ImageCollection.
var sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD');
var polygon_half = ee.Geometry.Polygon( // Budrio_campo_safe_half
    [[[11.53262979564736,44.57084254751062],
    [11.53232810024896,44.57044573201654],
    [11.53264162483709,44.57033969429463],
    [11.53295082827744,44.57073804075184]]]);

// Filter by metadata properties.
var vv = sentinel1
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.eq('resolution_meters', 10))
  .filter(ee.Filter.bounds(polygon_half));
  
var vh = sentinel1
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.eq('resolution_meters', 10))
  .filter(ee.Filter.bounds(polygon_half));
  
var SP17_1 = ee.ImageCollection('COPERNICUS/S1_GRD')
                .filterBounds(polygon_half)
                .filterDate('2017-04-04', '2017-04-05')
                .sort('VV')
                //.first();

var S1 = vv.select('VV');

//Create timeseries for SP17
var S1y = S1.filterDate('2017-04-04', '2017-05-23');
var S2y = S1.filterDate('2017-09-15', '2017-11-03');
var S3y = S1.filterDate('2020-03-07', '2020-03-22');

var start_date = ee.Date('2017-04-04');
var end_date = start_date.advance(49, 'days');
var start_date2 = ee.Date('2017-09-15');
var start_date3 = ee.Date('2020-03-07');

var S1_1_mean =  S1y.filterDate(start_date, start_date.advance(1, 'days')).mean();
var S1_1_median =  S1y.filterDate(start_date, start_date.advance(1,'days')).median();
var S1_2_mean =  S1y.filterDate(start_date.advance(1, 'days'), start_date.advance(2, 'days')).mean();
var S1_2_median =  S1y.filterDate(start_date.advance(1, 'days'), start_date.advance(2, 'days')).median();
var S1_3_mean =  S1y.filterDate(start_date.advance(2, 'days'), start_date.advance(3, 'days')).mean();
var S1_3_median =  S1y.filterDate(start_date.advance(2, 'days'), start_date.advance(3, 'days')).median();
var S1_4_mean =  S1y.filterDate(start_date.advance(6, 'days'), start_date.advance(7, 'days')).mean();
var S1_4_median =  S1y.filterDate(start_date.advance(6, 'days'), start_date.advance(7, 'days')).median();
var S1_5_mean =  S1y.filterDate(start_date.advance(7, 'days'), start_date.advance(8, 'days')).mean();
var S1_5_median =  S1y.filterDate(start_date.advance(7, 'days'), start_date.advance(8, 'days')).median();
var S1_6_mean =  S1y.filterDate(start_date.advance(8, 'days'), start_date.advance(9, 'days')).mean();
var S1_6_median =  S1y.filterDate(start_date.advance(8, 'days'), start_date.advance(9, 'days')).median();
var S1_7_mean =  S1y.filterDate(start_date.advance(12, 'days'), start_date.advance(13, 'days')).mean();
var S1_7_median =  S1y.filterDate(start_date.advance(12, 'days'), start_date.advance(13, 'days')).median();
var S1_8_mean =  S1y.filterDate(start_date.advance(13, 'days'), start_date.advance(14, 'days')).mean();
var S1_8_median =  S1y.filterDate(start_date.advance(13, 'days'), start_date.advance(14, 'days')).median();
var S1_9_mean =  S1y.filterDate(start_date.advance(14, 'days'), start_date.advance(15, 'days')).mean();
var S1_9_median =  S1y.filterDate(start_date.advance(14, 'days'), start_date.advance(15, 'days')).median();
var S1_10_mean =  S1y.filterDate(start_date.advance(18, 'days'), start_date.advance(19, 'days')).mean();
var S1_10_median =  S1y.filterDate(start_date.advance(18, 'days'), start_date.advance(19, 'days')).median();
var S1_11_mean =  S1y.filterDate(start_date.advance(19, 'days'), start_date.advance(20, 'days')).mean();
var S1_11_median =  S1y.filterDate(start_date.advance(19, 'days'), start_date.advance(20, 'days')).median();
var S1_12_mean =  S1y.filterDate(start_date.advance(20, 'days'), start_date.advance(21, 'days')).mean();
var S1_12_median =  S1y.filterDate(start_date.advance(20, 'days'), start_date.advance(21, 'days')).median();
var S1_13_mean =  S1y.filterDate(start_date.advance(24, 'days'), start_date.advance(25, 'days')).mean();
var S1_13_median =  S1y.filterDate(start_date.advance(24, 'days'), start_date.advance(25, 'days')).median();
var S1_14_mean =  S1y.filterDate(start_date.advance(25, 'days'), start_date.advance(26, 'days')).mean();
var S1_14_median =  S1y.filterDate(start_date.advance(25, 'days'), start_date.advance(26, 'days')).median();
var S1_15_mean =  S1y.filterDate(start_date.advance(26, 'days'), start_date.advance(27, 'days')).mean();
var S1_15_median =  S1y.filterDate(start_date.advance(26, 'days'), start_date.advance(27, 'days')).median();
var S1_16_mean =  S1y.filterDate(start_date.advance(30, 'days'), start_date.advance(31, 'days')).mean();
var S1_16_median =  S1y.filterDate(start_date.advance(30, 'days'), start_date.advance(31, 'days')).median();
var S1_17_mean =  S1y.filterDate(start_date.advance(31, 'days'), start_date.advance(32, 'days')).mean();
var S1_17_median =  S1y.filterDate(start_date.advance(31, 'days'), start_date.advance(32, 'days')).median();
var S1_18_mean =  S1y.filterDate(start_date.advance(32, 'days'), start_date.advance(33, 'days')).mean();
var S1_18_median =  S1y.filterDate(start_date.advance(32, 'days'), start_date.advance(33, 'days')).median();
var S1_19_mean =  S1y.filterDate(start_date.advance(36, 'days'), start_date.advance(37, 'days')).mean();
var S1_19_median =  S1y.filterDate(start_date.advance(36, 'days'), start_date.advance(37, 'days')).median();
var S1_20_mean =  S1y.filterDate(start_date.advance(37, 'days'), start_date.advance(38, 'days')).mean();
var S1_20_median =  S1y.filterDate(start_date.advance(37, 'days'), start_date.advance(38, 'days')).median();
var S1_21_mean =  S1y.filterDate(start_date.advance(38, 'days'), start_date.advance(39, 'days')).mean();
var S1_21_median =  S1y.filterDate(start_date.advance(38, 'days'), start_date.advance(39, 'days')).median();
var S1_22_mean =  S1y.filterDate(start_date.advance(42, 'days'), start_date.advance(43, 'days')).mean();
var S1_22_median =  S1y.filterDate(start_date.advance(42, 'days'), start_date.advance(43, 'days')).median();
var S1_23_mean =  S1y.filterDate(start_date.advance(43, 'days'), start_date.advance(44, 'days')).mean();
var S1_23_median =  S1y.filterDate(start_date.advance(43, 'days'), start_date.advance(44, 'days')).median();
var S1_24_mean =  S1y.filterDate(start_date.advance(44, 'days'), start_date.advance(45, 'days')).mean();
var S1_24_median =  S1y.filterDate(start_date.advance(44, 'days'), start_date.advance(45, 'days')).median();
var S1_25_mean =  S1y.filterDate(start_date.advance(48, 'days'), start_date.advance(49, 'days')).mean();
var S1_25_median =  S1y.filterDate(start_date.advance(48, 'days'), start_date.advance(49, 'days')).median();
var S1_1_std =  S1y.filterDate(start_date, start_date.advance(1, 'days')).mean();
var S1_2_std =  S1y.filterDate(start_date.advance(1, 'days'), start_date.advance(2, 'days')).mean();
var S1_3_std =  S1y.filterDate(start_date.advance(2, 'days'), start_date.advance(3, 'days')).mean();
var S1_4_std =  S1y.filterDate(start_date.advance(6, 'days'), start_date.advance(7, 'days')).mean();
var S1_5_std =  S1y.filterDate(start_date.advance(7, 'days'), start_date.advance(8, 'days')).mean();
var S1_6_std =  S1y.filterDate(start_date.advance(8, 'days'), start_date.advance(9, 'days')).mean();
var S1_7_std =  S1y.filterDate(start_date.advance(12, 'days'), start_date.advance(13, 'days')).mean();
var S1_8_std =  S1y.filterDate(start_date.advance(13, 'days'), start_date.advance(14, 'days')).mean();
var S1_9_std =  S1y.filterDate(start_date.advance(14, 'days'), start_date.advance(15, 'days')).mean();
var S1_10_std =  S1y.filterDate(start_date.advance(18, 'days'), start_date.advance(19, 'days')).mean();
var S1_11_std =  S1y.filterDate(start_date.advance(19, 'days'), start_date.advance(20, 'days')).mean();
var S1_12_std =  S1y.filterDate(start_date.advance(20, 'days'), start_date.advance(21, 'days')).mean();
var S1_13_std =  S1y.filterDate(start_date.advance(24, 'days'), start_date.advance(25, 'days')).mean();
var S1_14_std =  S1y.filterDate(start_date.advance(25, 'days'), start_date.advance(26, 'days')).mean();
var S1_15_std =  S1y.filterDate(start_date.advance(26, 'days'), start_date.advance(27, 'days')).mean();
var S1_16_std =  S1y.filterDate(start_date.advance(30, 'days'), start_date.advance(31, 'days')).mean();
var S1_17_std =  S1y.filterDate(start_date.advance(31, 'days'), start_date.advance(32, 'days')).mean();
var S1_18_std =  S1y.filterDate(start_date.advance(32, 'days'), start_date.advance(33, 'days')).mean();
var S1_19_std =  S1y.filterDate(start_date.advance(36, 'days'), start_date.advance(37, 'days')).mean();
var S1_20_std =  S1y.filterDate(start_date.advance(37, 'days'), start_date.advance(38, 'days')).mean();
var S1_21_std =  S1y.filterDate(start_date.advance(38, 'days'), start_date.advance(39, 'days')).mean();
var S1_22_std =  S1y.filterDate(start_date.advance(42, 'days'), start_date.advance(43, 'days')).mean();
var S1_23_std =  S1y.filterDate(start_date.advance(43, 'days'), start_date.advance(44, 'days')).mean();
var S1_24_std =  S1y.filterDate(start_date.advance(44, 'days'), start_date.advance(45, 'days')).mean();
var S1_25_std =  S1y.filterDate(start_date.advance(48, 'days'), start_date.advance(49, 'days')).mean();

var S2_1 =  S2y.filterDate(start_date2, start_date2.advance(1, 'days')).mean();
var S2_2 =  S2y.filterDate(start_date2.advance(4, 'days'), start_date2.advance(5, 'days')).mean();
var S2_3 =  S2y.filterDate(start_date2.advance(5, 'days'), start_date2.advance(6, 'days')).mean();
var S2_4 =  S2y.filterDate(start_date2.advance(6, 'days'), start_date2.advance(7, 'days')).mean();
var S2_5 =  S2y.filterDate(start_date2.advance(10, 'days'), start_date2.advance(11, 'days')).mean();
var S2_6 =  S2y.filterDate(start_date2.advance(11, 'days'), start_date2.advance(12, 'days')).mean();
var S2_7 =  S2y.filterDate(start_date2.advance(12, 'days'), start_date2.advance(13, 'days')).mean();
var S2_8 =  S2y.filterDate(start_date2.advance(16, 'days'), start_date2.advance(17, 'days')).mean();
var S2_9 =  S2y.filterDate(start_date2.advance(17, 'days'), start_date2.advance(18, 'days')).mean();
var S2_10 =  S2y.filterDate(start_date2.advance(18, 'days'), start_date2.advance(19, 'days')).mean();
var S2_11 =  S2y.filterDate(start_date2.advance(22, 'days'), start_date2.advance(23, 'days')).mean();
var S2_12 =  S2y.filterDate(start_date2.advance(23, 'days'), start_date2.advance(24, 'days')).mean();
var S2_13 =  S2y.filterDate(start_date2.advance(24, 'days'), start_date2.advance(25, 'days')).mean();
var S2_14 =  S2y.filterDate(start_date2.advance(28, 'days'), start_date2.advance(29, 'days')).mean();
var S2_15 =  S2y.filterDate(start_date2.advance(29, 'days'), start_date2.advance(30, 'days')).mean();
var S2_16 =  S2y.filterDate(start_date2.advance(30, 'days'), start_date2.advance(31, 'days')).mean();
var S2_17 =  S2y.filterDate(start_date2.advance(34, 'days'), start_date2.advance(35, 'days')).mean();
var S2_18 =  S2y.filterDate(start_date2.advance(35, 'days'), start_date2.advance(36, 'days')).mean();
var S2_19 =  S2y.filterDate(start_date2.advance(36, 'days'), start_date2.advance(37, 'days')).mean();
var S2_20 =  S2y.filterDate(start_date2.advance(40, 'days'), start_date2.advance(41, 'days')).mean();
var S2_21 =  S2y.filterDate(start_date2.advance(41, 'days'), start_date2.advance(42, 'days')).mean();
var S2_22 =  S2y.filterDate(start_date2.advance(42, 'days'), start_date2.advance(43, 'days')).mean();
var S2_23 =  S2y.filterDate(start_date2.advance(46, 'days'), start_date2.advance(47, 'days')).mean();
var S2_24 =  S2y.filterDate(start_date2.advance(47, 'days'), start_date2.advance(48, 'days')).mean();
var S2_25 =  S2y.filterDate(start_date2.advance(48, 'days'), start_date2.advance(49, 'days')).mean();

var S3_1 =  S3y.filterDate(start_date3, start_date3.advance(1, 'days')).mean();
var S3_2 =  S3y.filterDate(start_date3.advance(1, 'days'), start_date3.advance(2, 'days')).mean();
var S3_3 =  S3y.filterDate(start_date3.advance(2, 'days'), start_date3.advance(3, 'days')).mean();
var S3_4 =  S3y.filterDate(start_date3.advance(6, 'days'), start_date3.advance(7, 'days')).mean();
var S3_5 =  S3y.filterDate(start_date3.advance(7, 'days'), start_date3.advance(8, 'days')).mean();
var S3_6 =  S3y.filterDate(start_date3.advance(8, 'days'), start_date3.advance(9, 'days')).mean();
var S3_7 =  S3y.filterDate(start_date3.advance(12, 'days'), start_date3.advance(13, 'days')).mean();
var S3_8 =  S3y.filterDate(start_date3.advance(13, 'days'), start_date3.advance(14, 'days')).mean();
var S3_9 =  S3y.filterDate(start_date3.advance(14, 'days'), start_date3.advance(15, 'days')).mean();




var compall_1 =    
ee.Image.cat([S1_1_mean,S1_2_mean,S1_3_mean,S1_4_mean,S1_5_mean,S1_6_mean,S1_7_mean,S1_8_mean,S1_9_mean,S1_10_mean,S1_11_mean,S1_12_mean,S1_13_mean,S1_14_mean,S1_15_mean,S1_16_mean,S1_17_mean,S1_18_mean,S1_19_mean,S1_20_mean,S1_21_mean,S1_22_mean,S1_23_mean,S1_24_mean,S1_25_mean]);
var compall_2 =    
ee.Image.cat([S2_1,S2_2,S2_3,S2_4,S2_5,S2_6,S2_7,S2_8,S2_9,S2_10,S2_11,S2_12,S2_13,S2_14,S2_15,S2_16,S2_17,S2_18,S2_19,S2_20,S2_21,S2_22,S2_23,S2_24,S2_25]);
var compall_3 =    
ee.Image.cat([S3_1,S3_2,S3_3,S3_4,S3_5,S3_6,S3_7,S3_8,S3_9]);

//var compall = compall.reproject({
//crs: compall.projection().crs(),
//scale: 10
//});

var vis_params1 = {bands:'VH', min:-50, max:50};
var label1 = 'ciao'

//var map1_median = configureMap(S1_1_median, vis_params1, label1, 'middle-left');
var map1_median = configureMap(S1_14_median, vis_params1, label1, 'middle-left');
map1_median.setCenter(11.532666, 44.570584, 16);

var medianDictionary = compall_3.reduceRegion({
  reducer: ee.Reducer.median(),
  geometry: polygon_half,
  scale: 10,
  maxPixels: 1e9
});

var meanDictionary = compall_3.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: polygon_half,
  scale: 10,
  maxPixels: 1e9
});

var stdDictionary = compall_3.reduceRegion({
  reducer: ee.Reducer.stdDev(),
  geometry: polygon_half,
  scale: 10,
  maxPixels: 1e9
});

print(meanDictionary);
// print(stdDictionary);
// print(SP17_1.median())

ui.root.widgets().reset([map1_median]);

// Map.addLayer(S1_1_median, vis_params1, 'ciao');
Map.addLayer((compall), {min: -25, max: 0});