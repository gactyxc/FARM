/*------------
  related Sentinel-1 preprocession function
  -------------*/
// get the Sentinel-1 image
function Sentinel_1_collection(start_data,end_data,roi){
  var s1GRD = ee.ImageCollection("COPERNICUS/S1_GRD");
  // define the filter constraints
  var criteria = ee.Filter.and(
     ee.Filter.bounds(roi), ee.Filter.date(start_data, end_data));
  //sentinel-1 data collection
  var s1GRD_imgcol = s1GRD.filter(criteria)
                   // Filter to get images collected in interferometric wide swath mode.
                   .filter(ee.Filter.eq('instrumentMode', 'IW'))
                   //.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
                   .select('VV','VH');
  var s1GRD_imgcol = s1GRD_imgcol.map(function(image){
    image = image.where(image.gt(1), 1);
    return image;
  })
  return s1GRD_imgcol;
}
// The Refined Lee speckle filter
function refinedLee(image) {
  var image_info = image.get('system:time_start');
  var bandNames = image.bandNames();
  image = ee.Image(10).pow(image.divide(10));
  
  var result = ee.ImageCollection(bandNames.map(function(b){
    var img = image.select([b]);
    
    // img must be in natural units, i.e. not in dB!
    // Set up 3x3 kernels 
    var weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
    var kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, false);
  
    var mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3);
    var variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3);
  
    // Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
    var sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]]);
  
    var sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, false);
  
    // Calculate mean and variance for the sampled windows and store as 9 bands
    var sample_mean = mean3.neighborhoodToBands(sample_kernel); 
    var sample_var = variance3.neighborhoodToBands(sample_kernel);
  
    // Determine the 4 gradients for the sampled windows
    var gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs();
    gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs());
    gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs());
    gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs());
  
    // And find the maximum gradient amongst gradient bands
    var max_gradient = gradients.reduce(ee.Reducer.max());
  
    // Create a mask for band pixels that are the maximum gradient
    var gradmask = gradients.eq(max_gradient);
  
    // duplicate gradmask bands: each gradient represents 2 directions
    gradmask = gradmask.addBands(gradmask);
  
    // Determine the 8 directions
    var directions = sample_mean.select(1).subtract(sample_mean.select(4))
                                .gt(sample_mean.select(4).subtract(sample_mean.select(7)))
                                .multiply(1);
    directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4))
                           .gt(sample_mean.select(4).subtract(sample_mean.select(2)))
                           .multiply(2));
    directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4))
                           .gt(sample_mean.select(4).subtract(sample_mean.select(5)))
                           .multiply(3));
    directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4))
                           .gt(sample_mean.select(4).subtract(sample_mean.select(8)))
                           .multiply(4));
    // The next 4 are the not() of the previous 4
    directions = directions.addBands(directions.select(0).not().multiply(5));
    directions = directions.addBands(directions.select(1).not().multiply(6));
    directions = directions.addBands(directions.select(2).not().multiply(7));
    directions = directions.addBands(directions.select(3).not().multiply(8));
  
    // Mask all values that are not 1-8
    directions = directions.updateMask(gradmask);
  
    // "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
    directions = directions.reduce(ee.Reducer.sum());  
  
    //var pal = ['ffffff','ff0000','ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff', '000000'];
    //Map.addLayer(directions.reduce(ee.Reducer.sum()), {min:1, max:8, palette: pal}, 'Directions', false);
  
    var sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));
  
    // Calculate localNoiseVariance
    var sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);
  
    // Set up the 7*7 kernels for directional statistics
    var rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4));
  
    var diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0], 
      [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]]);
  
    var rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, false);
    var diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, false);
  
    // Create stacks for mean and variance using the original kernels. Mask with relevant direction.
    var dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
    var dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1));
  
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)));
  
    // and add the bands for rotated kernels
    for (var i=1; i<4; i++) {
      dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
      dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
      dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
      dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
    }
  
    // "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
    dir_mean = dir_mean.reduce(ee.Reducer.sum());
    dir_var = dir_var.reduce(ee.Reducer.sum());
  
    // A finally generate the filtered value
    var varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0));
  
    var b = varX.divide(dir_var);
  
    return dir_mean.add(b.multiply(img.subtract(dir_mean)))
      .arrayProject([0])
      // Get a multi-band image bands.
      .arrayFlatten([['sum']])
      .float();
  })).toBands().rename(bandNames);
  var result_image = ee.Image(10).multiply(result.log10());
  return image.addBands(
    result_image.select(['VV', 'VH'], ['VV', 'VH']),
    null,
    true
  ).set('system:time_start',image_info);
}
/*------------
  related Sentinel-2 preprocession function
  -------------*/
// function to exclude bad data at scene edges
function maskEdges(s2_img) {
  return s2_img.updateMask(
      s2_img.select('B8A').mask().updateMask(s2_img.select('B9').mask()));
}
// Function to mask clouds in Sentinel-2 imagery.
function maskClouds(img) {
  var max_cloud_probabiltly = 65;
  var clouds = ee.Image(img.get('cloud_mask')).select('probability');
  var isNotCloud = clouds.lt(max_cloud_probabiltly);
  return img.updateMask(isNotCloud);
}
function Sentinel_2_collection(start_data,end_data,roi){
  var s2Sr = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED");
  var s2Clouds = ee.ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY")
  //define the filter constraints
  var criteria = ee.Filter.and(
     ee.Filter.bounds(roi), ee.Filter.date(start_data, end_data));
  //sentinel-2 data collection 
  var sentinel2_bands = ['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B11','B12'],
      new_bands = ['Aerosols','B','G','R','RE1','RE2','RE3','NIR','RE4','Water Vapor','SWIR1','SWIR2'];
  // ---Filter input collections by desired data range and region.
  s2Sr = s2Sr.filter(criteria).map(maskEdges);
  s2Clouds = s2Clouds.filter(criteria);
  // ---Join S2 SR with cloud probability dataset to add cloud mask.
  var s2SrWithCloudMask = ee.Join.saveFirst('cloud_mask').apply({
    primary: s2Sr,
    secondary: s2Clouds,
    condition:
        ee.Filter.equals({leftField: 'system:index', rightField: 'system:index'})
    });
  // ---collect the images without cloud
  var s2CloudMasked =ee.ImageCollection(s2SrWithCloudMask).map(maskClouds)
            .select(sentinel2_bands,new_bands);
  var s2SR_imgCol = s2CloudMasked//.select(sentinel2_bands,new_bands);//s2CloudMasked   s2Sr;
  return s2SR_imgCol;
}

// function to add the Legend
function addLegend(palette, names) {
 var legend = ui.Panel({
   style: {
     position: 'bottom-right',
     padding: '5px 10px'
   }
 });
 var title = ui.Label({
   value: 'Legend',
   style: {
     fontWeight: 'bold',
     color: "black",
     fontSize: '16px'
   }
 });
 legend.add(title);
 var addLegendLabel = function(color, name) {
       var showColor = ui.Label({
         style: {
           backgroundColor: '#' + color,
           padding: '8px',
           margin: '0 0 4px 0'
         }
       });
       var desc = ui.Label({
         value: name,
         style: {margin: '0 0 4px 8px'}
       });
       return ui.Panel({
         widgets: [showColor, desc],
         layout: ui.Panel.Layout.Flow('horizontal')
       });
 };
 for (var i = 0; i < palette.length; i++) {
   var label = addLegendLabel(palette[i], names[i]);
   legend.add(label);
 }  
 Map.add(legend);
}

/*------------  main producer of FARM  -------------*/

//******* step 1: get the SAR features  ------ 
var taishan = ee.FeatureCollection("users/feature_selected/Taishan");
var test_region = ee.Geometry.Polygon(
        [[[132.62810612657142, 47.50495716091717],
           [132.62810612657142, 47.38073137689873],
           [132.81212712266517, 47.38073137689873],
           [132.81212712266517, 47.50495716091717]]], null, false);
var roi = test_region;
Map.centerObject(roi,13);
var palette = ['32CD32'];
var name = ['paddy rice'];
// study year and 研究区对应的物候起止时间
var year = 2021
var grow_start = [4,20]
var grow_trans = [6,10]
var grow_end = [11,1]
var start_allYear = ee.Date.fromYMD(year,1,1)
var end_allYear =  ee.Date.fromYMD(year,12,31)
var start_grow = ee.Date.fromYMD(year,grow_start[0],grow_start[1])
var trans_grow = ee.Date.fromYMD(year,grow_trans[0],grow_trans[1])
var end_grow = ee.Date.fromYMD(year,grow_end[0],grow_end[1])
var s1GRD_imgCol = Sentinel_1_collection(start_allYear,end_allYear,roi);
var s1GRD_imgCol_rl = s1GRD_imgCol.map(refinedLee); 
var s1VH_std = s1GRD_imgCol_rl.select(['VH']).reduce(ee.Reducer.stdDev()).clip(roi).rename('VHyear_std');
var s1VH_median = s1GRD_imgCol_rl.select(['VH']).reduce(ee.Reducer.percentile([50])).clip(roi).rename('VHyear_median');
var s1VV_min = s1GRD_imgCol_rl.select(['VV']).reduce(ee.Reducer.percentile([5])).clip(roi).rename('VVyear_min');
var s1GRD_imgCol_grow = s1GRD_imgCol_rl.filterDate(start_grow,end_grow);
var s1VH_grow_min = s1GRD_imgCol_grow.select(['VH']).reduce(ee.Reducer.percentile([5])).clip(roi).rename('VHgrow_min');
/* slope of the VH */
var VH_imgCol_slope = s1GRD_imgCol_rl.map(function(image){
      var date = image.get('system:time_start');
      var time = ee.Image(ee.Date(date).difference(ee.Date('2021-01-01'),'day')).float();
      image = image.select('VH');
      return image.addBands(time.rename('DOY'));
  })
var VH_imgCol_trans = VH_imgCol_slope.filterDate(start_grow,trans_grow);
var VH_imgCol_grow  = VH_imgCol_slope.filterDate(trans_grow,end_grow);
var VHslope_trans = VH_imgCol_trans.reduce(ee.Reducer.linearFit()).select('scale').rename('VH_slopeST');
var VHslope_grow  = VH_imgCol_grow.reduce(ee.Reducer.linearFit()).select('scale').rename('VH_slopeTG');
var VH_slopeST = VHslope_trans.lte(1);
var VH_slopeTG = VHslope_grow.gte(1);
var features_slope = VH_slopeST.and(VH_slopeTG).rename('slope').toInt16();
var s1_features = s1VH_std.addBands(s1VH_median).addBands(s1VV_min).addBands(s1VH_grow_min).addBands(features_slope);
// Export the S1 features for the SNIC and samples selection
Export.image.toAsset({
  image: s1_features,
  description: 's1features',
  assetId: 's1features',
  region: roi,
  scale: 10,
  //crs: 'EPSG:32633',
  maxPixels: 1e10
});
//******* step 2: get the rice and non-rice objects based on SNIC  ------ 
// SNIC
var s1_features = ee.Image("users/feature_selected/s1features");
var SNIC_features = s1_features.select('VHyear_std','VHyear_median','VVyear_min');
var seeds = ee.Algorithms.Image.Segmentation.seedGrid(36);
var inputsSNIC = {
  image: SNIC_features, 
  size: 36,
  compactness: 5,
  connectivity: 8,
  neighborhoodSize:256,
  seeds: seeds
}
var snic = ee.Algorithms.Image.Segmentation.SNIC(inputsSNIC)
                  .reproject(SNIC_features.select(0).projection(),null,10);
var SNIC_objects = snic.select('clusters');
// get the rice and non-rice objects 
var features_VH = s1_features.select('VHyear_median','VHgrow_min');
var features_slope = s1_features.select('slope');
var object_VH = features_VH.addBands(SNIC_objects).reduceConnectedComponents(ee.Reducer.mean(),'clusters',256);
var object_slope_rice = features_slope.addBands(SNIC_objects).reduceConnectedComponents(ee.Reducer.allNonZero(),'clusters',256);
var object_slope_nonRice = features_slope.addBands(SNIC_objects).reduceConnectedComponents(ee.Reducer.anyNonZero(),'clusters',256);
// constract the rules to get the rice and non-rice objects mask
// the code for getting the T1 and T2 value was provided by stpe A1 in the final 
var T1 = -20;
var T2 = -20;
var criterion1 = object_VH.select('VHyear_median').gte(T1).and(object_VH.select('VHgrow_min').lte(T2));
var criterion2 = object_slope_rice.select('slope').eq(1);
var rice_mask = criterion1.and(criterion2).rename('rice_mask');
rice_mask = rice_mask.updateMask(rice_mask);
var criterion3 = object_slope_nonRice.select('slope').eq(0);
var non_riceMask = (criterion1.not()).and(criterion3).rename('non_rice_mask');
non_riceMask = non_riceMask.updateMask(non_riceMask);
// revert the raster to polygon for the classification 
var rice_polygon = rice_mask.reduceToVectors({
  reducer: ee.Reducer.countEvery(), 
  geometry: roi, 
  scale: 10,
  maxPixels: 1e8
});
var nonRice_polygon = non_riceMask.reduceToVectors({
  reducer: ee.Reducer.countEvery(), 
  geometry: roi, 
  scale: 10,
  maxPixels: 1e8
});
//******* step 3: classification using the multi-RF ------ 
// Generate random points within the sample area.
var rice_samples = ee.FeatureCollection.randomPoints({
  region: rice_polygon,
  points:100,
  seed: 1234}
  );
var rice_samples = rice_samples.map(function(feature){
  return feature.set('code',1)
});
var nonRice_samples = ee.FeatureCollection.randomPoints({
  region: nonRice_polygon,
  points:100,
  seed: 1234}
  );
var nonRice_samples = nonRice_samples.map(function(feature){
  return feature.set('code',0)
});
var trainSamples = rice_samples.merge(nonRice_samples);

var addVIs = function(image){
  var ndvi = image.normalizedDifference(['NIR', 'R']).rename('NDVI');
  var evi = image.expression('EVI = 2.5 * (NIR-Red) / (NIR + 6 * Red - 7.5 * Blue + 1)',{
                                        'NIR':image.select('NIR'),
                                        'Red':image.select('R'),
                                        'Blue':image.select('B')});
  var lswi = image.normalizedDifference(['NIR', 'SWIR1']).rename('LSWI');
  ndvi = ndvi.multiply(ee.Image.constant(100))
  evi = evi.multiply(ee.Image.constant(100))
  lswi = lswi.multiply(ee.Image.constant(100))
  return image.addBands(ndvi).addBands(evi).addBands(lswi);
}
// multi-RF
var year = 2021
var final_result  = ee.Image.constant(0);
var image_count = ee.Image.constant(0);
for(var i=7;i<12;i=i+1) {
  for(var j=1; j<4; j=j+1){
    //step 1: sentinel-2 data preprocession
    var start_day = ((j-1)*10+1);
    var start_data = ee.Date.fromYMD(year,i,start_day);
    var end_day = j*10;
    var end_data =  ee.Date.fromYMD(year,i,end_day);
    var s2sr_imgCol = Sentinel_2_collection(start_data,end_data,roi);
    var s2sr_imgCol = s2sr_imgCol.map(addVIs);
    var s2sr_img = s2sr_imgCol.median();
    var image_mask = s2sr_img.select('R').mask();
    var image_one  = ee.Image.constant(1).mask(image_mask);
    image_one = image_one.unmask(0);
    //step 2: classification
    var compositeFeature = s2sr_img.select('R','NIR','NDVI','EVI','LSWI');
    var training = compositeFeature.sampleRegions({
      collection: trainSamples,
      properties: ['code'],
      scale: 10
    });
    
    var trained = ee.Classifier.smileRandomForest(100).train({
      features: training,
      classProperty: 'code',
    }).setOutputMode('PROBABILITY');
    var classified = compositeFeature.classify(trained);
    var result = classified.multiply(ee.Image.constant(100)).toUint8();
    result = result.unmask(0);
    final_result = final_result.add(result);
    image_count = image_count.add(image_one);
  }
};
var final_prop = final_result.divide(image_count);
var final_rice = final_prop.gte(50).mask(final_prop.gte(50));
addLegend(palette, name);
Map.addLayer(final_rice.clip(roi),{palette:palette,min:1,max:1},'final_rice');