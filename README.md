
# Population density trends in India (2010 vs 2015)

Gridded Population of the world is a rich datasource that combines census information from the world and uses land use/ water masks to allot them to spatial grids across the world. Let us explore this dataset for India, visualizing patterns such as spatial clusters and trends over time.


```
import georasters as gr
import rasterio
from rasterio.features import shapes
import operator
from osgeo import gdal, gdalnumeric, ogr, osr
import ogr
import subprocess
import geopandas as gp
import os
import subprocess
import glob
import pysal
import gmaps
from IPython.display import Latex
import matplotlib.pyplot as plt
import gmaps.geojson_geometries
import json
gdal.UseExceptions()

%matplotlib inline 
```


```
DATA_FOLDER = os.path.join('..', 'data')
gmaps.configure(api_key="AIzaSyBVXmDt-Mip1lZc6UOOa3S_UUQCc5ZXPR8")
```

# Get gridded population data

I [downloaded](http://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-adjusted-to-2015-unwpp-country-totals-rev10/data-download) data for the years 2010 and 2015, since the basis census data for India is from 2011, and extrapolation for 1 census period is reasonable. 


```
rstr_2015_world = os.path.join(DATA_FOLDER, '2015-world-grid', 'gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev10_2015_30_sec.tif')
rstr_2010_world = os.path.join(DATA_FOLDER, '2010-world-grid', 'gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev10_2010_30_sec.tif')
```

Since this dataset includes rasters for the whole, it must be clipped with a rasters that correspond to the country of interest. Getting the correct India boundary shape file from [this](https://github.com/datameet/maps) community sourced repository and cleaning it, to remove self-intersecting features,


```
def clean_geom(shp_file):
    gk500 = ogr.Open(
    shp_file, 1) # <====== 1 = update mode
    gk_lyr = gk500.GetLayer()

    for feature in gk_lyr:
        geom = feature.GetGeometryRef()
        if not geom.IsValid():
            feature.SetGeometry(geom.Buffer(0)) # <====== SetGeometry
            gk_lyr.SetFeature(feature) # <====== SetFeature
            assert feature.GetGeometryRef().IsValid() # Doesn't fail

    gk_lyr.ResetReading()
    assert all(feature.GetGeometryRef().IsValid() for feature in gk_lyr)
```


```
clean_geom(os.path.join(DATA_FOLDER, 'india-boundary', 'Admin2.shp'))
```

Clipping the world population grid files,


```
def clip_raster(input_file, mask_shp_file):
    output_file = input_file.replace('.tif', '_clipped.tif')
    subprocess.call('gdalwarp -cutline ' + mask_shp_file + ' -crop_to_cutline ' + input_file + ' ' + output_file, shell=True)

    # Build image overviews
    subprocess.call('gdaladdo --config COMPRESS_OVERVIEW DEFLATE ' + output_file + ' 2 4 8 16 32 64', shell=True)

    print('Done.' + '\n')
```


```
clip_raster(rstr_2010_world, os.path.join(DATA_FOLDER, 'india-boundary', 'Admin2.shp'))
clip_raster(rstr_2015_world, os.path.join(DATA_FOLDER, 'india-boundary', 'Admin2.shp'))
```

    Done.
    
    Done.
    



```
rstr_2015_world = gr.from_file(rstr_2015_world.replace('.tif', '_clipped.tif'))
rstr_2015_world = gr.from_file(rstr_2010_world.replace('.tif', '_clipped.tif'))
```

Plotting the file to check,


```
rstr_2015_world.plot()
```




    <matplotlib.image.AxesImage at 0x10a398748>




![png](output_15_1.png)


Converting the rasters to polygons for analyses,


```
def raster_to_polygons(raster):
    with rasterio.drivers():
        with rasterio.open(raster) as src:
            image = src.read(1) # first band
            results = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for i, (s, v) 
            in enumerate(
                shapes(image, mask=None, transform=src.affine)))
    return (results)
```


```
gdf_2015_world  = gp.GeoDataFrame.from_features(list(raster_to_polygons(os.path.join(DATA_FOLDER, '2015-world-grid', 'gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev10_2015_30_sec_clipped.tif'))))
gdf_2010_world  = gp.GeoDataFrame.from_features(list(raster_to_polygons(os.path.join(DATA_FOLDER, '2010-world-grid', 'gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev10_2010_30_sec_clipped.tif'))))
```


```
gdf_2015_world['ID'] = gdf_2015_world.index
gdf_2010_world['ID'] = gdf_2010_world.index
```


```
gdf_2015_world.to_file(os.path.join(DATA_FOLDER, '2015-world-grid', 'gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev10_2015_30_sec.shp'))
gdf_2010_world.to_file(os.path.join(DATA_FOLDER, '2010-world-grid', 'gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev10_2010_30_sec.shp'))
```

# Local Indicators of Spatial Association (LISA) clusters

A global measure of Spatial autocorrelation, Moran's I, describes how clustered the overall dataset is. This can be represented as, $$I = \frac N W \frac {\sum_i \sum_j w_{ij}(x_i-\bar x) (x_j-\bar x)} {\sum_i (x_i-\bar x)^2}$$

where N is the number of spatial units indexed by i and j, x  is the variable of interest, w is a matrix of spatial weights with zeroes on the diagonal (i.e., $w_{ii} = 0$) and W is the sum of $w_{ij}$. This measure however assumes homogeneity of the dataset, and does not reveal variations in this statistic across regions. Localized variations can be derived from this global measure by expressing this as a sum of the individual spatial units' I statistics,  $$I = \sum_i \frac {(x_i-\bar x)} {\sum_i (x_i-\bar x)^2} N  \sum_j w_{ij} (x_j-\bar x) $$ $$I = \sum_i \frac {(x_i-\bar x) N} {\sum_i (x_i-\bar x)^2} \sum_j w_{ij} (x_j-\bar x) $$
$$I = \sum_i \frac {I_i} {N} $$
where ${I_i} = \frac {Z_i} {m_2} \sum_j w_{ij} (x_j-\bar x)$ and $m_2 = \frac {\sum_i (x_i-\bar x)^2} {N}$

This statistic can be used to classify the observations into four clusters, those with high values surrounded by high values, low values surrounded by high values and vice versa. 

Creating queens weights and creating spatially pop density variables,


```
w_2015 = pysal.weights.Queen.from_dataframe(gdf_2015_world)
w_2014 = pysal.weights.Queen.from_dataframe(gdf_2014_world)
```

Since this will attempt to create an NXN matrix, it will take a while. I instead choose to use GeoDa to create the weights files and load them. 


```
gdf_2015_world = gp.read_file(os.path.join(DATA_FOLDER, '2015-world-grid', 'filtered.shp'))
gdf_2010_world = gp.read_file(os.path.join(DATA_FOLDER, '2010-world-grid', 'filtered.shp'))
W_2015 = pysal.open(os.path.join(DATA_FOLDER, '2015-world-grid', 'filtered_q.gal')).read()
W_2010 = pysal.open(os.path.join(DATA_FOLDER, '2010-world-grid', 'filtered_q.gal')).read()
```

    /Users/manu/anaconda3/envs/gdal2/lib/python3.6/site-packages/pysal/weights/weights.py:189: UserWarning: There are 165 disconnected observations
      warnings.warn("There are %d disconnected observations" % ni)
    /Users/manu/anaconda3/envs/gdal2/lib/python3.6/site-packages/pysal/weights/weights.py:190: UserWarning: Island ids: 73, 137, 892, 1596, 1686, 2310, 20107, 29322, 65380, 80807, 82696, 83099, 101426, 102323, 103666, 111271, 122609, 136759, 167792, 169479, 181696, 183361, 196218, 207540, 226494, 228570, 228899, 230787, 231022, 231525, 233040, 233594, 235228, 246263, 253825, 257677, 258827, 260241, 260493, 262167, 262413, 263246, 265156, 265166, 265498, 266053, 266062, 269907, 274726, 276745, 279952, 283222, 283754, 285530, 288558, 289071, 289353, 289664, 290363, 290497, 290683, 291258, 291267, 291518, 295949, 296718, 297957, 298769, 305081, 320514, 321705, 323296, 339532, 342820, 396215, 397028, 401773, 404291, 404383, 409123, 409374, 413402, 421202, 429101, 429883, 433455, 435234, 438825, 443007, 443954, 444199, 444491, 445295, 445375, 445458, 445532, 445604, 445828, 448417, 449357, 449508, 451140, 453126, 453127, 455020, 456747, 457183, 457184, 457240, 457614, 458343, 458555, 459353, 459724, 460055, 460719, 461805, 462434, 462616, 462737, 462851, 463254, 464061, 464325, 464649, 465410, 465933, 466048, 467317, 467371, 467578, 467629, 467751, 468699, 468854, 468855, 470937, 471075, 473738, 473800, 474105, 474374, 478766, 478874, 478947, 479068, 479069, 479115, 479116, 479312, 479313, 479368, 480142, 480797, 481047, 481060, 481414, 481417, 481529, 481534, 481535, 481536, 481537, 481538, 481539
      warnings.warn("Island ids: %s" % ', '.join(str(island) for island in self.islands))
    /Users/manu/anaconda3/envs/gdal2/lib/python3.6/site-packages/pysal/weights/weights.py:190: UserWarning: Island ids: 73, 137, 890, 1594, 1684, 2308, 20108, 29323, 65381, 80807, 82695, 83098, 101426, 102323, 103667, 111271, 122606, 136752, 167782, 169469, 181686, 183351, 196205, 207529, 226483, 228561, 228890, 230779, 231014, 231517, 233034, 233588, 235222, 246260, 253821, 257670, 258820, 260234, 260486, 262160, 262406, 263240, 265149, 265159, 265491, 266046, 266055, 269899, 274718, 276737, 279944, 283214, 283746, 285522, 288549, 289062, 289344, 289654, 290353, 290487, 290673, 291248, 291257, 291508, 295940, 296709, 297948, 298760, 305072, 320507, 321698, 323290, 339525, 342813, 396215, 397028, 401773, 404291, 404383, 409123, 409373, 413401, 421198, 429097, 429879, 433452, 435231, 438822, 443004, 443951, 444196, 444488, 445292, 445372, 445455, 445529, 445601, 445825, 448413, 449353, 449504, 451135, 453124, 453125, 455018, 456743, 457179, 457180, 457236, 457610, 458339, 458551, 459353, 459724, 460055, 460719, 461805, 462434, 462616, 462737, 462851, 463254, 464061, 464325, 464649, 465410, 465933, 466048, 467317, 467371, 467578, 467629, 467751, 468699, 468854, 468855, 470937, 471075, 473737, 473798, 474104, 474373, 478765, 478873, 478946, 479067, 479068, 479114, 479115, 479311, 479312, 479367, 480141, 480796, 481046, 481059, 481412, 481415, 481527, 481532, 481533, 481534, 481535, 481536, 481537
      warnings.warn("Island ids: %s" % ', '.join(str(island) for island in self.islands))



```
W_2015.transform = 'R'
W_2010.transform = 'R'
```

    ('WARNING: ', '73', ' is an island (no neighbors)')
    ('WARNING: ', '137', ' is an island (no neighbors)')
    ('WARNING: ', '892', ' is an island (no neighbors)')
    ('WARNING: ', '1596', ' is an island (no neighbors)')
    ('WARNING: ', '1686', ' is an island (no neighbors)')
    ('WARNING: ', '2310', ' is an island (no neighbors)')
    ('WARNING: ', '20107', ' is an island (no neighbors)')
    ('WARNING: ', '29322', ' is an island (no neighbors)')
    ('WARNING: ', '65380', ' is an island (no neighbors)')
    ('WARNING: ', '80807', ' is an island (no neighbors)')
    ('WARNING: ', '82696', ' is an island (no neighbors)')
    ('WARNING: ', '83099', ' is an island (no neighbors)')
    ('WARNING: ', '101426', ' is an island (no neighbors)')
    ('WARNING: ', '102323', ' is an island (no neighbors)')
    ('WARNING: ', '103666', ' is an island (no neighbors)')
    ('WARNING: ', '111271', ' is an island (no neighbors)')
    ('WARNING: ', '122609', ' is an island (no neighbors)')
    ('WARNING: ', '136759', ' is an island (no neighbors)')
    ('WARNING: ', '167792', ' is an island (no neighbors)')
    ('WARNING: ', '169479', ' is an island (no neighbors)')
    ('WARNING: ', '181696', ' is an island (no neighbors)')
    ('WARNING: ', '183361', ' is an island (no neighbors)')
    ('WARNING: ', '196218', ' is an island (no neighbors)')
    ('WARNING: ', '207540', ' is an island (no neighbors)')
    ('WARNING: ', '226494', ' is an island (no neighbors)')
    ('WARNING: ', '228570', ' is an island (no neighbors)')
    ('WARNING: ', '228899', ' is an island (no neighbors)')
    ('WARNING: ', '230787', ' is an island (no neighbors)')
    ('WARNING: ', '231022', ' is an island (no neighbors)')
    ('WARNING: ', '231525', ' is an island (no neighbors)')
    ('WARNING: ', '233040', ' is an island (no neighbors)')
    ('WARNING: ', '233594', ' is an island (no neighbors)')
    ('WARNING: ', '235228', ' is an island (no neighbors)')
    ('WARNING: ', '246263', ' is an island (no neighbors)')
    ('WARNING: ', '253825', ' is an island (no neighbors)')
    ('WARNING: ', '257677', ' is an island (no neighbors)')
    ('WARNING: ', '258827', ' is an island (no neighbors)')
    ('WARNING: ', '260241', ' is an island (no neighbors)')
    ('WARNING: ', '260493', ' is an island (no neighbors)')
    ('WARNING: ', '262167', ' is an island (no neighbors)')
    ('WARNING: ', '262413', ' is an island (no neighbors)')
    ('WARNING: ', '263246', ' is an island (no neighbors)')
    ('WARNING: ', '265156', ' is an island (no neighbors)')
    ('WARNING: ', '265166', ' is an island (no neighbors)')
    ('WARNING: ', '265498', ' is an island (no neighbors)')
    ('WARNING: ', '266053', ' is an island (no neighbors)')
    ('WARNING: ', '266062', ' is an island (no neighbors)')
    ('WARNING: ', '269907', ' is an island (no neighbors)')
    ('WARNING: ', '274726', ' is an island (no neighbors)')
    ('WARNING: ', '276745', ' is an island (no neighbors)')
    ('WARNING: ', '279952', ' is an island (no neighbors)')
    ('WARNING: ', '283222', ' is an island (no neighbors)')
    ('WARNING: ', '283754', ' is an island (no neighbors)')
    ('WARNING: ', '285530', ' is an island (no neighbors)')
    ('WARNING: ', '288558', ' is an island (no neighbors)')
    ('WARNING: ', '289071', ' is an island (no neighbors)')
    ('WARNING: ', '289353', ' is an island (no neighbors)')
    ('WARNING: ', '289664', ' is an island (no neighbors)')
    ('WARNING: ', '290363', ' is an island (no neighbors)')
    ('WARNING: ', '290497', ' is an island (no neighbors)')
    ('WARNING: ', '290683', ' is an island (no neighbors)')
    ('WARNING: ', '291258', ' is an island (no neighbors)')
    ('WARNING: ', '291267', ' is an island (no neighbors)')
    ('WARNING: ', '291518', ' is an island (no neighbors)')
    ('WARNING: ', '295949', ' is an island (no neighbors)')
    ('WARNING: ', '296718', ' is an island (no neighbors)')
    ('WARNING: ', '297957', ' is an island (no neighbors)')
    ('WARNING: ', '298769', ' is an island (no neighbors)')
    ('WARNING: ', '305081', ' is an island (no neighbors)')
    ('WARNING: ', '320514', ' is an island (no neighbors)')
    ('WARNING: ', '321705', ' is an island (no neighbors)')
    ('WARNING: ', '323296', ' is an island (no neighbors)')
    ('WARNING: ', '339532', ' is an island (no neighbors)')
    ('WARNING: ', '342820', ' is an island (no neighbors)')
    ('WARNING: ', '396215', ' is an island (no neighbors)')
    ('WARNING: ', '397028', ' is an island (no neighbors)')
    ('WARNING: ', '401773', ' is an island (no neighbors)')
    ('WARNING: ', '404291', ' is an island (no neighbors)')
    ('WARNING: ', '404383', ' is an island (no neighbors)')
    ('WARNING: ', '409123', ' is an island (no neighbors)')
    ('WARNING: ', '409374', ' is an island (no neighbors)')
    ('WARNING: ', '413402', ' is an island (no neighbors)')
    ('WARNING: ', '421202', ' is an island (no neighbors)')
    ('WARNING: ', '429101', ' is an island (no neighbors)')
    ('WARNING: ', '429883', ' is an island (no neighbors)')
    ('WARNING: ', '433455', ' is an island (no neighbors)')
    ('WARNING: ', '435234', ' is an island (no neighbors)')
    ('WARNING: ', '438825', ' is an island (no neighbors)')
    ('WARNING: ', '443007', ' is an island (no neighbors)')
    ('WARNING: ', '443954', ' is an island (no neighbors)')
    ('WARNING: ', '444199', ' is an island (no neighbors)')
    ('WARNING: ', '444491', ' is an island (no neighbors)')
    ('WARNING: ', '445295', ' is an island (no neighbors)')
    ('WARNING: ', '445375', ' is an island (no neighbors)')
    ('WARNING: ', '445458', ' is an island (no neighbors)')
    ('WARNING: ', '445532', ' is an island (no neighbors)')
    ('WARNING: ', '445604', ' is an island (no neighbors)')
    ('WARNING: ', '445828', ' is an island (no neighbors)')
    ('WARNING: ', '448417', ' is an island (no neighbors)')
    ('WARNING: ', '449357', ' is an island (no neighbors)')
    ('WARNING: ', '449508', ' is an island (no neighbors)')
    ('WARNING: ', '451140', ' is an island (no neighbors)')
    ('WARNING: ', '453126', ' is an island (no neighbors)')
    ('WARNING: ', '453127', ' is an island (no neighbors)')
    ('WARNING: ', '455020', ' is an island (no neighbors)')
    ('WARNING: ', '456747', ' is an island (no neighbors)')
    ('WARNING: ', '457183', ' is an island (no neighbors)')
    ('WARNING: ', '457184', ' is an island (no neighbors)')
    ('WARNING: ', '457240', ' is an island (no neighbors)')
    ('WARNING: ', '457614', ' is an island (no neighbors)')
    ('WARNING: ', '458343', ' is an island (no neighbors)')
    ('WARNING: ', '458555', ' is an island (no neighbors)')
    ('WARNING: ', '459353', ' is an island (no neighbors)')
    ('WARNING: ', '459724', ' is an island (no neighbors)')
    ('WARNING: ', '460055', ' is an island (no neighbors)')
    ('WARNING: ', '460719', ' is an island (no neighbors)')
    ('WARNING: ', '461805', ' is an island (no neighbors)')
    ('WARNING: ', '462434', ' is an island (no neighbors)')
    ('WARNING: ', '462616', ' is an island (no neighbors)')
    ('WARNING: ', '462737', ' is an island (no neighbors)')
    ('WARNING: ', '462851', ' is an island (no neighbors)')
    ('WARNING: ', '463254', ' is an island (no neighbors)')
    ('WARNING: ', '464061', ' is an island (no neighbors)')
    ('WARNING: ', '464325', ' is an island (no neighbors)')
    ('WARNING: ', '464649', ' is an island (no neighbors)')
    ('WARNING: ', '465410', ' is an island (no neighbors)')
    ('WARNING: ', '465933', ' is an island (no neighbors)')
    ('WARNING: ', '466048', ' is an island (no neighbors)')
    ('WARNING: ', '467317', ' is an island (no neighbors)')
    ('WARNING: ', '467371', ' is an island (no neighbors)')
    ('WARNING: ', '467578', ' is an island (no neighbors)')
    ('WARNING: ', '467629', ' is an island (no neighbors)')
    ('WARNING: ', '467751', ' is an island (no neighbors)')
    ('WARNING: ', '468699', ' is an island (no neighbors)')
    ('WARNING: ', '468854', ' is an island (no neighbors)')
    ('WARNING: ', '468855', ' is an island (no neighbors)')
    ('WARNING: ', '470937', ' is an island (no neighbors)')
    ('WARNING: ', '471075', ' is an island (no neighbors)')
    ('WARNING: ', '473738', ' is an island (no neighbors)')
    ('WARNING: ', '473800', ' is an island (no neighbors)')
    ('WARNING: ', '474105', ' is an island (no neighbors)')
    ('WARNING: ', '474374', ' is an island (no neighbors)')
    ('WARNING: ', '478766', ' is an island (no neighbors)')
    ('WARNING: ', '478874', ' is an island (no neighbors)')
    ('WARNING: ', '478947', ' is an island (no neighbors)')
    ('WARNING: ', '479068', ' is an island (no neighbors)')
    ('WARNING: ', '479069', ' is an island (no neighbors)')
    ('WARNING: ', '479115', ' is an island (no neighbors)')
    ('WARNING: ', '479116', ' is an island (no neighbors)')
    ('WARNING: ', '479312', ' is an island (no neighbors)')
    ('WARNING: ', '479313', ' is an island (no neighbors)')
    ('WARNING: ', '479368', ' is an island (no neighbors)')
    ('WARNING: ', '480142', ' is an island (no neighbors)')
    ('WARNING: ', '480797', ' is an island (no neighbors)')
    ('WARNING: ', '481047', ' is an island (no neighbors)')
    ('WARNING: ', '481060', ' is an island (no neighbors)')
    ('WARNING: ', '481414', ' is an island (no neighbors)')
    ('WARNING: ', '481417', ' is an island (no neighbors)')
    ('WARNING: ', '481529', ' is an island (no neighbors)')
    ('WARNING: ', '481534', ' is an island (no neighbors)')
    ('WARNING: ', '481535', ' is an island (no neighbors)')
    ('WARNING: ', '481536', ' is an island (no neighbors)')
    ('WARNING: ', '481537', ' is an island (no neighbors)')
    ('WARNING: ', '481538', ' is an island (no neighbors)')
    ('WARNING: ', '481539', ' is an island (no neighbors)')
    ('WARNING: ', '73', ' is an island (no neighbors)')
    ('WARNING: ', '137', ' is an island (no neighbors)')
    ('WARNING: ', '890', ' is an island (no neighbors)')
    ('WARNING: ', '1594', ' is an island (no neighbors)')
    ('WARNING: ', '1684', ' is an island (no neighbors)')
    ('WARNING: ', '2308', ' is an island (no neighbors)')
    ('WARNING: ', '20108', ' is an island (no neighbors)')
    ('WARNING: ', '29323', ' is an island (no neighbors)')
    ('WARNING: ', '65381', ' is an island (no neighbors)')
    ('WARNING: ', '80807', ' is an island (no neighbors)')
    ('WARNING: ', '82695', ' is an island (no neighbors)')
    ('WARNING: ', '83098', ' is an island (no neighbors)')
    ('WARNING: ', '101426', ' is an island (no neighbors)')
    ('WARNING: ', '102323', ' is an island (no neighbors)')
    ('WARNING: ', '103667', ' is an island (no neighbors)')
    ('WARNING: ', '111271', ' is an island (no neighbors)')
    ('WARNING: ', '122606', ' is an island (no neighbors)')
    ('WARNING: ', '136752', ' is an island (no neighbors)')
    ('WARNING: ', '167782', ' is an island (no neighbors)')
    ('WARNING: ', '169469', ' is an island (no neighbors)')
    ('WARNING: ', '181686', ' is an island (no neighbors)')
    ('WARNING: ', '183351', ' is an island (no neighbors)')
    ('WARNING: ', '196205', ' is an island (no neighbors)')
    ('WARNING: ', '207529', ' is an island (no neighbors)')
    ('WARNING: ', '226483', ' is an island (no neighbors)')
    ('WARNING: ', '228561', ' is an island (no neighbors)')
    ('WARNING: ', '228890', ' is an island (no neighbors)')
    ('WARNING: ', '230779', ' is an island (no neighbors)')
    ('WARNING: ', '231014', ' is an island (no neighbors)')
    ('WARNING: ', '231517', ' is an island (no neighbors)')
    ('WARNING: ', '233034', ' is an island (no neighbors)')
    ('WARNING: ', '233588', ' is an island (no neighbors)')
    ('WARNING: ', '235222', ' is an island (no neighbors)')
    ('WARNING: ', '246260', ' is an island (no neighbors)')
    ('WARNING: ', '253821', ' is an island (no neighbors)')
    ('WARNING: ', '257670', ' is an island (no neighbors)')
    ('WARNING: ', '258820', ' is an island (no neighbors)')
    ('WARNING: ', '260234', ' is an island (no neighbors)')
    ('WARNING: ', '260486', ' is an island (no neighbors)')
    ('WARNING: ', '262160', ' is an island (no neighbors)')
    ('WARNING: ', '262406', ' is an island (no neighbors)')
    ('WARNING: ', '263240', ' is an island (no neighbors)')
    ('WARNING: ', '265149', ' is an island (no neighbors)')
    ('WARNING: ', '265159', ' is an island (no neighbors)')
    ('WARNING: ', '265491', ' is an island (no neighbors)')
    ('WARNING: ', '266046', ' is an island (no neighbors)')
    ('WARNING: ', '266055', ' is an island (no neighbors)')
    ('WARNING: ', '269899', ' is an island (no neighbors)')
    ('WARNING: ', '274718', ' is an island (no neighbors)')
    ('WARNING: ', '276737', ' is an island (no neighbors)')
    ('WARNING: ', '279944', ' is an island (no neighbors)')
    ('WARNING: ', '283214', ' is an island (no neighbors)')
    ('WARNING: ', '283746', ' is an island (no neighbors)')
    ('WARNING: ', '285522', ' is an island (no neighbors)')
    ('WARNING: ', '288549', ' is an island (no neighbors)')
    ('WARNING: ', '289062', ' is an island (no neighbors)')
    ('WARNING: ', '289344', ' is an island (no neighbors)')
    ('WARNING: ', '289654', ' is an island (no neighbors)')
    ('WARNING: ', '290353', ' is an island (no neighbors)')
    ('WARNING: ', '290487', ' is an island (no neighbors)')
    ('WARNING: ', '290673', ' is an island (no neighbors)')
    ('WARNING: ', '291248', ' is an island (no neighbors)')
    ('WARNING: ', '291257', ' is an island (no neighbors)')
    ('WARNING: ', '291508', ' is an island (no neighbors)')
    ('WARNING: ', '295940', ' is an island (no neighbors)')
    ('WARNING: ', '296709', ' is an island (no neighbors)')
    ('WARNING: ', '297948', ' is an island (no neighbors)')
    ('WARNING: ', '298760', ' is an island (no neighbors)')
    ('WARNING: ', '305072', ' is an island (no neighbors)')
    ('WARNING: ', '320507', ' is an island (no neighbors)')
    ('WARNING: ', '321698', ' is an island (no neighbors)')
    ('WARNING: ', '323290', ' is an island (no neighbors)')
    ('WARNING: ', '339525', ' is an island (no neighbors)')
    ('WARNING: ', '342813', ' is an island (no neighbors)')
    ('WARNING: ', '396215', ' is an island (no neighbors)')
    ('WARNING: ', '397028', ' is an island (no neighbors)')
    ('WARNING: ', '401773', ' is an island (no neighbors)')
    ('WARNING: ', '404291', ' is an island (no neighbors)')
    ('WARNING: ', '404383', ' is an island (no neighbors)')
    ('WARNING: ', '409123', ' is an island (no neighbors)')
    ('WARNING: ', '409373', ' is an island (no neighbors)')
    ('WARNING: ', '413401', ' is an island (no neighbors)')
    ('WARNING: ', '421198', ' is an island (no neighbors)')
    ('WARNING: ', '429097', ' is an island (no neighbors)')
    ('WARNING: ', '429879', ' is an island (no neighbors)')
    ('WARNING: ', '433452', ' is an island (no neighbors)')
    ('WARNING: ', '435231', ' is an island (no neighbors)')
    ('WARNING: ', '438822', ' is an island (no neighbors)')
    ('WARNING: ', '443004', ' is an island (no neighbors)')
    ('WARNING: ', '443951', ' is an island (no neighbors)')
    ('WARNING: ', '444196', ' is an island (no neighbors)')
    ('WARNING: ', '444488', ' is an island (no neighbors)')
    ('WARNING: ', '445292', ' is an island (no neighbors)')
    ('WARNING: ', '445372', ' is an island (no neighbors)')
    ('WARNING: ', '445455', ' is an island (no neighbors)')
    ('WARNING: ', '445529', ' is an island (no neighbors)')
    ('WARNING: ', '445601', ' is an island (no neighbors)')
    ('WARNING: ', '445825', ' is an island (no neighbors)')
    ('WARNING: ', '448413', ' is an island (no neighbors)')
    ('WARNING: ', '449353', ' is an island (no neighbors)')
    ('WARNING: ', '449504', ' is an island (no neighbors)')
    ('WARNING: ', '451135', ' is an island (no neighbors)')
    ('WARNING: ', '453124', ' is an island (no neighbors)')
    ('WARNING: ', '453125', ' is an island (no neighbors)')
    ('WARNING: ', '455018', ' is an island (no neighbors)')
    ('WARNING: ', '456743', ' is an island (no neighbors)')
    ('WARNING: ', '457179', ' is an island (no neighbors)')
    ('WARNING: ', '457180', ' is an island (no neighbors)')
    ('WARNING: ', '457236', ' is an island (no neighbors)')
    ('WARNING: ', '457610', ' is an island (no neighbors)')
    ('WARNING: ', '458339', ' is an island (no neighbors)')
    ('WARNING: ', '458551', ' is an island (no neighbors)')
    ('WARNING: ', '459353', ' is an island (no neighbors)')
    ('WARNING: ', '459724', ' is an island (no neighbors)')
    ('WARNING: ', '460055', ' is an island (no neighbors)')
    ('WARNING: ', '460719', ' is an island (no neighbors)')
    ('WARNING: ', '461805', ' is an island (no neighbors)')
    ('WARNING: ', '462434', ' is an island (no neighbors)')
    ('WARNING: ', '462616', ' is an island (no neighbors)')
    ('WARNING: ', '462737', ' is an island (no neighbors)')
    ('WARNING: ', '462851', ' is an island (no neighbors)')
    ('WARNING: ', '463254', ' is an island (no neighbors)')
    ('WARNING: ', '464061', ' is an island (no neighbors)')
    ('WARNING: ', '464325', ' is an island (no neighbors)')
    ('WARNING: ', '464649', ' is an island (no neighbors)')
    ('WARNING: ', '465410', ' is an island (no neighbors)')
    ('WARNING: ', '465933', ' is an island (no neighbors)')
    ('WARNING: ', '466048', ' is an island (no neighbors)')
    ('WARNING: ', '467317', ' is an island (no neighbors)')
    ('WARNING: ', '467371', ' is an island (no neighbors)')
    ('WARNING: ', '467578', ' is an island (no neighbors)')
    ('WARNING: ', '467629', ' is an island (no neighbors)')
    ('WARNING: ', '467751', ' is an island (no neighbors)')
    ('WARNING: ', '468699', ' is an island (no neighbors)')
    ('WARNING: ', '468854', ' is an island (no neighbors)')
    ('WARNING: ', '468855', ' is an island (no neighbors)')
    ('WARNING: ', '470937', ' is an island (no neighbors)')
    ('WARNING: ', '471075', ' is an island (no neighbors)')
    ('WARNING: ', '473737', ' is an island (no neighbors)')
    ('WARNING: ', '473798', ' is an island (no neighbors)')
    ('WARNING: ', '474104', ' is an island (no neighbors)')
    ('WARNING: ', '474373', ' is an island (no neighbors)')
    ('WARNING: ', '478765', ' is an island (no neighbors)')
    ('WARNING: ', '478873', ' is an island (no neighbors)')
    ('WARNING: ', '478946', ' is an island (no neighbors)')
    ('WARNING: ', '479067', ' is an island (no neighbors)')
    ('WARNING: ', '479068', ' is an island (no neighbors)')
    ('WARNING: ', '479114', ' is an island (no neighbors)')
    ('WARNING: ', '479115', ' is an island (no neighbors)')
    ('WARNING: ', '479311', ' is an island (no neighbors)')
    ('WARNING: ', '479312', ' is an island (no neighbors)')
    ('WARNING: ', '479367', ' is an island (no neighbors)')
    ('WARNING: ', '480141', ' is an island (no neighbors)')
    ('WARNING: ', '480796', ' is an island (no neighbors)')
    ('WARNING: ', '481046', ' is an island (no neighbors)')
    ('WARNING: ', '481059', ' is an island (no neighbors)')
    ('WARNING: ', '481412', ' is an island (no neighbors)')
    ('WARNING: ', '481415', ' is an island (no neighbors)')
    ('WARNING: ', '481527', ' is an island (no neighbors)')
    ('WARNING: ', '481532', ' is an island (no neighbors)')
    ('WARNING: ', '481533', ' is an island (no neighbors)')
    ('WARNING: ', '481534', ' is an island (no neighbors)')
    ('WARNING: ', '481535', ' is an island (no neighbors)')
    ('WARNING: ', '481536', ' is an island (no neighbors)')
    ('WARNING: ', '481537', ' is an island (no neighbors)')



```
lisa_2015 = pysal.Moran_Local(gdf_2015_world['raster_val'].values, W_2015)
lisa_2010 = pysal.Moran_Local(gdf_2010_world['raster_val'].values, W_2010)
```

    /Users/manu/anaconda3/envs/gdal2/lib/python3.6/site-packages/pysal/esda/moran.py:851: RuntimeWarning: invalid value encountered in true_divide
      self.z_sim = (self.Is - self.EI_sim) / self.seI_sim
    /Users/manu/anaconda3/envs/gdal2/lib/python3.6/site-packages/scipy/stats/_distn_infrastructure.py:879: RuntimeWarning: invalid value encountered in greater
      return (self.a < x) & (x < self.b)
    /Users/manu/anaconda3/envs/gdal2/lib/python3.6/site-packages/scipy/stats/_distn_infrastructure.py:879: RuntimeWarning: invalid value encountered in less
      return (self.a < x) & (x < self.b)
    /Users/manu/anaconda3/envs/gdal2/lib/python3.6/site-packages/scipy/stats/_distn_infrastructure.py:1735: RuntimeWarning: invalid value encountered in greater_equal
      cond2 = (x >= self.b) & cond0



```
gdf_2015_world['C_ID'] = lisa_2015.q
gdf_2010_world['C_ID'] = lisa_2010.q
gdf_2015_world['PP_VAL'] = lisa_2015.p_sim 
gdf_2010_world['PP_VAL'] = lisa_2010.p_sim
```


```
gdf_2015_world['significance'] = gdf_2015_world['PP_VAL'] < 0.5
gdf_2010_world['significance'] = gdf_2010_world['PP_VAL'] < 0.5
```

# 2010 LISA Clusters


```
gdf_india = gp.read_file(os.path.join(DATA_FOLDER, 'india-boundary', 'Admin2.shp'))
```


```
base = gdf_india.plot(figsize=(50,50), color='grey')
gdf_2010_world.loc[(gdf_2010_world['C_ID']==0) & (gdf_2010_world['significance']==True), ].plot(ax=base, color='red')
gdf_2010_world.loc[(gdf_2010_world['C_ID']==2) & (gdf_2010_world['significance']==True), ].plot(ax=base, color='blue')
# gdf_2010_world.loc[(gdf_2010_world['C_ID']==1) & (gdf_2010_world['significance']==True), ].plot(ax=base, color='#83cef4')
# gdf_2010_world.loc[(gdf_2010_world['C_ID']==3) & (gdf_2010_world['significance']==True), ].plot(ax=base, color='#e59696')
```




    <matplotlib.axes._subplots.AxesSubplot at 0x18c6dc748>




![png](output_35_1.png)


# 2015 LISA Clusters


```
base = gdf_india.plot(figsize=(50,50), color='grey')
gdf_2015_world.loc[(gdf_2015_world['C_ID']==0) & (gdf_2015_world['significance']==True), ].plot(ax=base, color='red')
gdf_2015_world.loc[(gdf_2015_world['C_ID']==2) & (gdf_2015_world['significance']==True), ].plot(ax=base, color='blue')
# gdf_2015_world.loc[(gdf_2015_world['C_ID']==1) & (gdf_2015_world['significance']==True), ].plot(ax=base, color='#83cef4')
# gdf_2015_world.loc[(gdf_2015_world['C_ID']==3) & (gdf_2015_world['significance']==True), ].plot(ax=base, color='#e59696')
```




    <matplotlib.axes._subplots.AxesSubplot at 0x195812940>




![png](output_37_1.png)


While the overall trend of high density clusters remain centered around big cities, notable observations include,

1. New high density clusters along West Bengal/ Bangladesh border in 2015.
2. High density clusters around big cities are slightly bigger in 2015, implying more suburban developments. 

Using google maps as base,


```
fig = gmaps.figure()
hh_layer = gmaps.geojson_layer(
              json.loads(gdf_2010_world.loc[(gdf_2010_world['C_ID']==0) & (gdf_2010_world['significance']==True), ].geometry.to_json()),
              fill_color='red',
              stroke_color='red',
              fill_opacity=0.8)
ll_layer = gmaps.geojson_layer(
              json.loads(gdf_2010_world.loc[(gdf_2010_world['C_ID']==2) & (gdf_2010_world['significance']==True), ].geometry.to_json()),
              fill_color='blue',
              stroke_color='blue',
              fill_opacity=0.8)
lh_layer = gmaps.geojson_layer(
              json.loads(gdf_2010_world.loc[(gdf_2010_world['C_ID']==1) & (gdf_2010_world['significance']==True), ].geometry.to_json()),
              fill_color='#83cef4',
              stroke_color='#83cef4',
              fill_opacity=0.8)
hl_layer = gmaps.geojson_layer(
              json.loads(gdf_2010_world.loc[(gdf_2010_world['C_ID']==3) & (gdf_2010_world['significance']==True), ].geometry.to_json()),
              fill_color='#e59696',
              stroke_color='#e59696',
              fill_opacity=0.8)
fig.add_layer(hh_layer)
fig.add_layer(ll_layer)
fig.add_layer(lh_layer)
fig.add_layer(hl_layer)
fig
```

# 2015 LISA Clusters


```
fig = gmaps.figure()
hh_layer = gmaps.geojson_layer(
              json.loads(gdf_2015_world.loc[(gdf_2015_world['C_ID']==0) & (gdf_2015_world['significance']==True), ].geometry.to_json()),
              fill_color='red',
              stroke_color='red',
              fill_opacity=0.8)
ll_layer = gmaps.geojson_layer(
              json.loads(gdf_2015_world.loc[(gdf_2015_world['C_ID']==2) & (gdf_2015_world['significance']==True), ].geometry.to_json()),
              fill_color='blue',
              stroke_color='blue',
              fill_opacity=0.8)
lh_layer = gmaps.geojson_layer(
              json.loads(gdf_2015_world.loc[(gdf_2015_world['C_ID']==1) & (gdf_2015_world['significance']==True), ].geometry.to_json()),
              fill_color='#83cef4',
              stroke_color='#83cef4',
              fill_opacity=0.8)
hl_layer = gmaps.geojson_layer(
              json.loads(gdf_2015_world.loc[(gdf_2015_world['C_ID']==3) & (gdf_2015_world['significance']==True), ].geometry.to_json()),
              fill_color='#e59696',
              stroke_color='#e59696',
              fill_opacity=0.8)
fig.add_layer(hh_layer)
fig.add_layer(ll_layer)
fig.add_layer(lh_layer)
fig.add_layer(hl_layer)
fig
```
