# This file is derived from the rasterstats module, 
# https://github.com/perrygeo/python-rasterstats

# The following licence notice applies to rasterstats:

#Copyright (c) 2013 Matthew Perry
#All rights reserved.
#Redistribution and use in source and binary forms, with or without modification,
#are permitted provided that the following conditions are met:

#* Redistributions of source code must retain the above copyright notice, this
#  list of conditions and the following disclaimer.

#* Redistributions in binary form must reproduce the above copyright notice, this
#  list of conditions and the following disclaimer in the documentation and/or
#  other materials provided with the distribution.

#* Neither the name of the software nor the names of its
#  contributors may be used to endorse or promote products derived from
#  this software without specific prior written permission.

#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
#ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
#ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# This modification (raster_stats_multi) is (c) Harry Gibson

try:
    from osgeo import gdal,ogr
    from osgeo.gdalnumeric import *
except ImportError:
    import gdal,ogr
    from gdalnumeric import *
from rasterstats_utils import *
from shapely.geometry import shape, box, MultiPolygon
from shapely import wkt

def raster_stats_multi(vectors, rasterlist, geom_attr='GeomWKT', id_attr='fid', 
                        band_num=1, nodata_value=None, 
                        global_src_extent=False, categorical=False, stats=None, 
                        copy_properties=False, all_touched = False):
    '''
    Multi-raster version of the raster_stats (zonal_stats) function found in rasterstats package.
    
    When running zonal stats using the rasterstats package each feature (zone) must first 
    be rasterized. These are then used to mask the input raster. 
    However we often need to run raster stats on many (thousands) of input rasters 
    (all with identical geotransforms) for the same zones. 
    
    In this scenario the rasterization of the zones is a major overhead.
    This version rasterizes once and then runs the overlay against all rasters (which must have 
    the same resolution / extent as one another). It returns a generator so the stats for 
    each raster are generated when the calling code is ready for them.
    '''
    DEFAULT_STATS = ['count', 'min', 'max', 'mean']
    VALID_STATS = DEFAULT_STATS + \
        ['sum', 'std', 'median', 'majority', 'minority', 'unique', 'range']
    if not stats:
        if not categorical:
            stats = DEFAULT_STATS
        else:
            stats = []
    else:
        if isinstance(stats, basestring):
            if stats in ['*', 'ALL']:
                stats = VALID_STATS
            else:
                stats = stats.split()
    for x in stats:
        if x not in VALID_STATS:
            raise RasterStatsError("Stat `%s` not valid;" \
                " must be one of \n %r" % (x, VALID_STATS))

    run_count = False
    if categorical or 'majority' in stats or 'minority' in stats or \
       'unique' in stats:
        # run the counter once, only if needed
        run_count = True
    
    # open the first raster and use this, we will assume they are all the same size / bounds etc
    initrast = rasterlist[0]
    rds = gdal.Open(initrast, gdal.GA_ReadOnly)
    if not rds:
        raise RasterStatsError("Cannot open %r as GDAL raster" % raster)
    rb = rds.GetRasterBand(band_num)
    rgt = rds.GetGeoTransform()
    rsize = (rds.RasterXSize, rds.RasterYSize)
    rbounds = raster_extent_as_bounds(rgt, rsize)

    if nodata_value is not None:
        nodata_value = float(nodata_value)
        rb.SetNoDataValue(nodata_value)
    else:
        nodata_value = rb.GetNoDataValue()

    mem_drv = ogr.GetDriverByName('Memory')
    driver = gdal.GetDriverByName('MEM')

    results = []

    # in order to avoid re-rasterizing the zones for every values raster we've moved the rasterization out of the loop 
    # and will save the rasterized zone arrays into a dictionary (so we need enough memory to hold that)
    zoneFeatureRasters = {}
    globL = inf
    globB = inf
    globT = -inf
    globR = -inf
    
    for i,feat in enumerate(vectors):
    #for i,feat in vectors.iteritems():
        try:
            geomWKT = feat[geom_attr]
        except KeyError:
            print "No geom attr found in feature!"
            continue
        geom = wkt.loads(geomWKT)
        
        # Point and MultiPoint don't play well with GDALRasterize
        # convert them into box polygons the size of a raster cell
        buff = rgt[1] / 2.0
        if geom.type == "MultiPoint":
            geom = MultiPolygon([box(*(pt.buffer(buff).bounds)) 
                                for pt in geom.geoms])
        elif geom.type == 'Point':
            geom = box(*(geom.buffer(buff).bounds))

        ogr_geom_type = shapely_to_ogr_type(geom.type)

        # "Clip" the geometry bounds to the overall raster bounding box
        # This should avoid any rasterIO errors for partially overlapping polys
        geom_bounds = list(geom.bounds)
        if geom_bounds[0] < rbounds[0]:
            geom_bounds[0] = rbounds[0]
        if geom_bounds[1] < rbounds[1]:
            geom_bounds[1] = rbounds[1]
        if geom_bounds[2] > rbounds[2]:
            geom_bounds[2] = rbounds[2]
        if geom_bounds[3] > rbounds[3]:
            geom_bounds[3] = rbounds[3]
        
        # Record the overall bounds of the features
        if geom_bounds[0] < globL:
            globL = geom_bounds[0]
        if geom_bounds[1] < globB:
            globB = geom_bounds[1]
        if geom_bounds[2] > globR:
            globR = geom_bounds[2]
        if geom_bounds[3] > globT:
            globT = geom_bounds[3]
            
        # calculate new geotransform of the feature subset
       
        src_offset = bbox_to_pixel_offsets(rgt, geom_bounds, rsize)

        new_gt = (
            (rgt[0] + (src_offset[0] * rgt[1])),
            rgt[1],
            0.0,
            (rgt[3] + (src_offset[1] * rgt[5])),
            0.0,
            rgt[5]
        )
        fid = None
        try:
            fid= feat[id_attr]
        except KeyError:
            fid = i
        if src_offset[2] < 0 or src_offset[3] < 0:
                # we're off the raster completely, no overlap at all
                # so there's no need to even bother trying to calculate
                print "Feature "+fid+" is off raster extent - skipping!"
                zoneFeatureRasters[fid] = None
            
        else: # Create a temporary vector layer in memory
            mem_ds = mem_drv.CreateDataSource('out')
            mem_layer = mem_ds.CreateLayer('out', None, ogr_geom_type)
            ogr_feature = ogr.Feature(feature_def=mem_layer.GetLayerDefn())
            ogr_geom = ogr.CreateGeometryFromWkt(geom.wkt)
            ogr_feature.SetGeometryDirectly(ogr_geom)
            mem_layer.CreateFeature(ogr_feature)

            # Rasterize it
            rvds = driver.Create('rvds', src_offset[2], src_offset[3], 1, gdal.GDT_Byte)
            rvds.SetGeoTransform(new_gt)
            #(raster_dataset, [1], shape_layer, None, None, burn_values=[1], ['ALL_TOUCHED=TRUE']
            gdal.RasterizeLayer(rvds, [1], mem_layer, None, None, [1], ['ALL_TOUCHED='+str(all_touched)])
            rv_array = rvds.ReadAsArray()
            zoneFeatureRasters[fid] = {
                 "zonearray":rv_array,
                 "src_offset":src_offset
            }
            
    initrast=None      
    if global_src_extent:
        # outside the loop: everything except actually reading the raster data
        # create an in-memory numpy array of the source raster data
        # covering the whole extent of the vector layer
        #if strategy != "ogr":
        #    raise RasterStatsError("global_src_extent requires OGR vector")

        # find extent of ALL features
        #ds = ogr.Open(vectors)
        #layer = ds.GetLayer(layer_num)
        #ex = layer.GetExtent()
        # transform from OGR extent to xmin, ymin, xmax, ymax
        #layer_extent = (ex[0], ex[2], ex[1], ex[3])
        
        layer_extent = (globL, globB, globR, globT)
        global_src_offset = bbox_to_pixel_offsets(rgt, layer_extent, rsize)
            
    # now do the raster calculation aspects of the original task once for each input raster but getting the zone rasters from the populated dictionary
    # rather than re-rasterizing each time
    for rast in rasterlist:
        rastresults = []
        rds = gdal.Open(rast, gdal.GA_ReadOnly)
        if not rds:
           # raise RasterStatsError("Cannot open %r as GDAL raster" % rast)
           print
           print ("Cannot open %r as GDAL raster" % rast)
           print
           continue
        rb = rds.GetRasterBand(band_num)
        # we have to assume the raster size and transform are the same 
        thisRgt = rds.GetGeoTransform()
        thisRsize = (rds.RasterXSize, rds.RasterYSize)
        thisRbounds = raster_extent_as_bounds(rgt, rsize)
        if (thisRgt != rgt or thisRsize != rsize or thisRbounds != rbounds):
            print "Raster " + rast +" has differing size or geotransform from others - skipping!"
            continue

        if global_src_extent:
            global_src_array = rb.ReadAsArray(*global_src_offset)

        if nodata_value is not None:
            nodata_value = float(nodata_value)
            rb.SetNoDataValue(nodata_value)
        else:
            nodata_value = rb.GetNoDataValue()
       
        #for i, feat in enumerate(features_iter):
        # for i,feat in vectors.iteritems():
        for i, feat in enumerate(vectors):
            fid = None
            try:
                fid = feat[id_attr]
            except:
                fid = i
            if zoneFeatureRasters[fid] is None:
                # this happens when the feature was outside the raster extent so rasterizing it was skipped
                #feature_stats = dict([(s,None) for s in stats])
                continue
            else:
                zone_array = zoneFeatureRasters[fid]["zonearray"]
                src_offset = zoneFeatureRasters[fid]["src_offset"]
                if not global_src_extent:
                    # use feature's source extent and read directly from source
                    # fastest option when you have fast disks and well-indexed raster
                    # advantage: each feature uses the smallest raster chunk
                    # disadvantage: lots of disk reads on the source raster
                    src_array = rb.ReadAsArray(*src_offset)
                else:
                    # derive array from global source extent array
                    # useful *only* when disk IO or raster format inefficiencies are your limiting factor
                    # advantage: reads raster data in one pass before loop
                    # disadvantage: large vector extents combined with big rasters need lotsa memory
                    xa = src_offset[0] - global_src_offset[0]
                    ya = src_offset[1] - global_src_offset[1]
                    xb = xa + src_offset[2]
                    yb = ya + src_offset[3]
                    src_array = global_src_array[ya:yb, xa:xb]
                
                # Mask the source data array with our current feature
                # we take the logical_not to flip 0<->1 to get the correct mask effect
                # we also mask out nodata values explictly
                masked = numpy.ma.MaskedArray(
                    src_array,
                    mask=numpy.logical_or(
                        src_array == nodata_value,
                        numpy.logical_not(zone_array)
                    )
                )

                if run_count:
                    pixel_count = Counter(masked.compressed())

                if categorical:  
                    feature_stats = dict(pixel_count)
                else:
                    feature_stats = {}

                if 'min' in stats:
                    feature_stats['min'] = float(masked.min())
                if 'max' in stats:
                    feature_stats['max'] = float(masked.max())
                if 'mean' in stats:
                    feature_stats['mean'] = float(masked.mean())
                if 'count' in stats:
                    feature_stats['count'] = int(masked.count())
                # optional
                if 'sum' in stats:
                    feature_stats['sum'] = float(masked.sum())
                if 'std' in stats:
                    feature_stats['std'] = float(masked.std())
                if 'median' in stats:
                    feature_stats['median'] = float(numpy.median(masked.compressed()))
                if 'majority' in stats:
                    try:
                        feature_stats['majority'] = pixel_count.most_common(1)[0][0]
                    except IndexError:
                        feature_stats['majority'] = None
                if 'minority' in stats:
                    try:
                        feature_stats['minority'] = pixel_count.most_common()[-1][0]
                    except IndexError:
                        feature_stats['minority'] = None
                if 'unique' in stats:
                    feature_stats['unique'] = len(pixel_count.keys())
                if 'range' in stats:
                    try:
                        rmin = feature_stats['min']
                    except KeyError:
                        rmin = float(masked.min())
                    try:
                        rmax = feature_stats['max']
                    except KeyError:
                        rmax = float(masked.max())
                    feature_stats['range'] = rmax - rmin
        
            try:
                # Use the provided feature id as __fid__
                feature_stats[id_attr] = feat[id_attr]
            except:
                # use the enumerator
                feature_stats[id_attr] = i 

            if copy_properties:
                for key, val in feat.iteritems():
                    if key == id_attr or key == geom_attr:
                        continue
                    feature_stats[key] = val
            rastresults.append(feature_stats)
        yield {'rastername':rast,'stats':rastresults}
    rb = None
    rds = None
    zoneFeatureRasters = None
    ds = None
    
    