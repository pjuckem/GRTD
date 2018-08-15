"""
Functions for use with general model notebooks
"""
import gdal
import ogr

def process_raster_data(src, method, NCOL, NROW, gt, shapeproj, hnoflo, conversion=1.0):
    '''
    Takes a raster data source (ESRI grid, GeoTiff, .IMG and many other formats)
    and returns a numpy array. Arrangment of pixels is given as input and may 
    correspond to a MODFLOW grid.
    
    src : string
        complete path to raster data source
    method : string
        gdal method for interpolation. Choices are:
            gdal.GRA_NearestNeighbour 
                Nearest neighbour (select on one input pixel)
            gdal.GRA_Bilinear
                Bilinear (2x2 kernel)
            gdal.GRA_Cubic
                Cubic Convolution Approximation (4x4 kernel)
            gdal.GRA_CubicSpline
                Cubic B-Spline Approximation (4x4 kernel)
            gdal.GRA_Lanczos
                Lanczos windowed sinc interpolation (6x6 kernel)
            gdal.GRA_Average
                Average (computes the average of all non-NODATA contributing pixels)
            gdal.GRA_Mode
                Mode (selects the value which appears most often of all the sampled points)
            gdal.GRA_Max
                Max (selects maximum of all non-NODATA contributing pixels)
            gdal.GRA_Min
                Min (selects minimum of all non-NODATA contributing pixels)
            gdal.GRA_Med
                Med (selects median of all non-NODATA contributing pixels)
            gdal.GRA_Q1
                Q1 (selects first quartile of all non-NODATA contributing pixels)
            gdal.GRA_Q3
                Q3 (selects third quartile of all non-NODATA contributing pixels)

    conversion : float
        factor to be applied to raw data values to change units

    requires global variables (for now):
    NCOL, NROW : number of rows and columns
    gt : geotransform list
    shapeproj : coordinate reference system of NHDPlus (or other desired projection)
    hnoflo : to be used as missing data value (from model_spec.py)

    returns:
    2D array of raster data source projected onto model grid. 
    Returns a zero array with the correct shape if the source does not exist.
    '''
    NCOL = int(NCOL)
    NROW = int(NROW)
    hnoflo = float(hnoflo)
    if os.path.exists(src):
        rast = gdal.Open(src)

        dest = make_grid(NCOL, NROW, gt, shapeproj, hnoflo)
        gdal.ReprojectImage(rast, dest, rast.GetProjection(), shapeproj, method)

        grid = dest.GetRasterBand(1).ReadAsArray()

        grid = grid * conversion

        dest = None
        rast = None
    else:
        grid = np.ones((NROW, NCOL)) * hnoflo
        print('Data not processed for\n{}\n Check that the file exists and path is correct'.format(src))

    return grid

def process_vector_data(src, attribute, NCOL, NROW, gt, shapeproj, hnoflo):
    '''
    Takes a vector data source (ESRI shapefile) and returns a numpy array.
    Arrangment of pixels is given as input and may correspond to a MODFLOW grid.

    src : complete path to vector data source
    attribute : field in data table to assign to rasterized pixels
    
    requires global variables:
    NCOL, NROW : number of rows and columns
    gt : geotransform list
    shapeproj : coordinate reference system of NHDPlus
    hnoflo : to be used as missing data value (from model_spec.py)
    
    returns:
    2D array of vector data source projected onto model grid.
    Returns a zero array with the correct shape if the source does not exist.
    '''
    if os.path.exists(src):

        datasource = ogr.Open(src)
        layer = datasource.GetLayer()

        src = make_grid(NCOL, NROW, gt, shapeproj, 0)
        args = 'ATTRIBUTE={}'.format(attribute)
        gdal.RasterizeLayer(src, [1], layer, options = [args])

        grid = src.GetRasterBand(1).ReadAsArray()

        src = None
        dst = None        
    else:
        grid = np.ones((NROW, NCOL)) * hnoflo
        print('Data not processed for\n{}\n Check that the file exists and path is correct'.format(src))

    return grid

def make_raster(dst_file, data, NCOL, NROW, gt, proj, nodata):
    '''
    Writes numpy array to a GeoTiff file.
    
    dst_file : name of file to write
    data : 2D numpy array
    NCOL, NROW : number of rows and columns. These may coincide with a MODFLOW grid.
    gt : 6-element geotransform list [C, A, B, F, E, D]. Gives the coordinates of one pixel
        (the upper left pixel). If there is no rotation, B=D=0. If cells are square, A=-E.   
        Letter designations come from the original documentation.
        
        C = x coordinate in map units of the upper left corner of the upper left pixel
        A = distance from C along x axis to upper right pixel corner of the upper left pixel
        B = distance from C along x axis to lower left pixel corner of the upper left pixel,
        F = y coordinate in map units of the upper left corner of the upper left pixel
        E = distance from C along y axis to lower left pixel corner of the upper left pixel
        D = distance from C along y axis to upper right pixel corner of the upper left pixel
        
    proj : projection of the GeoTiff
    nodata : value to use as missing data in the GeoTiff
    '''
    import gdal
    driver = gdal.GetDriverByName("GTiff")
    dst = driver.Create(dst_file, NCOL, NROW, 1, gdal.GDT_Float32)
    dst.SetGeoTransform(gt)
    dst.SetProjection(proj)
    band = dst.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.WriteArray(data)
    dst = None

def make_grid(NCOL, NROW, gt, proj, nodata):  # NOTE: in NB1_v3 nodata defaults to hnoflow; here, it must be provided
    '''
    Creates a blank raster image in memory.
        
    NCOL, NROW : number of rows and columns. These may coincide with a MODFLOW grid.
    gt : 6-element geotransform list [C, A, B, F, E, D]. Gives the coordinates of one pixel
        (the upper left pixel). If there is no rotation, B=D=0. If cells are square, A=-E.   
        Letter designations come from the original documentation.
        
        C = x coordinate in map units of the upper left corner of the upper left pixel
        A = distance from C along x axis to upper right pixel corner of the upper left pixel
        B = distance from C along x axis to lower left pixel corner of the upper left pixel,
        F = y coordinate in map units of the upper left corner of the upper left pixel
        E = distance from C along y axis to lower left pixel corner of the upper left pixel
        D = distance from C along y axis to upper right pixel corner of the upper left pixel
        
    proj : projection of the GeoTiff
    nodata : value to use as missing data in the GeoTiff
    '''
    import gdal
    mem_drv = gdal.GetDriverByName('MEM')
    grid_ras = mem_drv.Create('', NCOL, NROW, 1, gdal.GDT_Float32)
    grid_ras.SetGeoTransform(gt)
    grid_ras.SetProjection(proj)
    band = grid_ras.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    #band.SetNoDataValue(-999.0)
    array = np.zeros((NROW,NCOL))
    band.WriteArray(array)
    return grid_ras

def process_mohp_data(tif_files):
    '''
    Loops a list of MOHP tif files. The rest of the algorithm is similar to the function
    "process_raster_data" except that a transformation from the ESRI WKT format to a generic
    format is needed. When MOHP data source is finalized, this function can be modified
    to work with the final format.
    
    src : complete path to raster data source
    method : gdal method for interpolation
    conversion : factor to be applied to raw data values to change units

    requires global variables (for now):
    NCOL, NROW : number of rows and columns
    gt : geotransform list
    shapeproj : coordinate reference system of NHDPlus (or other desired projection)
    hnoflo : to be used as missing data value (from model_spec.py)

    returns:
    2D array of raster data source projected onto model grid. Each column contains
    a different stream order MOHP. Each row corresponds to a model cell. 
    Number of rows is NCOL x NCOL. Number of columns is number of stream orders present.
    '''
    import gdal
    gdal.UseExceptions()
    import ogr
    import osr
    
    arr = np.zeros((NCOL * NROW, len(tif_files)))
    if tif_files != []:
        for col, src in enumerate(tif_files):
            hp = gdal.Open(src)

            dest = make_grid(NCOL, NROW, gt, shapeproj)

            srs = osr.SpatialReference()
            srs.ImportFromWkt(hp.GetProjection())
            srs.MorphFromESRI()
            hp_prj = srs.ExportToWkt()
            hp.SetProjection(hp_prj)

            gdal.ReprojectImage(hp, dest, hp.GetProjection(), shapeproj, gdal.GRA_NearestNeighbour)

            hp_grd = dest.GetRasterBand(1).ReadAsArray()
            hp_grd = hp_grd / 10000.

            dst = None
            hp = None
            
            arr[:, col] = hp_grd.ravel()
    return arr

def make_clockwise(coords):
    '''
    Function to determine direction of vertices of a polygon (clockwise or CCW).
    Probably not needed, but here just in case. 
    
    coords : array with dim (n, 2)
            n is number of vertices in the polygon. The last vertex is the same 
            as the first to close the polygon. The first column is x and the second is y.
    '''
    # if the points are counterclockwise, reverse them
    x1 = coords[:-1, 0]
    x2 = coords[1:, 0]
    y1 = coords[:-1, 1]
    y2 = coords[1:, 1]
    ccw = np.sum((x2 - x1) * (y2 + y1)) < 0
    if ccw:
        coords = np.flipud(coords)
        print('yup, coordinates are ccw')
        print("let's change them to CW")
    return coords

# test data for make_clockwise

# print('clockwise')
# x = np.array([1, 1, 2, 2, 1])
# y = np.array([1, 2, 2, 1, 1])
# coords = np.array(zip(x, y))
# c = make_clockwise(coords)
# print( c)
# print('\n')
# print('CCW')
# x = np.array([1, 2, 2, 1, 1])
# y = np.array([1, 1, 2, 2, 1])
# coords = np.array(zip(x, y))
# c = make_clockwise(coords)
# print( c)

import pysal as ps
import numpy as np
import pandas as pd
import os
import ast
from shutil import copyfile

def dbf2df(dbf_path, index=None, cols=False, incl_index=False):
    '''
    Read a dbf file as a pandas.DataFrame, optionally selecting the index
    variable and which columns are to be loaded.

    __author__  = "Dani Arribas-Bel <darribas@asu.edu> "
    ...

    Arguments
    ---------
    dbf_path    : str
                  Path to the DBF file to be read
    index       : str
                  Name of the column to be used as the index of the DataFrame
    cols        : list
                  List with the names of the columns to be read into the
                  DataFrame. Defaults to False, which reads the whole dbf
    incl_index  : Boolean
                  If True index is included in the DataFrame as a
                  column too. Defaults to False

    Returns
    -------
    df          : DataFrame
                  pandas.DataFrame object created
    '''
    db = ps.open(dbf_path)
    if cols:
        if incl_index:
            cols.append(index)
        vars_to_read = cols
    else:
        vars_to_read = db.header
    data = dict([(var, db.by_col(var)) for var in vars_to_read])
    if index:
        index = db.by_col(index)
        db.close()
        return pd.DataFrame(data, index=index)
    else:
        db.close()
        return pd.DataFrame(data)