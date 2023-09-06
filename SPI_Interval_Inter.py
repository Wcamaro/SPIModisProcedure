import os,sys
import re
import numpy as np
from osgeo import gdal
from osgeo import gdal_array
from osgeo import osr
import datetime
from datetime import timedelta, date


def reshape_tiff(filename):
    raw_file = gdal.Open(filename)
    reshape_array = raw_file.ReadAsArray()
    
    return reshape_array



def cut(raw_matrix,bbox,cells_Size):
    """
    Function for slicing the given matrix based on the passed bounding box.
    """
    cells_degree = 1/cells_Size
    original_bbox = {
                'top': 50,
                'bottom': -50,
                'left': -180,
                'right': 180
        }
             
    cell_bbox_size = {
        'x': abs(bbox['left']-bbox['right'])*cells_degree,
        'y': abs(bbox['top']-bbox['bottom'])*cells_degree,
    }
    
    slice_start = {
        'x': int(abs(original_bbox['left']-bbox['left'])*cells_degree),
        'y': int(abs(original_bbox['top']-bbox['top'])*cells_degree),
    }
    slice_end = {
        'x': int(slice_start['x']+cell_bbox_size['x']),
        'y': int(slice_start['y']+cell_bbox_size['y']),
    }

    matrix_sliced_y = raw_matrix[slice_start['y']:slice_end['y']]
    matrix_sliced = [row[slice_start['x']:slice_end['x']] for row in matrix_sliced_y]
    
    return np.asarray(matrix_sliced)


def WriteGTiff_2(dict_x,folder,name,xmin,xmax,ymin,ymax):
    """
    Write a Gtiff, from a dictionary (Elements = Numpy Array type) defining each
    element from the dictionary like a Raster.
    """
    gdal.AllRegister()
    driver = gdal.GetDriverByName('Gtiff')
    nrows,ncols = np.shape(dict_x[dict_x.keys()[0]])
    xres = (xmax-xmin)/float(ncols)
    yres = (ymax-ymin)/float(nrows)
    geotransform = (xmin,xres,0,ymax,0,-yres)
    dict_x.keys().sort()
    for keydate in dict_x.keys():
        filename = r'%s/Trmm3B42_%s_%s.tif' % (folder,name,str(keydate))
        outDataset = driver.Create(filename,ncols,nrows,1,gdal.GDT_Float32)
        outDataset.SetGeoTransform(geotransform)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        outDataset.SetProjection(srs.ExportToWkt())
        nband = 1
        array = np.empty([nrows,ncols])
        array = dict_x[keydate]
        array[np.isnan(array)] = -99
        array[np.isinf(array)] = -99
        outband = outDataset.GetRasterBand(nband)
        outband.SetNoDataValue(-99)
        outband.WriteArray(array)
        outband.GetStatistics(0,1)
        outband = None
        outDataset = None


root_directory = r'C:\Users\ITHACA\DROUGHT\PYTHON_W\LINKS\TRMM'
out_folder = root_directory+r'\SPI_MODIS_2017\INT_SPI'
if not os.path.exists(out_folder):
                os.makedirs(out_folder)

bbox_Africa = {
                'top': 40,
                'bottom': -40,
                'left': -18,
                'right': 60
        }
date_interval = 16 ##Days to cumulate
dates = range(1,365,date_interval)
years = range(1998,2018)
filename_intervals = root_directory+r'\INTER_SPI1_SPI3\LC_Climate_Combine_SMI_MaxCumuIntervalWorld.tif'
listSPI1 = os.listdir(root_directory+r'\SPI_MODIS_2017\1')
listSPI3 = os.listdir(root_directory+r'\SPI_MODIS_2017\3')
cells_Size = 0.25
Intervals_World_array = reshape_tiff(filename_intervals)
Intervals_Africa_array = cut(Intervals_World_array,bbox_Africa,cells_Size)

for dat in dates:
    dict_SPI = {}
    for y in years:
        end_date = (date.fromordinal(date.toordinal(datetime.date(y,01,01))+dat-2))
        if end_date.year < 1998:
            print "No %s" % end_date
            continue 
        if end_date.year == 1998 and end_date.month <= 3:
            print "No %s" % end_date
            continue
        initial_date = (date.fromordinal(date.toordinal(end_date) - (30*3)+1))
        if initial_date.year < 1998:
            print "No %s" % initial_date
            continue
        if len(str(end_date.month))==2:
            mm = str(end_date.month)
        else:
            mm = '0%s' % str(end_date.month)    
        if len(str(end_date.day))==2:
            dd = str(end_date.day)
        else:
            dd = '0%s' % str(end_date.day)
        if date.toordinal(end_date) > date.toordinal(date.today()):
            print initial_date, 'NO - end_date'
            continue  
        yyyy = end_date.year
        keydate = int('%s%s%s' %(yyyy,mm,dd))
        print keydate
        matchingSPI1 = [s for s in listSPI1 if str(keydate) in s]
        matchingSPI3 = [s for s in listSPI3 if str(keydate) in s]
        SPI1_WORLD_array = reshape_tiff(root_directory+r'\SPI_MODIS_2017\1\\'+matchingSPI1[0])
        SPI3_WORLD_array = reshape_tiff(root_directory+r'\SPI_MODIS_2017\3\\'+matchingSPI3[0])
        SPI1_Africa_array = cut(SPI1_WORLD_array,bbox_Africa,cells_Size)
        SPI3_Africa_array = cut(SPI3_WORLD_array,bbox_Africa,cells_Size)
        SPI_integration = np.full(Intervals_Africa_array.shape,-99,dtype=np.float32)
        SPI_integration[Intervals_Africa_array==1] = SPI1_Africa_array[Intervals_Africa_array==1]
        SPI_integration[Intervals_Africa_array==3] = SPI3_Africa_array[Intervals_Africa_array==3]
        dict_SPI[keydate] = SPI_integration
    WriteGTiff_2(dict_SPI,out_folder,'SPI_Integration_1_3_Monthly',bbox_Africa['left'],bbox_Africa['right'],bbox_Africa['bottom'],bbox_Africa['top'])
print 'Thank You, Procedure Finished. By W.Camaro'
            
