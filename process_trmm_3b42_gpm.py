"""
Author: Walther C.A. Camaro G. : walther.camaro@gmail.com
"""
"""
Import process,gets data from hdf and creates matrix to calculate SPI,
calculate media and ln media from a array dictionary. 
"""

from osgeo import gdal
import numpy as np
from exceptions import TypeError, ImportError
import os,sys
import datetime
import Image
from datetime import timedelta, date
from osgeo import gdal_array
from osgeo import osr
global bbox_trmm2_daily
global cells_degree
global cells_size
import scipy as sp
from scipy import stats


"""
bbox bin
"""

bbox_trmm2_daily = {
        'top': 50,
        'bottom': -50,
        'left': -180,
        'right': 180
}


"""
General Parameters
"""

cells_size = 0.25
cells_degree = 1/cells_size

def init(filename,bbox=None,rt=None):
    """
    Loads a TRMM data into a matrix. The Open method is set to load the :0 band of the file
    which corresponds to the precipitation data. To load the error switch to the :1. By. Simone Dalmasso
    """
    ext = os.path.splitext(filename)[-1]
    bbox = checkBbox(bbox)
    reshaped = reshape(filename,ext,rt)
    rearranged = rearrange(reshaped)
    matrix = cut(rearranged,bbox)
    return matrix, bbox 
    
    

def checkBbox(bbox):
    """
    Returns the correct bbox. By Simone Dalmasso. 
    """
    if not bbox:
                    bbox = bbox_trmm2_daily
    elif not isinstance(bbox,dict):
            raise TypeError('BBOX must be a dictionary in key: values pairs !!!!')
    return bbox

def reshape(filename,ext,rt=None):
    """
    Opens and reshapes an image file. By Simone Dalmasso.
    """
    fp = open(filename, 'rb')
    data_string = fp.read()
    fp.close()
    raw_file = np.fromstring(data_string, np.float32)
    raw_file = raw_file.byteswap()
    raw_file = np.asarray(raw_file, np.float32)
    if rt == 1:
        reshaped_matrix = raw_file.reshape(480, 1440)
        reshaped_matrix = reshaped_matrix[40:440]
    else:
        reshaped_matrix = raw_file.reshape(400, 1440)
    reshaped_matrix = np.flipud(reshaped_matrix)
    return reshaped_matrix

def mask(matrix):
    """
    Applies a mask to image in order to remove the water cells.
    By Simone Dalmasso.
    """
    mask_image = '/media/sf_Share_VM/Mask_image/trmm_land_sea_wide_bin.png'
    img = Image.open(mask_image)
    imgarray = np.array(img)
    matrix = np.array(matrix)
    masked = matrix / imgarray
    return masked

def rearrange(matrix):
    b = np.split(matrix,2,axis=1)[0]
    a = np.split(matrix,2,axis=1)[1]
    rearranged = np.concatenate((a,b),axis=1)
    return rearranged
    
def cut(raw_matrix,bbox):
    """
    Fuction for slicing the given matrix based on the passed bounding box
    By Simone Dalmasso.
    """

    bbox_matrix = bbox_trmm2_daily
    cell_bbox_size = {
            'x': abs(bbox['left']-bbox['right'])*cells_degree,
            'y': abs(bbox['top']-bbox['bottom'])*cells_degree
    }
    slice_start = {
            'x': abs(bbox_matrix['left']-bbox['left'])*cells_degree,
            'y': abs(bbox_matrix['top']-bbox['top'])*cells_degree
    }
    slice_end = {
            'x': slice_start['x']+cell_bbox_size['x'],
            'y': slice_start['y']+cell_bbox_size['y'],
    }
    matrix_sliced_y = raw_matrix[slice_start['y']:slice_end['y']]
    matrix_sliced = [row[slice_start['x']:slice_end['x']] for row in matrix_sliced_y]
    return matrix_sliced

def media(dict_x):
    """
    Calculate average value between differents arrays inside at an dictionary,
    for values bigger or equal to zero.
    """
    n_elements = np.empty([len(dict_x[dict_x.keys()[0]]),len((dict_x[dict_x.keys()[0]])[0])])
    cumulate_value = np.empty([len(dict_x[dict_x.keys()[0]]),len((dict_x[dict_x.keys()[0]])[0])])

    for key in dict_x.keys():
        n_elements = n_elements +(dict_x[key]>0)*1
        n_elements[np.isnan(dict_x[key])] = np.nan
        n_elements[np.isinf(dict_x[key])] = np.inf
        cumulate_value = cumulate_value +(dict_x[key]*(dict_x[key]>0))
    mean = cumulate_value / n_elements
    return mean,n_elements

def medialn(dict_x,n_elements):
    """
    Calculate average value (lnvalue) between differents arrays inside at an dictionary,
    for values differents to zero.
    """
    cumulate_value = np.empty([len(dict_x[dict_x.keys()[0]]),len((dict_x[dict_x.keys()[0]])[0])])

    for key in dict_x.keys():
        cumulate_value = cumulate_value +(dict_x[key])
    mean = cumulate_value / n_elements
    return mean

def probnorain(dict_x):
    """
    Calculate zero discrete prob.
    """
    n_elements = np.empty([len(dict_x[dict_x.keys()[0]]),len((dict_x[dict_x.keys()[0]])[0])])
    NR = np.empty([len(dict_x[dict_x.keys()[0]]),len((dict_x[dict_x.keys()[0]])[0])])

    for key in dict_x.keys():
        n_elements = n_elements +(dict_x[key]==0)*1
        n_elements[np.isnan(dict_x[key])] = np.nan
        n_elements[np.isinf(dict_x[key])] = np.inf
        NR = NR +(dict_x[key]*(dict_x[key]>0))
    NR = n_elements / len(dict_x.keys())
    return NR

def cumulatedict(dict_x):

    cumulate_value = np.empty([len(dict_x[dict_x.keys()[0]]),len((dict_x[dict_x.keys()[0]])[0])])

    for key in dict_x.keys():
        cumulate_value = cumulate_value +(dict_x[key]*(dict_x[key]>=0))
    return cumulate_value

def calendardays(year,month):
    if (year % 4 == 0 and not year % 100 == 0)or year % 400 == 0:
        days = [31,29,31,30,31,30,30,31,30,31,30,31]
    else:
        days = [31,28,31,30,31,30,30,31,30,31,30,31]
    return days[month-1]

def WriteGTiff(dict_x,folder,month,xmin,xmax,ymin,ymax,cumulated,name):
    """
    Write a Gtiff, from a dictionary (Elements = Numpy Array type) defining each
    element from the dictionary like a band in the Raster.
    """
    gdal.AllRegister()
    driver = gdal.GetDriverByName('Gtiff')
    nyears = len(dict_x.keys())
    nrows,ncols = np.shape(dict_x[dict_x.keys()[0]])
    xres = (xmax-xmin)/float(ncols)
    yres = (ymax-ymin)/float(nrows)
    geotransform = (xmin,xres,0,ymax,0,-yres)
    filename = r'%s/Trmm3B42_%s_%s.%s.tif' % (folder,str(cumulated),month,name)
    outDataset = driver.Create(filename,ncols,nrows,nyears,gdal.GDT_Float32)
    outDataset.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    outDataset.SetProjection(srs.ExportToWkt())
    a = (map(int,dict_x.keys()))
    yearmin = min(a)
    yearmax = max(a)
    a = sorted(a)
    array = np.empty([nrows,ncols])
    nband = 1
    file = open(r'%s_Band_year_relations.csv' % filename,'w')
    file.write('Year,Gtiff_Band\n')
    for y in a:
        year = str(y)
        array = dict_x[year]
        array[np.isnan(array)] = -99
        array[np.isinf(array)] = -99
        outband = outDataset.GetRasterBand(nband)
        outband.SetNoDataValue(-99)
        outband.WriteArray(array)
        file.write('%s,%s\n' %(year,str(nband)))
        nband = nband+1
        outband.GetStatistics(0,1)
        outband = None
    outDataset = None
    file.close()
        

def WriteGTiff_2(dict_x,folder,name,xmin,xmax,ymin,ymax,cumulated):
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
    for year in dict_x.keys():
        endd = (date.fromordinal(date.toordinal(datetime.date(int(year),01,01))+cumulated-2))
        if len(str(endd.day))==1:
            end_day = '0%s' %str(endd.day)
        else:
            end_day = str(endd.day)
        if len(str(endd.month))==1:
            end_month = '0%s' %str(endd.month)
        else:
            end_month = str(endd.month)
        end_date = '%s%s%s' %(endd.year,end_month,end_day)
        filename = r'%s/Trmm3B42_%s_%s.tif' % (folder,name,end_date)
        outDataset = driver.Create(filename,ncols,nrows,1,gdal.GDT_Float32)
        outDataset.SetGeoTransform(geotransform)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        outDataset.SetProjection(srs.ExportToWkt())
        nband = 1
        array = np.empty([nrows,ncols])
        array = dict_x[year]
        array[np.isnan(array)] = -99
        array[np.isinf(array)] = -99
        outband = outDataset.GetRasterBand(nband)
        outband.SetNoDataValue(-99)
        outband.WriteArray(array)
        outband.GetStatistics(0,1)
        outband = None
        outDataset = None

def SPI_DICT(dict_matrix,P_NoRain,a,B,c0,c1,c2,d1,d2,d3):
    dict_gamma = {}
    dict_t = {}
    dict_SPI = {}
    for year in dict_matrix.keys():
        dict_gamma[year] = P_NoRain + ((1-P_NoRain)*(sp.stats.gamma.cdf(dict_matrix[year],a,loc=0,scale=B)))
        matrix_t1 = np.empty([len(dict_matrix[year]),len(dict_matrix[year][1])])
        matrix_t2 = np.empty([len(dict_matrix[year]),len(dict_matrix[year][1])])
        matrix_t = np.empty([len(dict_matrix[year]),len(dict_matrix[year][1])])
        matrix_t1=(np.log(1/((dict_gamma[year])**2)))**(0.5)
        matrix_t1[dict_gamma[year]>0.5]=0
        matrix_t1[dict_gamma[year]<0]=0
        matrix_t2=(np.log(1/((1-dict_gamma[year])**2)))**(0.5)
        matrix_t2[dict_gamma[year]<=0.5]=0
        matrix_t= matrix_t1 + matrix_t2
        matrix_t[np.isnan(dict_gamma[year])]=np.nan
        matrix_t[np.isinf(dict_gamma[year])]=np.nan
        matrix_t[B == 0]=0.0
        matrix_SPI1 = np.empty([len(dict_matrix[year]),len(dict_matrix[year][1])])
        matrix_SPI2 = np.empty([len(dict_matrix[year]),len(dict_matrix[year][1])])
        matrix_SPI = np.empty([len(dict_matrix[year]),len(dict_matrix[year][1])])
        matrix_SPI1 = -(matrix_t-((c0+c1*matrix_t+c2*(matrix_t**2))/(1+d1*matrix_t+d2*(matrix_t**2)+d3*(matrix_t**3))))
        matrix_SPI1[dict_gamma[year]>0.5]=0
        matrix_SPI1[dict_gamma[year]<=0]=0
        matrix_SPI1[np.isinf(matrix_t1)]=-30
        matrix_SPI2 = (matrix_t-((c0+c1*matrix_t+c2*(matrix_t**2))/(1+d1*matrix_t+d2*(matrix_t**2)+d3*(matrix_t**3))))
        matrix_SPI2[dict_gamma[year]<=0.5]=0
        matrix_SPI2[np.isinf(matrix_t2)]=30
        matrix_SPI= matrix_SPI1 + matrix_SPI2
        matrix_SPI[np.isnan(dict_gamma[year])]=np.nan
        matrix_SPI[np.isinf(dict_gamma[year])]=np.nan
        matrix_SPI[B == 0]=0.0
        dict_t[year]=matrix_t
        dict_SPI[year]=matrix_SPI
    return dict_SPI
