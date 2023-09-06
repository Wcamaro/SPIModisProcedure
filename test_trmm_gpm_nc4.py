import numpy as np
from netCDF4 import Dataset
from osgeo import gdal
from osgeo import gdal_array
from osgeo import osr
import process_trmm_3b42_gpm 
from process_trmm_3b42_gpm import init


def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)



def WriteGTiff_2(array_x,folder,name,xmin,xmax,ymin,ymax,fill_value):
    """
    Write a Gtiff, from a dictionary (Elements = Numpy Array type) defining each
    element from the dictionary like a Raster.
    """
    gdal.AllRegister()
    driver = gdal.GetDriverByName('Gtiff')
    nrows,ncols = np.shape(array_x)
    xres = (xmax-xmin)/float(ncols)
    yres = (ymax-ymin)/float(nrows)
    geotransform = (xmin,xres,0,ymax,0,-yres)
    filename = r'%s/%s.tif' % (folder,name)
    outDataset = driver.Create(filename,ncols,nrows,1,gdal.GDT_Float32)
    outDataset.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    outDataset.SetProjection(srs.ExportToWkt())
    nband = 1
    outband = outDataset.GetRasterBand(nband)
    outband.SetNoDataValue(fill_value)
    outband.WriteArray(array_x)
    outband.GetStatistics(0,1)
    outband = None
    outDataset = None




trmm = Dataset(r'/media/sf_Share_VM/Nc4_test/3B42RT_Daily.20161031.7.nc4',"r")
gpm = Dataset(r'/media/sf_Share_VM/Nc4_test/3B-DAY-L.MS.MRG.3IMERG.20161031-S000000-E235959.V03.nc4',"r")
trmm_bin_filename = (r'/media/sf_Share_VM/TRMM_3B42RT_daily/2016/10/3B42RT_daily.2016.10.31.bin')

(trmm_bin_array, bbox) = init(trmm_bin_filename, rt=1)
trmm_bin_array = np.asarray(trmm_bin_array)
trmm_bin_array [trmm_bin_array == -9999.9] = -99

trmm_array = trmm['precipitation'][:]
trmm_lat = trmm['lat'][:]
trmm_lon = trmm['lon'][:]
gpm_array = gpm['precipitationCal'][:].data
trmm_fillValue = trmm['precipitation']._FillValue
gpm_fillValue = gpm['precipitationCal']._FillValue


a = trmm_array.repeat(5, axis=0)
b = a.repeat(5, axis=1)

trmm_reshaped = rebin(b, (3600,1200))
trmm_reshaped[trmm_reshaped < 0] = -99
trmm_reshaped[trmm_reshaped == trmm_fillValue] = -99


trmm_array[trmm_array < 0] = -99
gpm_array[gpm_array == gpm_fillValue] = -99

folder = r'/media/sf_Share_VM/Nc4_test'



WriteGTiff_2(np.rot90(trmm_array),folder,'trmm3b42_20161031',-180,180,-60,60,-99)
WriteGTiff_2(np.rot90(trmm_reshaped),folder,'trmmReshape_20161031',-180,180,-60,60,-99)
WriteGTiff_2(trmm_bin_array,folder,'trmm3b42_20161031_binformat',-180,180,-50,50,-99)
WriteGTiff_2(np.rot90(gpm_array),folder,'GPM_20161031',-180,180,-90,90,-99)




