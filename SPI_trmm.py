"""
Author: Walther C.A. Camaro G. : walther.camaro@gmail.com
"""

import process_trmm_3b42_daily
from process_trmm_3b42_daily import init,  WriteGTiff_3, cumulatedict, media, medialn, probnorain, SPI_DICT, calendardays
import datetime
from datetime import timedelta, date
import numpy as np
import os


months = range(1,3)
folder_root = r'D:\Python_Projects\Data\TRMM\TRMM_3B42_daily'  #Original Not Real Time DataSet
folder_rootRT = r'D:\Python_Projects\Data\TRMM\TRMM_3B42RT_daily' #Original Real Time Dataset
cumulateds = [1,3]
for cumulated in cumulateds: # Number of months to cumulate; i.e: 1 Monthly_Spi; 3 Tri_Monthly_SPI; 6 Semestral_SPI...

    c0 = 2.515517
    c1 = 0.802853
    c2 = 0.010328
    d1 = 1.432788
    d2 = 0.189269
    d3 = 0.001308

    for dat in months:
        dict_matrix = {}
        dict_ln = {}
        for y in range(1998,2019):
            ##end_date = datetime.date(int(y),dat,15) #X-th day of the month datetime.date(int(y),dat,X)
            end_date = datetime.date(int(y), dat, calendardays(int(y), dat))  # end of the month
            if end_date.year < 1998:
                print "No %s" % end_date
                continue 
            if end_date.year == 1998 and end_date.month < cumulated:
                print "No %s" % end_date
                continue
            initial_date = (date.fromordinal(date.toordinal(end_date) - (30*cumulated)+1))
            if initial_date.year < 1998:
                print "No %s" % initial_date
                continue
            if len(str(end_date.month))==2:
                end_month = str(end_date.month)
            else:
                end_month = '0%s' % str(end_date.month)    
            if len(str(end_date.day))==2:
                    end_day = str(end_date.day)
            else:
                    end_day = '0%s' % str(end_date.day)
            if date.toordinal(end_date) > date.toordinal(date.today()):
                print initial_date, 'NO - end_date'
                continue
            dict_cumulate = {}
            for kcum in range(date.toordinal(initial_date),date.toordinal(end_date)+1):
                year = str(date.fromordinal(kcum).year)
                if len(str(date.fromordinal(kcum).month))==2:
                    month = str(date.fromordinal(kcum).month)
                else:
                    month = '0%s' % str(date.fromordinal(kcum).month)
                if len(str(date.fromordinal(kcum).day))==2:
                    day = str(date.fromordinal(kcum).day)
                else:
                    day = '0%s' % str(date.fromordinal(kcum).day)
                if os.path.exists(r'%s\%s\%s\3B42_daily.%s.%s.%s.7.bin' % (folder_root, year, month, year, month, day)):
                    filename = r'%s\%s\%s\3B42_daily.%s.%s.%s.7.bin' % (folder_root, year, month, year, month, day)
                    RT = None
                elif os.path.exists(r'%s\%s\%s\3B42RT_daily.%s.%s.%s.bin' % (folder_rootRT, year, month, year, month, day)):
                    filename = r'%s\%s\%s\3B42RT_daily.%s.%s.%s.bin' % (folder_rootRT, year, month, year, month, day)
                    RT = 1
                elif os.path.exists(r'%s\%s\%s\3B42_Daily.%s%s%s.7.nc4' % (folder_root, year, month, year, month, day)):
                    filename = r'%s\%s\%s\3B42_Daily.%s%s%s.7.nc4' % (folder_root, year, month, year, month, day)
                    RT = 0
                elif os.path.exists(r'%s\%s\%s\3B42RT_Daily.%s%s%s.7.nc4' % (folder_rootRT, year, month, year, month, day)):
                    filename = r'%s\%s\%s\3B42RT_Daily.%s%s%s.7.nc4' % (folder_rootRT, year, month, year, month, day)
                    RT = 1
                else: 
                    filename = None
                    continue
                (cumulatematrix,bbox) = init(filename, rt=RT)
                cumulatematrix = np.asarray(cumulatematrix)
                cumulatematrix [cumulatematrix == -9999.9] = -99
                dict_cumulate[kcum] = cumulatematrix
            matrix = cumulatedict(dict_cumulate)
            dict_matrix[str(y)] = matrix
            dict_ln[str(y)] = np.log(dict_matrix[str(y)])
            dict_ln[str(y)][dict_matrix[str(y)]<=0] = 0
        if dict_matrix == {}:
            continue
        else:
            (media_01,n_01) = media(dict_matrix) # media without zeros.
            media_ln01 = medialn(dict_ln,n_01) # media ln without zeros.
            lnmedia = np.log((media_01)) # ln(media)
            lnmedia[media_01 == 0] = 0 # Giving a 0 where there isn't rain in all the period
            A = (lnmedia - media_ln01) 
            A[A < 0] = 0 # Correcting small values of rain in all the season.
            a = 1/(4*A)*(1+(1+(4*A/3))**0.5) # Shape Parameter using Maximum Likelihood.
            a[A == 0] = 0 # Correcting small values of rain in all the season.
            B = media_01/a # Scale parameter using Maximum Likelihood.
            B[a == 0] = 0 # Correcting small values of rain in all the season.
            P_NoRain = probnorain(dict_matrix)
            dict_SPI = SPI_DICT(dict_matrix,P_NoRain,a,B,c0,c1,c2,d1,d2,d3)
            folder_output = r'D:\Python_Projects\Results\TRMM\SPI_2017\%s' % (cumulated)
            if not os.path.exists(folder_output):
                os.makedirs(folder_output)
            WriteGTiff_3(dict_SPI,folder_output,'SPI_%smonthly' % cumulated,bbox['left'],bbox['right'],bbox['bottom'],bbox['top'],dat,end_date.day)
print 'Thank You, Procedure Finished. By W.Camaro'