# ------------------------------------------------------------------------------
#
# Copyright 2011, 2012, 2013 Brent L. Brock and the Craighead Institute
#
# This file is part of Wild Planner.
#
# Wild Planner is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Wild Planner is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Wild Planner in the file named LICENSE.TXT.  If not, see <http://www.gnu.org/licenses/>.
#
# ------------------------------------------------------------------------------
#
# corridor_mapper.py
# Created on: Jan 15 03:58:21 PM
# Author: Brent L. Brock, HoloScene Wildlife Services LLC
# Usage: corridor_mapper <vshdLayer> <lcLayer> 
# ---------------------------------------------------------------------------


# Import system modules
import sys, string, os, arcpy
from arcpy.sa import *
# from functionlib import CreateTempWorkspace, CleanFiles, ConvertMetersToOther, Rescale

# Create the Geoprocessor object

# Set the necessary product code
arcpy.SetProduct("ArcInfo")
arcpy.overwriteoutput = 1

arcpy.AddMessage("\tChecking availability of spatial analyst license...")  
if arcpy.CheckExtension("spatial") == "Available":
    arcpy.CheckOutExtension("spatial")         
else:
    arcpy.AddWarning("\tSpatial Analyst is needed to run this tool.  Processing aborted...")
    raise "LicenseError"    

# Script arguments...
raster = sys.argv[1]
outRaster = sys.argv[2]

outras = Con(IsNull(raster),0,raster)
outras.save(outRaster) 

#Add outputs to display
arcpy.SetParameterAsText(1, outRaster)









    
    






