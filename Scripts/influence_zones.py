# ------------------------------------------------------------------------------
#
# Copyright 2011, 2012, 2013 Brent L. Brock and HoloScene Wildlife Services LLC
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
# circuit_mapper.py
# Created on: January 28 2015 03:58:21 PM
# Author: Brent L. Brock
# ---------------------------------------------------------------------------


# Import system modules
import sys, string, os, fileinput, arcpy
import local_params as lp
from functionlib import CreateTempWorkspace, CleanFiles, ConvertMetersToOther, ConvertDistanceToMeters, InvertRaster, Rescale

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
pntLayer = sys.argv[1]
threshDist = sys.argv[2]
viewShed = sys.argv[3]
elevLayer = sys.argv[4]
CellSize = sys.argv[5]
extent = sys.argv[6]
outZones = sys.argv[7]
scriptpath = sys.path[0]

### Set the extent and cell size to link Layer
if extent == "" or extent == "#":
    arcpy.env.extent = arcpy.Describe(pntLayer).extent
    arcpy.AddMessage("extent blank")
else:
    arcpy.env.extent = extent
##if CellSize == "" or CellSize == "#":
##    try:
##        arcpy.env.cellSize = '"%s"' % arcpy.Describe(source).catalogPath
##    except:
##        try:
##            arcpy.env.cellSize = '"%s"' % arcpy.Describe(lcLayer).catalogPath
##        except:
##            arcpy.AddError("Cell size cannot be set. Please specify an output cell size and try again...")
##else:
##    arcpy.env.cellSize = CellSize
##CellSize = str(arcpy.env.cellSize)

arcpy.env.cellSize = CellSize
##arcpy.env.cellSize = round(float(arcpy.env.cellSize), 12)
CellSize = str(arcpy.env.cellSize)

# Get threshold distance
distance = threshDist.split(" ")
threshDist = str(ConvertDistanceToMeters(distance[0], distance[1]))
del distance

### Set path for density reclass table
##dFile = sys.path[0] + os.sep + "ReclassTables" + os.sep + "density.txt"
##
### Allow specification of custom habitat layer
##if not lcSource == "Custom":
##
##    # Locate path for correct reclass tables
##    lcSource = lcSource.replace("'", "")
##    if lcSource == "Montana Landcover;National Landcover Dataset (NLCD)":
##        lcSource = lcSource.split(";")[0]
##        
##    if lcSource == "Montana Landcover":
##        fPath = sys.path[0] + os.sep + "ReclassTables" + os.sep + "MT_landcov"
##        
##    elif lcSource == "National Landcover Dataset (NLCD)":
##        fPath = sys.path[0] + os.sep + "ReclassTables" + os.sep + "NLCD"
##
### Find the correct reclass table
##if spPref == "NO-PREFERENCE" or lcSource == "Custom":
##    rFile = "none"
##elif spPref == "FOREST-SHRUB":
##    rFile = fPath + os.sep + "ForestClasses.txt"    
##elif spPref == "GRASSLAND-OPEN":
##    rFile = fPath + os.sep + "ForestClasses.txt"
##elif spPref == "WETLAND-MARSH":
##    rFile = fPath + os.sep + "WetlandClasses.txt"
##
### Local variables...
##lcProj = "xxlcproj"
        
# Create cost surface
arcpy.AddMessage("Creating cost surface raster...")
desc = arcpy.Describe(pntLayer) #.SpatialReference.LinearUnitName
r = desc.SpatialReference.LinearUnitName
costExtent = desc.extent
if not r == "Meter":
    d = str(ConvertMetersToOther("908", r))
else:
    d = "908"

pntProjName = desc.spatialreference.name

## ------------------------------------------------------------------------------------------------------------------------
## Version 2 modifications. Structure costs calculated as maximum
## cost to threshold distance and then multiplied by inverse distance
## squared beyond threshold distance.

#Calculate threshold buffers and inversed distance
maxThresh = float(threshDist) * 4
arcpy.AddMessage("\tBuffering structures by threshold distance: " + threshDist +"...")
threshRaster = arcpy.sa.EucDistance(pntLayer, float(threshDist), CellSize)
threshRaster = arcpy.sa.Con(threshRaster, 1, "", "Value >= 0")
arcpy.AddMessage("\tCalculating inverse distance from thresholds...")
threshRaster = arcpy.sa.EucDistance(threshRaster, maxThresh, CellSize)
Dmax = arcpy.GetRasterProperties_management (threshRaster, "MAXIMUM").getOutput(0)
threshRaster = arcpy.sa.Divide(threshRaster, float(Dmax))
#threshRaster = Rescale(pntDensity, "in_memory\\euc_rescale.tif", 1, 100)
#threshRaster = arcpy.sa.Times(threshRaster, -1)

# Calculate structure density    
arcpy.AddMessage("\tCalculating structure density...")
pntDensity = arcpy.sa.PointDensity(pntLayer, "NONE", CellSize, "CIRCLE " + d + " MAP", "SQUARE_MILES")
arcpy.AddMessage("\tAssigning house density costs...")

## This is new ##################################################################
## Instead of forcing arbitrary scores to structure density, use log transformed
## structure density as exponent for distance decay function. Areas of high density decay linearly
## while areas of low density decay exponentially

arcpy.AddMessage("\t\tLog transforming structure density...")
pntDensity = arcpy.sa.Log10(pntDensity)
pntDensity = Rescale(pntDensity, "#", 1, 2)
pntDensity = InvertRaster(pntDensity, "#")
pntDensity = arcpy.sa.Con(arcpy.sa.IsNull(pntDensity),2,pntDensity)
pntDensity = arcpy.sa.Negate(pntDensity)

# This is a clunky way to calculate distance with a negative slope m=density index 
pntDensity = arcpy.sa.Times(threshRaster, pntDensity)
pntDensity = arcpy.sa.Plus(pntDensity, 1)
# Truncate values to zero
pntDensity = arcpy.sa.Con(pntDensity,0,pntDensity, "Value < 0")


# Rescale distance values
pntDensity = Rescale(pntDensity,"#", 1, 20)
pntDensity = arcpy.sa.Con(arcpy.sa.IsNull(pntDensity),1,pntDensity)

# Perform viewshed Analysis
if viewShed == 'true':
    arcpy.AddMessage("Initiating viewshed analysis for travel areas...")

    arcpy.AddMessage("\tFinding areas hidden from view...")
    outVShed = arcpy.sa.Viewshed(elevLayer, pntLayer)
    pntDensity = arcpy.sa.Con(outVShed, "1", pntDensity, "Value = 0")

### Make sure inputs are in same projection
##costProjName = arcpy.Describe(pntDensity).spatialreference.name
##if not spPref == "NO-PREFERENCE":
##    lcProjName = arcpy.Describe(lcLayer).spatialreference.name
##
##    if not lcProjName == costProjName:
##        arcpy.AddWarning("Landcover and Linkage Layer projections do not match. Attempting to reproject " + lcLayer + "...")
##        arcpy.ProjectRaster_management(lcLayer, lcProj, pntDensity)
##        lcLayer = lcProj
##
### Create cost surface
##if not rFile == "none":
##    arcpy.AddMessage("\tAssigning habitat costs...")
##    reclass = arcpy.sa.ReclassByASCIIFile(lcLayer, rFile)
##
##    arcpy.AddMessage("\tCombining habitat and house density costs...")
##    
##    xxplus = arcpy.sa.Times(pntDensity, reclass)
##    costRaster = arcpy.sa.Con(pntDensity, pntDensity, xxplus, "Value > 4")
##    
##elif lcSource == "Custom":
##    arcpy.AddMessage("\tApplying custom cost surface...")
##    xxplus = arcpy.sa.Times(pntDensity, lcLayer)
##    costRaster = arcpy.sa.Con(pntDensity, pntDensity, xxplus, "Value > 4")
##    
##else:        
##    arcpy.AddMessage("\tNo habitat preference specified... \n\tCalculating movement resistance from housing density only...")
##    costRaster = pntDensity

# Save output influence zone maps if specified
##if  iZones == "true":
arcpy.AddMessage("\tSaving influence zone raster...")
##outZones = outName + "_iZones.img"
pntDensity = Rescale(pntDensity, "#", 0, 20)
pntDensity.save(outZones)

### Make sure Min value of costRaster > 0
##Habmin = arcpy.GetRasterProperties_management (costRaster, "MINIMUM").getOutput(0)
##if Habmin == 0:
##    arcpy.AddWarning("\tZero values in cost raster found! Adjusting minimum value to 1")
##    costRaster = arcpy.sa.Con(costRaster, "1", costRaster, "Value = 0")
##
### Check costRaster and throw error if values are out of range
##Habmax = arcpy.GetRasterProperties_management (costRaster, "MAXIMUM").getOutput(0)
##if int(float(Habmin)) < 1 or int(float(Habmax)) > 20:
##    arcpy.AddError("Cost raster values are out of range. Check habitat preference input parameters.")
##    arcpy.AddMessage("Min: " + str(Habmin) + ", Max: " + str(Habmax))
##    sys.exit(-1)
##
###Set up inputs for Circuitscape
##arcpy.AddMessage("Preparing Inputs for Circuitscape Modeling...")
##
##desc = arcpy.Describe(source)
##dType = desc.dataType
##sourceName = outName + "_source.asc"
##
### Check if source patches are raster and convert if necessary
##if dType == "ShapeFile" or dType == "FeatureLayer":
##    arcpy.AddMessage("\tConverting Source Patches to ASCII...")
##    arcpy.FeatureToRaster_conversion(source, sourceID, "in_memory\\xxsrcrast", CellSize) # Tried to use .getOutput(0) to set sourceRaster but get inconsistent results and errors.
##    sourceRaster = "in_memory\\xxsrcrast"
##
##    if arcpy.Describe(sourceRaster).SpatialReference.name == 'Unknown':
##        sr = arcpy.Describe(source).SpatialReference
##        arcpy.DefineProjection_management(sourceRaster, sr)
##else:
##    sourceRaster = source
##    cell = arcpy.GetRasterProperties_management(sourceRaster, "CELLSIZEX").getOutput(0)
##
##    if not str(cell) == CellSize:
##        arcpy.AddWarning("\tResampling Cost Surface to " + CellSize + " meter cells...")
##        xxResample = arcpy.Resample_management(sourceRaster, 'in_memory\\xxResample', CellSize, "NEAREST").getOutput(0)
##        sourceRaster = xxResample
##
### Make sure source patches are in correct projection
##SrcPrj = arcpy.Describe(sourceRaster).spatialreference.name
##
##if not SrcPrj == costProjName:
##    arcpy.AddMessage("\tReprojecting Source Patch Layer to match cost raster...")
##    arcpy.ProjectRaster_management (sourceRaster, "in_memory\\xxsrcrast2", costRaster) # Tried to use .getOutput(0) to set sourceRaster but get inconsistent results and errors.
##    sourceRaster = "in_memory\\xxsrcrast2"
##
### Make sure source and cost rasters are same extent    
##SrcExt = arcpy.Describe(sourceRaster).extent
##if not SrcExt == costExtent:
##    arcpy.AddMessage("\tExtracting source raster to cost raster extent...")
##    sourceRaster = arcpy.sa.ExtractByMask (sourceRaster, costRaster)
##
### Convert inputs to ASCII
##arcpy.AddMessage("\tConverting Cost Surface to ASCII...")
##costName = outName + "_cost.asc"
##
##### BUG!!!! This step is removing some nodata cells from the output!!!! Tried setting a snapraster to the input but no luck!!!
##### This step works correctly when used manually from ArcToolbox. 
##arcpy.RasterToASCII_conversion(costRaster, costName)
##arcpy.RasterToASCII_conversion(sourceRaster, sourceName)
##
##arcpy.AddMessage("Performing Circuitscape Modeling...")
###Copy circuitscape.ini file to temp workspace and modify
##arcpy.AddMessage("\tGenerating Circuitscape setup file...")
##cs_ini = os.path.dirname(outName) + os.sep + "cs_ini.ini"
##fileToSearch  = scriptpath + os.sep + "circuitscape.ini"
## 
##tempFile = open( cs_ini, 'w' ) 
## 
##for line in fileinput.input( fileToSearch ):
##    tempFile.write(line.replace( "%NODES%", sourceName ) )
##tempFile.close() 
##
##tempFile = open( cs_ini, 'r+' )
##                                 
##for line in fileinput.input( cs_ini): 
##    tempFile.write( line.replace( "%OUTPUT%", outName ) )
##tempFile.close()
##
##tempFile = open( cs_ini, 'r+' )
##for line in fileinput.input( cs_ini):
##    tempFile.write( line.replace( "%HABITAT%", costName ) )
##tempFile.close()
##
### Run Circuitscape
##arcpy.AddMessage("\tExecuting Circuitscape...")
##cs_ini = cs_ini.replace("\\", "\\\\")
##
### Check location of cs_run.exe and change path if found in a different location
##filename = "cs_run.exe"
##n = False
##
##if os.path.exists(lp.cs_path):
##    arcpy.AddMessage('\tcs_run.exe found...')
##else:
##    arcpy.AddWarning("\tCircuitscape executable not found at default location. \n\tSearching C: drive for file...")
##    for root, dirs, names in os.walk("c:\\"):
##        if filename in names:
##            n = os.path.join(root, filename)
##            n = n.replace('\\', '\\\\')
##            arcpy.AddMessage("\tCircuitscape executable found at: " + n)
##            arcpy.AddMessage("\tUpdating local copy of script with new Circuitscape path...")
##            print n
##
##            # paramFile  = sys.argv[0]
##            old = os.path.dirname(sys.argv[0]) + os.sep + "local_params.py"
##            new = os.path.dirname(sys.argv[0]) + os.sep + "local_params.txt"
##            f = open( old, 'r' ) 
##            temp = open(new, 'w')
##            for line in f:
##                if "cs_path = " in line:
##                    l = line.replace(line, "cs_path = '" + n) + "'\n"
##                else:
##                    l = line
##                print l
##                temp.write(l)
##            f.close()
##
##            temp.close()
##            os.remove(old)
##            os.rename(new, old)
##            break
##    else:
##        arcpy.AddError('Circuitscape executable not found. Make sure Circuitscape is insalled \nor manually edit the cs_path variable in local_params.py')
##        
##inASCII = outName + "_cum_curmap.asc"
##outRaster = outName + "_cum_curmap.img"
##curTemp = os.path.dirname(outName) + os.sep + "curTemp.img"
##if os.path.exists(curTemp):
##    os.remove(curTemp)
##
### Execute Circuitscape
##os.system('"' + lp.cs_path + '" ' + cs_ini)
##
##arcpy.AddMessage("\tSetting Source Areas to NULL...")
##sourceName = arcpy.sa.Con(arcpy.sa.IsNull(sourceRaster),-9999, sourceRaster)
###curTemp = arcpy.sa.Con(sourceName,inASCII, "", "Value = 0")
##curTemp = arcpy.sa.SetNull(sourceName, inASCII, "Value <> -9999")
##
##if arcpy.Describe(curTemp).SpatialReference.name == 'Unknown':
##    sr = arcpy.Describe(costRaster).SpatialReference
##    arcpy.DefineProjection_management(curTemp, sr)
##
### Save output raster
##arcpy.CopyRaster_management(curTemp, outRaster)
##
##arcpy.AddMessage("\tComputing Final Raster Statistics...")
##arcpy.CalculateStatistics_management(outRaster)
##
### delete resistance output file if it exists
##name = outName + "_resistances_3columns.out"
##if os.path.exists(name):
##    os.remove(name)
##n = name.split(".")
##os.rename(n[0], name)
###Add outputs to display
##arcpy.SetParameterAsText(13, outRaster)


arcpy.SetParameterAsText(6, outZones)

# Clean up temporary workspace
##os.remove(cs_ini)
##arcpy.Delete_management(curTemp)
arcpy.Delete_management("in_memory")
