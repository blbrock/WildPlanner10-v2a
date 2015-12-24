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
# Created on: Thu Apr 14 2011 03:58:21 PM
# Author: Brent L. Brock, Craighead Institute
# Usage: corridor_mapper <vshdLayer> <lcLayer> 
# ---------------------------------------------------------------------------


# Import system modules
import sys, string, os, arcpy
from functionlib import CreateTempWorkspace, CleanFiles, ConvertMetersToOther, Rescale

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
costRaster = sys.argv[2]
##lcLayer = sys.argv[2]
##lcSource = sys.argv[3]
##spPref = sys.argv[4]
linkLayer = sys.argv[3]
direction = sys.argv[4]
outRaster = sys.argv[5]

if direction == "Any Direction":
    direction = "East-West;North-South"

### Set up temporary workspace
##scratchWS = arcpy.env.scratchWorkspace
##if scratchWS:
##    tWorkspace = CreateTempWorkspace(scratchWS)
##else:
##    # Set output workspace
##    r = outRaster.rsplit(os.sep,1)
##    outWorkspace = r[0]
##    arcpy.Workspace = outWorkspace
##    tWorkspace = CreateTempWorkspace(outWorkspace)

# Set the extent to link Layer
# arcpy.Extent = arcpy.Describe(linkLayer).extent

arcpy.env.extent = arcpy.Describe(costRaster).catalogPath
arcpy.env.cellSize = arcpy.Describe(costRaster).catalogPath
CellSize = str(arcpy.env.cellSize)
##CostCellSize = arcpy.GetRasterProperties_management(costRaster, "CELLSIZEX")
##arcpy.AddMessage("CellSize = " + CellSize + ", CostCellSize = " + str(CostCellSize))
##if not arcpy.env.cellSize == CostCellSize:
##    arcpy.AddMessage("Resampling cost raster to match analysis cell size...")
##    costRaster = arcpy.Resample_management (costRaster, arcpy.env.scratchGDB + os.sep + "costRaster", CellSize, "NEAREST")
##
### Set path for density reclass table
##dFile = sys.path[0] + os.sep + "ReclassTables" + os.sep + "density.txt"
##
### Locate path for correct reclass tables
##lcSource = lcSource.replace("'", "")
##if lcSource == "Montana Landcover;National Landcover Dataset (NLCD)":
##    lcSource = lcSource.split(";")[0]
##    
##if lcSource == "Montana Landcover":
##    fPath = sys.path[0] + os.sep + "ReclassTables" + os.sep + "MT_landcov"
##    
##elif lcSource == "National Landcover Dataset (NLCD)":
##    fPath = sys.path[0] + os.sep + "ReclassTables" + os.sep + "NLCD"
##
### Find the correct reclass table
##if spPref == "FOREST-SHRUB":
##    rFile = fPath + os.sep + "ForestClasses.txt"
##elif spPref == "GRASSLAND-OPEN":
##    rFile = fPath + os.sep + "ForestClasses.txt"
##elif spPref == "WETLAND-MARSH":
##    rFile = fPath + os.sep + "WetlandClasses.txt"

# Local variables...
reclass = "in_memory" + os.sep + "xxrcls"
east_cst = "in_memory" + os.sep + "xxeastcst"
west_cst = "in_memory" + os.sep + "xxwestcst"
north_cst = "in_memory" + os.sep + "xxnorthcst"
south_cst = "in_memory" + os.sep + "xxsouthcst"
corTempNS = "in_memory" + os.sep + "xxcorns"
corTempEW = "in_memory" + os.sep + "xxcorew"
resclNS = "in_memory" + os.sep + "xxcornsr"
resclEW = "in_memory" + os.sep + "xxcorewr"
corTemp = "in_memory" + os.sep + "xxcor"
corTemp2 = "in_memory" + os.sep + "xxcor2"
lcProj = "in_memory" + os.sep + "xxlcproj"

# Make sure inputs are in same projection

##lcProjName = arcpy.Describe(lcLayer).spatialreference.name
##linkProjName = arcpy.Describe(linkLayer).spatialreference.name
##
##if not lcProjName == linkProjName:
##    arcpy.AddWarning("Landcover and Linkage Layer projections do not match. Attempting to reproject " + lcLayer + "...")
##    arcpy.ProjectRaster_management(lcLayer, lcProj, linkLayer)
##    lcLayer = lcProj

### Create cost surface
##arcpy.AddMessage("Creating cost surface raster...")
##r = arcpy.Describe(pntLayer).SpatialReference.LinearUnitName
##if not r == "Meter":
##    d = str(ConvertMetersToOther("908", r))
##else:
##    d = "908"
##    
##arcpy.AddMessage("\tCalculating structure density...")
##pntDensity = arcpy.sa.PointDensity(pntLayer, "NONE", "30", "CIRCLE " + d + " MAP", "SQUARE_MILES")
##arcpy.AddMessage("\tAssigning house density costs...")
##pntDensity = arcpy.sa.ReclassByASCIIFile(pntDensity, dFile)
###
##pntDensity.save(tWorkspace + os.sep + "xxdensity")
### Create cost surface
##arcpy.AddMessage("\tAssigning habitat costs...")
##reclass = arcpy.sa.ReclassByASCIIFile(lcLayer, rFile)
###
##reclass.save(tWorkspace + os.sep + "xxrcls")
##
##
##arcpy.AddMessage("\tCombining habitat and house density costs...")
##xxplus = arcpy.sa.Plus(pntDensity, reclass)
###
##reclass.save(tWorkspace + os.sep + "xxplus")
##
##xxCombine = arcpy.sa.Con(xxplus, "50", xxplus, "Value > 50")
##
####maxval = arcpy.GetRasterProperties_management (xxplus, "MAXIMUM")
####maxval = float(str(maxval))    
####maxcost = maxval * 10
##costRaster = arcpy.sa.Con(linkLayer, "50", xxplus, "Value = 0")
##costRaster.save(tWorkspace + os.sep + "costraster1") 

# Create source and destination polygons (NSEW)
arcpy.AddMessage("Creating corridor destination polygons...")
maxval = arcpy.GetRasterProperties_management (linkLayer, "MAXIMUM")
##xxblank = arcpy.sa.Reclassify(linkLayer, "VALUE", "0 " + str(maxval) + " 1", "DATA")
if arcpy.GetRasterProperties_management(costRaster,"VALUETYPE") > 8:
    costLayer = arcpy.sa.Int(costRaster)
else:
    costLayer = costRaster
    
costBndry = arcpy.RasterToPolygon_conversion(costLayer, "in_memory" + os.sep + "xxpoly", "SIMPLIFY", "VALUE")
costBndry = arcpy.Dissolve_management(costBndry, "in_memory" + os.sep + "cb_dissolved","", "", "SINGLE_PART","")
# Createfeature envelope of costRaster
costEnv = arcpy.MinimumBoundingGeometry_management(costBndry, "in_memory" + os.sep + "cb_env", "ENVELOPE", "ALL")
# Convert costRaster boundary to line
costBndry = arcpy.FeatureToLine_management(costBndry, "in_memory" + os.sep + "xxline", "", "NO_ATTRIBUTES")
costBndry = arcpy.Buffer_analysis(costBndry, "in_memory" + os.sep + "buf", CellSize + " Meters", "FULL", "ROUND", "NONE", "")
costEnv = arcpy.FeatureToLine_management(costEnv, "in_memory" + os.sep + "xxenv", "", "NO_ATTRIBUTES")
EnvPnt = arcpy.FeatureVerticesToPoints_management(costEnv, "in_memory" + os.sep + "xxpnt")
fc = arcpy.MakeFeatureLayer_management(EnvPnt, "pLayer")
oid_fieldname = arcpy.Describe(fc).OIDFieldName
arcpy.SelectLayerByAttribute_management(fc, "NEW_SELECTION", oid_fieldname + " = 1 OR " + oid_fieldname + " = 3")
diag_1 = arcpy.PointsToLine_management(fc, "in_memory" + os.sep + "diag_1")
arcpy.SelectLayerByAttribute_management(fc, "SWITCH_SELECTION")
diag_2 = arcpy.PointsToLine_management(fc, "in_memory" + os.sep + "diag_2")
costEnv = arcpy.Merge_management([costEnv, diag_1, diag_2], "in_memory" + os.sep + "merge_line")
costEnv = arcpy.FeatureToPolygon_management(costEnv, "in_memory" + os.sep + "xxenv", "", "NO_ATTRIBUTES")
costSource = arcpy.Intersect_analysis([costBndry, costEnv], "in_memory" + os.sep + "sources")



##tmpExtent = arcpy.env.extent
##arcpy.env.extent = 'MAXOF'
##
######
##fc = arcpy.Buffer_analysis(xxblank, "in_memory" + os.sep + "buf", CellSize + " Meters", "RIGHT", "FLAT", "NONE", "")
##fc = arcpy.MakeFeatureLayer_management(fc, "pLayer")
### arcpy.CopyFeatures_management (fc, arcpy.env.scratchGDB + "\\buffLayer")
##
### rows = arcpy.SearchCursor(fc)
##
##
##
##with arcpy.da.SearchCursor(fc, oid_fieldname) as cursor:
##
### for row in rows:
##    for row in cursor:
##        # p = row.getValue("OBJECTID")
##        p = row[0]
##        #p = row.getValue(oid_fieldname)
##        # arcpy.MakeFeatureLayer_management(fc, "fLayer")
##        # arcpy.CopyFeatures_management (fc, arcpy.env.scratchGDB + "\\fLayer")
##        arcpy.SelectLayerByAttribute_management(fc, "NEW_SELECTION", oid_fieldname + " = " + str(p))
##        # pLayer = arcpy.env.scratchGDB + "\\pLayer" + str(p)
##        # arcpy.AddMessage("pLayer = " + pLayer)
##        # arcpy.CopyFeatures_management (fc, arcpy.env.scratchGDB + "\\fLayer1")
##        # arcpy.AddMessage("number selected: " + str(int(arcpy.GetCount_management(fc).getOutput(0))))
##        
##        # arcpy.AddMessage("pLayer = " + pLayer)
##        
##        pLayer = arcpy.MakeFeatureLayer_management(fc, "\\pLayer" + str(p))
##        arcpy.CopyFeatures_management(pLayer, arcpy.env.scratchGDB + os.sep + "pLayerfc" + str(p))
##        if p == 1:
##            west = pLayer
##        elif p == 2:
##            north = pLayer
##        elif p == 3:
##            east = pLayer
##        elif p == 4:
##            south = pLayer
##        else:
##            arcpy.AddError("Origin layers out of range. Check geometry of linkage layer...")
##            arcpy.Delete_management("in_memory")
##            sys.exit(0)
##arcpy.AddMessage("east = " + str(east))
##    
####west = pLayer1
####north = pLayer2
####east = pLayer3
####south = pLayer4
##
###################
##try:
##    arcpy.CopyFeatures_management (west, arcpy.env.scratchGDB + os.sep + "west")
##    arcpy.CopyFeatures_management (north, arcpy.env.scratchGDB + os.sep + "north")
##    arcpy.CopyFeatures_management (east, arcpy.env.scratchGDB + os.sep + "east")
##    arcpy.CopyFeatures_management (south, arcpy.env.scratchGDB + os.sep + "south")
##except:
##    pass
##
##del row
##del cursor
##del fc
##
##arcpy.env.extent = tmpExtent

# Run cost-distance analysis
arcpy.AddMessage("Creating cost-distance rasters...")
sources = arcpy.MakeFeatureLayer_management(costSource, "sources")
oid_fieldname = arcpy.Describe(sources).OIDFieldName
# arcpy.ClearEnvironment("extent")
if direction == "East-West" or direction == "East-West;North-South":
    ##### ????
    # originEW = "in_memory" + os.sep + "xxoriginEW"
    # arcpy.Merge_management([east, west], originEW)
    arcpy.SelectLayerByAttribute_management(sources, "NEW_SELECTION", oid_fieldname + " = 2")
    east_cst = arcpy.sa.CostDistance(sources, costRaster)
    arcpy.SelectLayerByAttribute_management(sources, "NEW_SELECTION", oid_fieldname + " = 4")
    west_cst = arcpy.sa.CostDistance(sources, costRaster)
    arcpy.AddMessage("Modeling east-west linkage zones...")
    corTempEW = arcpy.sa.Corridor(east_cst, west_cst)
    ##########
    corTempEW.save(arcpy.env.scratchGDB + "\\corTempEW")
##    east_cst.save(arcpy.env.scratchGDB + "\\east_cst")
##    west_cst.save(arcpy.env.scratchGDB + "\\west_cst")
    if not direction == "East-West;North-South":
        corTemp = corTempEW
        # origin = originEW
        
if direction == "North-South" or direction == "East-West;North-South":
    ##### ????
##    originNS = "in_memory" + os.sep + "xxoriginNS"
##    arcpy.Merge_management([north, south], originNS)
    arcpy.SelectLayerByAttribute_management(sources, "NEW_SELECTION", oid_fieldname + " = 1")
    south_cst = arcpy.sa.CostDistance(sources, costRaster)
    arcpy.SelectLayerByAttribute_management(sources, "NEW_SELECTION", oid_fieldname + " = 3")
    north_cst = arcpy.sa.CostDistance(sources, costRaster)
    arcpy.AddMessage("Modeling north-south linkage zones...")
    corTempNS = arcpy.sa.Corridor(north_cst, south_cst)
    ##########
    corTempNS.save(arcpy.env.scratchGDB + "\\corTempNS")
    if not direction == "East-West;North-South":
        corTemp = corTempNS
        # origin = originNS

if direction == "East-West;North-South":
    arcpy.AddMessage("Combining linkage zones...")
##    corTempEW.save(arcpy.env.scratchGDB + os.sep + "xxcorTempEW")
##    corTempEW = "in_memory" + os.sep + "xxcorTempEW"
##    # corTempNS.save(tWorkspace + os.sep + "xxcorTempNS")
##    corTempNS = "in_memory" + os.sep + "xxcorTempNS"
##    resclEW = arcpy.env.scratchGDB + os.sep + "rescale_EW"
##    resclNS = arcpy.env.scratchGDB + os.sep + "rescale_NS"
    arcpy.AddMessage("Attempting to rescale " + arcpy.Describe(corTempEW).baseName)
    #Type = arcpy.GetRasterProperties (corTempEW, "MINIMUM")
    resclEW = Rescale(corTempEW, "in_memory" + os.sep + "rescale_EW", 0, 100)
    resclEW.save(arcpy.env.scratchGDB + os.sep + "rescale_EW_save")
    arcpy.AddMessage("Attempting to rescale " + arcpy.Describe(corTempNS).baseName)
##    Rescale(corTempNS, resclNS, 0, 100)
    resclNS = Rescale(corTempNS, "in_memory" + os.sep + "rescale_NS", 0, 100)
    corTemp = arcpy.Mosaic_management(resclEW, resclNS, "MINIMUM")
    resclNS.save(arcpy.env.scratchGDB + os.sep + "core_Temp")
    
    corTemp = resclNS
    # origin = arcpy.Merge_management([originEW, originNS], "in_memory" + os.sep + "xxorigin")

corTemp = arcpy.sa.SetNull(costRaster, corTemp, "Value = 0")
corTemp.save(outRaster)
##sliceRaster = arcpy.sa.Slice(corTemp, 5, "NATURAL_BREAKS")
##
### Create connected area mask
##arcpy.AddMessage("Creating connected area mask to remove isolated habitat patches...")
##
###
### sliceRaster.save(arcpy.env.scratchGDB + "sliceRaster")
##movementAreas = arcpy.sa.SetNull(linkLayer, linkLayer, "Value = 0")
###
##movementAreas.save(arcpy.env.scratchGDB + os.sep + "movement")
##region = arcpy.sa.RegionGroup(movementAreas, "EIGHT")
###
### region.save(tWorkspace + os.sep + "region")
##arcpy.RasterToPolygon_conversion(region, "in_memory" + os.sep + "xxregion", "NO_SIMPLIFY", "VALUE")
##arcpy.MakeFeatureLayer_management("in_memory" + os.sep + "xxregion", 'region_lyr')
###############################
##arcpy.SelectLayerByLocation_management('region_lyr', 'intersect', sources)
##arcpy.MakeFeatureLayer_management('region_lyr', 'connected')
##arcpy.AddMessage("Removing compromised areas from analysis area...")
##corNull = arcpy.sa.SetNull(linkLayer, corTemp, "Value = 0")
###
### corNull.save(tWorkspace + os.sep + "corNull")
##arcpy.AddMessage("Removing isolated patches from analysis area...")
##corTemp2 = arcpy.sa.ExtractByMask(corNull, "connected")
###
##corTemp2.save("in_memory" + os.sep + "corTemp2")
##arcpy.AddMessage("Ranking linkage zones...")
##sliceRaster = arcpy.sa.Slice(corTemp2, 5, "NATURAL_BREAKS")
##conRaster = arcpy.sa.Con(linkLayer, "0", sliceRaster, "Value = 0")
##outIsNull = arcpy.sa.IsNull(conRaster)
##rankRaster = arcpy.sa.Con(outIsNull, "0", conRaster, "Value = 1")
##rankRaster.save(outRaster)


#Add outputs to display
arcpy.SetParameterAsText(4, outRaster)

# Clean up temporary workspace
arcpy.Delete_management("pLayer")
arcpy.Delete_management("sources")
##arcpy.Delete_management(south)
##arcpy.Delete_management(west)
##arcpy.Delete_management(north)
##arcpy.Delete_management(east)
# CleanFiles(tWorkspace)
arcpy.Delete_management("in_memory")










    
    






