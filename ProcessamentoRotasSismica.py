
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon,LineString,asLineString,MultiLineString,Point
import pandas as pd

def CreateLines(
				AreaFile,
				OutputName='Linhas.shp',
				MajorAxis=50,
				MinorAxis=30,
				LineName=None
				):
	'''
	Creates a shapefile with sample lines
	
	AreaFile: str
		Shapefile name with area limits. This shapefile must be a Polygon or Point (Polygon preference)
	
	OutputName: str default 'Linhas.shp'
		String with name of outputfile if it doesn't have the file extention .shp, it will create a directory 
		with the files inside
	
	MajorAxis: float, int, default 50
		Distance in meters between major axis' lines. 
	
	MinorAxis: float, int, default 50
		Distance in meters between minor axis' lines 
	
	LineName: str, None default is None	
		Name of shapefile with major axis. If it's None the definition would calculates the major line based 
		in farest point of area.
	'''		
	# open file
	gdf = gpd.read_file(AreaFile)
				
	# checks
	if not isinstance(gdf.geometry[0], Polygon):
		# creates geometry
		geometry=Polygon([[p.x, p.y] for p in gdf.geometry])
		# gets x and y values
		x,y = geometry.exterior.coords.xy
		# Calculates polygon center
		xc = sum(x)/len(x)
		yc = sum(y)/len(y)
		# Calculates difference
		xt=np.asarray(x)-xc
		yt=np.asarray(y)-yc
		# Creates sorter flags
		gdf['sorter'] = np.arctan2(yt[:-1],xt[:-1])
		# Sorts values
		ngdf = gdf.sort_values('sorter')
		# Recreates geometry
		geometry=Polygon([[p.x, p.y] for p in ngdf.geometry])
		# Creates Polygon shape
		area=gpd.GeoDataFrame(['area'],geometry=[geometry],crs=gdf.crs)
		
	else:
		area=gdf

	# gets the polygon centers
	xc = area.centroid.x[0]
	yc = area.centroid.y[0]
	
	# gets coordinates 
	coords = area.geometry.apply(lambda x:x.exterior.coords.xy)
	x,y = coords.values[0]
	
	if not LineName:
		# calcualtes farest point from center
		distance = ((np.asarray(x)-xc)**2+(np.asarray(y)-yc)**2)**0.5
		# calculates ind
		ind = np.argmax(distance)
		# calculates individual distance
		dx = x[ind]-xc
		dy = y[ind]-yc
		# creates line
		line = LineString(
			[(xc-dx,yc-dy),
			(xc+dx,yc+dy)]
		)
		
		line = gpd.GeoDataFrame(['main'],geometry=[line],crs=gdf.crs)
	else:
		line=gpd.read_file(LineName)
	# interesects with polygon
	line = area.intersection(line)

	# calculates minot linecache
	# calculates line center
	xLineC = line.centroid.x
	yLineC = line.centroid.y
	# gets lne values
	minorX,minorY = line.geometry.values[0].coords.xy
	# calcunates new coordinates
	newx = (-np.asarray(minorY) + yLineC.values)+xLineC.values
	newy = (np.asarray(minorX) - xLineC.values)+yLineC.values
	# creates minor line
	minorLine = LineString([(xx,yy) for xx,yy in zip(newx,newy) ])
	minorLine = gpd.GeoDataFrame(geometry=[minorLine],crs=gdf.crs)
	minorLine = area.intersection(minorLine)
	
	# apply buffers
	# area bounds
	bnds = area.bounds
	# calculates number of lines
	if dy>dx:
		N = np.ceil(len(np.arange(bnds['minx'],bnds['maxx'],MajorAxis))/2)
		M = np.ceil(len(np.arange(bnds['miny'],bnds['maxy'],MinorAxis))/2)
			
	else:
		N = np.ceil(len(np.arange(bnds['miny'],bnds['maxy'],MajorAxis))/2)
		M = np.ceil(len(np.arange(bnds['minx'],bnds['maxx'],MinorAxis))/2)
    
	# creates buffers
	Majorbuffers=np.arange(MajorAxis,N*MajorAxis,MajorAxis)
	Minorbuffers=np.arange(MinorAxis,M*MinorAxis,MinorAxis)

	# definition tu buffer
	def LineBuffer(line,brange):

		nline = asLineString(line.geometry.values[0].buffer(brange,cap_style=3).exterior)
		
		return nline

	# creates Multilinestring 
	MajorGeometry = [MultiLineString([LineBuffer(line,buf) for buf in Majorbuffers])+line]
	MinorGeometry = [MultiLineString([LineBuffer(minorLine,buf) for buf in Minorbuffers])+minorLine]
	Geometry=MajorGeometry+MinorGeometry
	# creates geodataframes
	Lines = gpd.GeoDataFrame(['Major Lines','Minor Lines'],geometry=Geometry,crs=gdf.crs)
	Lines = area.intersection(Lines)
	# write file
	Lines.to_file(OutputName)
	
def CreateAreaFromCSV(fname,epsg):
	'''
	Creates a shapefile wi
	'''
	# loads file
	df = pd.read_csv(fname)
	# creates geometry
	df['geometry'] = df.apply(Point(df.x,df.y))
	# generates gdf file
	gdf=gpd.GeoDataFrame(df,crs={'init':'epsg:'+str(epsg)})
	
	return gdf

	
if __name__=='__main__':
	
	fname = 'c:/Users/tito.conte/Documents/OneDrive for Business/ferramentas_scripts/Sismica/BarraRiacho-Limites_SIRGAS24S.shp'
	
	CreateLines(fname)