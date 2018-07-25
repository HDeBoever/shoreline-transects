'''
Define shoreline transects by computing the depths along a midpoint normal for each reach
2018-07-09: For use in shoreline modelling by Cascadia Coast Research

Get a google earth .kmz path file to represet the reaches.
Convert the .kmz file to a UTM .txt file using this website: http://zonums.com/online/kml2x/

Program Usage: path/to/program/python3.6 reach_transects.py (UTM info file).txt

After the initial run, the program will create text files in the transect_files folder. One file pertaining to each reach.

The input text file and the orignal KMZ file should be formulated in  such a way that the reaches go  in a clockwise direction around
the area of interest. This ensures that the transects will be drawn from the sea inwards, rather than from the land outwards.

See shoreline_transectsREADME.txt for more information about the program, its functions, and limitations.

Author: Henri De Boever
'''

from scipy.interpolate import LinearNDInterpolator
from mpl_toolkits.basemap import Basemap
from collections import OrderedDict
from pyproj import Proj, transform
import matplotlib.pyplot as plt
import numpy as np, simplekml
import sys, os, urllib

class Transect():

	# Class initializer
	#
	# The attributes below are specific to each instance of a Reaches class object
	# This program currently handles coordinates in UTM coords in zone 10U
	def __init__(self, bathymetric_data_file, reaches_file):
		print("\nInstantiating a new Transect() abstraction layer.\n")
		self.data = OrderedDict()
		self.npz_file = np.load(bathymetric_data_file)
		self.reaches_file = reaches_file

		# Create the subdirectories in case they don't already exist
		if(os.path.exists('transect_files') is False):
			os.mkdir('transect_files')
			print("Created transect_files directory.")
		if(os.path.exists('reach_files') is False):
			os.mkdir('reach_files')
			print("Created reach_files directory. Please ensure that the .txt file containing the UTM Reaches is stored here as well.")
			sys.exit(0)

	# Class functions
	#
	# The text file is generated from a .kmz google earth path file on http://zonums.com/online/kml2x/
	# Read a text file containing the coords of each node that forms the reaches
	# Then create a list of elements containing 2 nodes which define a reach
	def read_reach_file(self):
		print("Reading reach file: %s\n" % self.reaches_file)
		if(os.path.exists('reach_files/' + self.reaches_file) is False):
			print("No such file in Directory. Place a file containing relevant reach data in %s/reach_files" % os.getcwd())
			sys.exit(0)
		with open('reach_files/' + self.reaches_file, 'r') as f:
			reach_coords_list = []
			# Skip the first 6 lines of the file
			node_data = f.readlines()[6:]
			# Set the current coord
			curr_coord = (float(node_data[0].split(' ')[ :-1][0]), float(node_data[0].split(' ')[ :-1][1]))
			for index, line in enumerate(node_data[1:]):
				next_coord = (float(line.split(' ')[ :-1][0]), float(line.split(' ')[ :-1][1]))
				reach_coords_list.append((curr_coord, next_coord))
				curr_coord = next_coord
		return reach_coords_list

	# Given a list of reahces, obtain the transects associated with each reach
	def get_transects(self, reach_coords, norm_length, epsilon, num_points, upper_elevation_limit, lower_elevation_limit, transect_plot_flag):
		if(transect_plot_flag is True):
			print("transect_plot_flag is turned on\n")
		else:
			print("transect_plot_flag is turned off\n")
		# Create an instance of the interpolator
		interpolator = self.instantiate_interpolator(self.npz_file)
		# Append the column names to the dict
		self.data[0] = ("reach", "normal", "midpoint_coord", "reach_length", "transect_UTM_coords")
		# The data list will contain coordinate tuples which will define each reach.
		# e.g. ((x1, y1),(x2, y2)) is an element in the list. The coordinates define a line on the plane.
		# Find the reach normal: get the coords for it.
		for index, reach in enumerate(reach_coords):
			print("\n----------Transect #%d----------\n" % (index + 1))
			print("Reach: " + str(reach))
			# compute the reach length, reach midpoint, and normal endpoints
			reach_length = self.get_reach_length(reach)
			midpoint_coord = self.get_reach_midpoint(reach)
			in_normal_coords = self.get_reach_normal(reach, midpoint_coord, norm_length/2, 'in')
			out_normal_coords = self.get_reach_normal(reach, midpoint_coord, norm_length/2, 'out')

			# Combine both norms into a single one
			normal = (in_normal_coords[0], out_normal_coords[0])
			# Generate points along the normal to use for the interpolation
			normal_points, facing_shore = self.get_points_along_norm(normal, num_points, reach)
			# Get the z value associated with each point
			z_values = self.get_z_values(normal_points, self.npz_file, interpolator)

			# Associate each z value with a corresponding horizontal value on the line
			points_distances = self.create_points_list(norm_length, z_values, norm_length/num_points, upper_elevation_limit, lower_elevation_limit, facing_shore)
			print("Initial List: %d" % len(points_distances))

			# Perform the Ramer Douglas Peucker Algorithm on the z values params: (list of points with associated distances, epsilon)
			filtered_list = self.ramer_douglas_peucker(points_distances, epsilon)
			print("Filtered List after RDP: %d" % len(filtered_list))
			if(facing_shore is False):
				# Flip the distances in the tuples
				distances = []
				final_list = []
				for tuple in filtered_list[::-1]:
					distances.append(tuple[0])
				for i, tuple in enumerate(filtered_list):
					final_list.append((distances[i], tuple[1]))
				if(transect_plot_flag is True):
					self.plot_list(final_list, (index+1))
			else:
				if(transect_plot_flag is True):
					self.plot_list(filtered_list, (index+1))
				pass
			# For each point in the filtered list, using the midpoints and distance from the midpoint, determine the coords (UTM) of each point
			# transect_UTM_coords will contain the UTM coordinates of each of the points in filtered list associated with depth at that point.
			transect_UTM_coords = self.get_UTM_coords(filtered_list, normal, midpoint_coord, facing_shore, norm_length)
			# Append all the computed values into the self.data dict
			self.data[(index+1)] = (reach, normal, midpoint_coord, reach_length, transect_UTM_coords)
			print("Finished generating transects for area of interest.\n")

	# Called by get_transects() to populate baseline parameters. Returns the length of the reach
	def get_reach_length(self, reach):
		print("Getting reach length...")
		x1, x2, y1, y2 = reach[0][0], reach[1][0], reach[0][1], reach[1][1]
		# Euclidian distance will be accurate enough for our purposes, no need for haversine or vincenty formulae
		length = np.sqrt((x1-x2)**2 + (y1-y2)**2)
		return length

	# Called by get_transects() to populate baseline parameters. Returns the midpoint_coords of a reach
	def get_reach_midpoint(self, reach):
		print("Getting reach midpoint...")
		x1, x2, y1, y2 = reach[0][0], reach[1][0], reach[0][1], reach[1][1]
		# Application of the midpoint of a line formula
		midpoint_coord = (((x1 + x2)/2), ((y1 + y2)/2))
		return midpoint_coord

	# Called by get_transects() to populate baseline parameters. Returns 2 nodes on either side of the reach line
	# Perpendicular bisector formula
	# https://math.stackexchange.com/questions/306468/perpendicular-line-passing-through-the-midpoint-of-another-line
	def get_reach_normal(self, reach, midpoint, distance, direction):
		x1, x2, y1, y2 = reach[0][0], reach[1][0], reach[0][1], reach[1][1]
		reach_slope = (y2-y1)/(x2-x1)
		norm_slope = -1/reach_slope
		m_x, m_y = midpoint[0], midpoint[1]

		# Use the point slope formula to get the perpendicular bisector
		# y - m_y = norm_slope*(x - m_x)
		#print("y = %fx - %f\n" % (norm_slope, ((norm_slope*(-m_x)) + m_y)))

		# We know that the midpoint is a point on the reach, as well as on the normal
		# http://www.wolframalpha.com/input/?i=solve+for+x+500+%3D+(((x-g)%5E2)+%2B+((mx+%2B+b+-+h)%5E2))%5E(1%2F2)
		if(direction == 'out'):
			print("Getting coords for seaward reach normal")
			# solve for x
			b = ((norm_slope*(-m_x)) + m_y)
			x = (1/((norm_slope**2) + 1)) * (np.sqrt(-(b**2) - (2*b*m_x*norm_slope) + (2*b*m_y) \
			  - (m_x**2 * norm_slope**2) + (2*m_x*m_y*norm_slope) - (m_y**2) + (distance**2 * norm_slope**2) + distance**2) \
			  - (b*norm_slope) + m_x + (m_y*norm_slope))
		elif(direction == 'in'):
			print("Getting coords for landward reach normal")
			# solve for x
			b = ((norm_slope*(-m_x)) + m_y)
			x = (1/((norm_slope**2) + 1)) * (-np.sqrt(-(b**2) - (2*b*m_x*norm_slope) + (2*b*m_y) \
			  - (m_x**2 * norm_slope**2) + (2*m_x*m_y*norm_slope) - (m_y**2) + ((distance)**2 * norm_slope**2) + (distance)**2) \
			  - (b*norm_slope) + m_x + (m_y*norm_slope))

		# plug x into equation for normal line and solve for y
		y = (norm_slope * x) + b
		normal_coords = ((x, y), (m_x, m_y))
		return normal_coords

	# Generate a list of equally distanced points along each normal
	# Ensure that the points are generated from the sea inwards. ie. Clockwise
	def get_points_along_norm(self, norm, num_points, reach):
		print("\nGenerating %s equidistant points along the normal defined by: %s\n" % (str(num_points), str(norm)))
		points = []
		facing_shore = None
		# Points that define the normal
		x1, x2, y1, y2 = norm[0][0], norm[1][0], norm[0][1], norm[1][1]

		# Points that define the reach
		rx1, rx2, ry1, ry2 = reach[0][0], reach[1][0], reach[0][1], reach[1][1]

		delta_y = ry2 - ry1
		line_slope = (y2-y1)/(x2-x1)
		# Find the step interval in the x axis between the 2 different coordinates
		# Then solve b = y - m*x; y = m*x + b
		b = y1 - (line_slope*x1)
		x_interval = np.abs((x1-x2))/num_points
		#print("Step size: %f" % (float(np.sqrt((x1-x2)**2 + (y1-y2)**2)/num_points)))
		# Step out on the normal num_points times until the end of itself.
		# Add the first point to the list before calculating the next one
		# Ensure that we initialize points from the sea inwards
		if(delta_y < 0):
			x_temp = max(x1, x2)
			x_interval = -np.abs((x1-x2))/num_points
			facing_shore = False
		elif(delta_y >= 0):
			x_temp = min(x1, x2)
			x_interval = np.abs((x1-x2))/num_points
			facing_shore = True

		y_temp = (line_slope*x_temp) + b
		for i in range(0, num_points):
			points.append((x_temp, y_temp))
			y_temp = (line_slope*x_temp) + b
			x_temp += x_interval

		return points, facing_shore

	# Given a list of (x,y) coordinates along a normal, interpolate the z value at each point
	# This function involves the loading of a scipy interpolator and a bathymetric data file
	def get_z_values(self, points, file, interpolator):
		print("\nGetting the z values for the given list of coords...")
		np_x = np.array([coords[0] for coords in points])
		np_y = np.array([coords[1] for coords in points])
		# Interpolate the given x and y values to a list of z values
		z_values = interpolator(np_x, np_y)
		print('Done!\n')
		# Returns a list of z values associated with each point along a line
		return z_values.tolist()

	# Since the interpolator only returns z values, we need to zip in the distance to each z value on the normal as well
	# Stop adding points to the list as soon as there is a point that is above the upper elevation threshold
	# This function is not currently using the min_elevation variable, but it is there in case we come up with a use for it in the future
	def create_points_list(self, norm_length, points, step, max_elevation, min_elevation, facing_shore):
		points_distances = []
		for index, z_value in enumerate(points):
			if(z_value > max_elevation):
				break
			if(z_value <= max_elevation and facing_shore is False):
				points_distances.append((norm_length - ((index+1)*step), z_value))
			elif(z_value <= max_elevation and facing_shore is True):
				points_distances.append((((index+1)*step), z_value))
		return points_distances

	# Implementation of the Ramer-Douglas-Peuker Algorithm for curve simplification
	# More info here: https://github.com/sebleier/RDP/blob/master/__init__.py
	def ramer_douglas_peucker(self, points_distances, epsilon):
		# Reduces a series of points to a simplified version that loses detail,
		# but maintains the general shape of the series.
		dmax = 0.0
		index = 0
		results = []
		for i in range(1, len(points_distances) - 1):
			distance = self.perpendicular_distance(points_distances[i], points_distances[0], points_distances[-1])
			if distance > dmax:
				index = i
				dmax = distance
		if (dmax >= epsilon):
			# Recursive calls subset of the line
			results = self.ramer_douglas_peucker(points_distances[:index+1], epsilon)[:-1] \
			+ self.ramer_douglas_peucker(points_distances[index:], epsilon)
		else:
			try:
				results = [points_distances[0], points_distances[-1]]
			except IndexError as e:
				print(e)
				print("Empty Transect. No points in the desired elevation range could be found for this transect.")
		return results

	# Helper function used in perpendicular_distance
	def distance(self, a, b):
		return  np.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)

	# Helper function used in ramer_douglas_peucker
	# https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Line_defined_by_two_points
	def perpendicular_distance(self, point, start, end):
		if (start == end):
			return self.distance(point, start)
		else:
			n = np.abs(((end[1] - start[1])*point[0]) - ((end[0] - start[0])*point[1]) + end[0]*start[1] - end[1]*start[0])
			d = np.sqrt(((end[1] - start[1]) ** 2) + ((end[0] - start[0]) ** 2))
			return (n / d)

	# Given points along the normal, find their UTMx and UTMy coordinates.
	def get_UTM_coords(self, list, norm, midpoint, facing_shore, norm_length):
		print("\nFor every point along the normal, finding the UTM coords and writing them to self...\n")
		UTM_with_elevations =[]
		# slope = (y2-y1)/(x2-x1)
		# m_y = norm_slope*m_x + b
		norm_slope = (norm[1][1] - norm[0][1])/(norm[1][0] - norm[0][0])
		m_x, m_y = midpoint[0], midpoint[1]
		b = m_y - (norm_slope*m_x)
		# Starting at the norm start point, find the UTM_x and UTM_y values at each stored distance in the list
		x_start = norm[0][0]
		y_start = norm[0][1]
		for i in range (0, len(list)):
			curr_distance = list[i][0]
			# Error regarding the sqrt of a negative number: Happens for empty reaches
			try:
				UTM_x = (1/((norm_slope**2) + 1)) * (np.sqrt(-(b**2) - (2*b*x_start*norm_slope) + (2*b*y_start) \
				  - (x_start**2 * norm_slope**2) + (2*x_start*y_start*norm_slope) - (y_start**2) + (curr_distance**2 * norm_slope**2) + curr_distance**2) \
				  - (b*norm_slope) + x_start + (y_start*norm_slope))
				UTM_y = norm_slope*UTM_x + b
				if(facing_shore is False):
					UTM_with_elevations.append((UTM_x, UTM_y, list[i][1], norm_length - list[i][0]))
				elif(facing_shore is True):
					UTM_with_elevations.append((UTM_x, UTM_y, list[i][1], list[i][0]))
			except RuntimeWarning as e:
				print(e)
		return UTM_with_elevations

	# Create the interpolator in a function to allow more efficient running of the program
	def instantiate_interpolator(self, file):
		print("Creating the interpolator...")
		try:
			interpolator_array = file['LI']
			interpolator = interpolator_array[()]
			return interpolator
		except UnicodeError as e:
			print(e)
			print("Please run the program again.")
			sys.exit(0)

	# Helper function to visualize the a 2d list (the transects in this case) in a matplotlib window
	def plot_list(self, list, transect_number):
		fig = plt.figure()
		plt.title("Transect #%d" % (transect_number))
		# plot the transect
		plt.plot(*zip(*list))
		# Display the figure to the user
		plt.show()
		# Ensure that the figure object is closed in order to save memory space.
		plt.close()

	# Looks at the contents of the self.data dict and plots the values
	def plot_data(self, projection, lats, lons, m_x, m_y, norm_lons, norm_lats, transect_points):
		print("Plotting the reaches, transect points, and normals on a basemap instance...\n")
		fig = plt.figure()
		# Ensure that the aspect ratios of the axes are the same to avoid distortion in the output
		plt.gca().set_aspect('equal', adjustable='box')
		plt.title('Reaches and Transects')

		lon_temp = [point[0] for point in transect_points]
		lat_temp = [point[1] for point in transect_points]
		projection.scatter(lon_temp, lat_temp, color = 'r', s = 5)

		# draw coastlines, parallels and meridians
		if('RDN' in self.reaches_file):
			frame_lat = [49, 49.75]
			frame_lon = [-125,-123.5]
			projection.drawparallels(np.arange(frame_lat[0],frame_lat[1],0.25),labels = [1,0,0,0], color = 'black')
			projection.drawmeridians(np.arange(frame_lon[0],frame_lon[1],0.25),labels = [0,0,0,1], color = 'black')
		elif('RDT' in self.reaches_file):
			frame_lat = [49, 49.25]
			frame_lon = [-126, -125.75]
			projection.drawparallels(np.arange(frame_lat[0],frame_lat[1],0.05),labels = [1,0,0,0], color = 'black')
			projection.drawmeridians(np.arange(frame_lon[0],frame_lon[1],0.05),labels = [0,0,0,1], color = 'black')

		try:
			# Try to get better satellite based imagery from the ArcGIS servers
			projection.arcgisimage(service = 'ESRI_Imagery_World_2D', xpixels = 10000, verbose = True)
			projection.scatter(norm_lons, norm_lats, color = "white", s = 2)
		except urllib.error.HTTPError as e:
			print(e)
			print("Server unresponsive. Plotting without satellite imagery.")
			projection.scatter(norm_lons, norm_lats, color = "black", s = 2)
		projection.scatter(m_x, m_y, color = 'red')
		projection.drawcoastlines()
		projection.plot(lons, lats)
		# Plot the norms one by one
		normals = list(zip(norm_lons, norm_lats))
		reach_nodes = []
		midpoints = []
		for key, value in self.data.items():
			reach_nodes.append(value[0][0])
			reach_nodes.append(value[0][0])
			reach_nodes.append(value[0][1])
			reach_nodes.append(value[0][1])
			midpoints.append(value[2])

		print('\n')
		plt.plot(*zip(*reach_nodes))
		plt.scatter(*zip(*midpoints), color = 'r', s = 5)
		# join the points with lines
		for key, value in self.data.items():
			# Each normal is being plotted individually here
			plt.plot(*zip(value[1][0], value[1][1]))
		plt.show()

	# Instantiates a Basemap projection object
	# The projection can then be passed to other functions which can populate it with data
	# Uses the self.data dict as a parameter to establish the mean lat and longs
	def create_projection(self):
		proj_in = Proj(init = "epsg:32610")
		proj_out = Proj(init = "epsg:4326")
		# Initialize lists
		UTM_y, UTM_x, UTM_normx, UTM_normy, UTM_m_x, UTM_m_y, transect_points = ([] for i in range(7))
		# Iterate through self.data
		for key, value in self.data.items():
			# Skip the first entry in the dict which contains the column names
			if key == 0:
				continue
			# Tuple containg reach info
			UTM_x.append(value[0][0][0])
			UTM_x.append(value[0][1][0])
			UTM_y.append(value[0][0][1])
			UTM_y.append(value[0][1][1])
			UTM_normx.append(value[1][0][0])
			UTM_normx.append(value[1][1][0])
			UTM_normy.append(value[1][0][1])
			UTM_normy.append(value[1][1][1])
			UTM_m_x.append(value[2][0])
			UTM_m_y.append(value[2][1])

			# Convert the UTM transect point coords to equivalent long lat format
			for tup in value[4]:
				lon_transect_point, lat_transect_point = transform(proj_in, proj_out, tup[0], tup[1])
				transect_points.append((lon_transect_point, lat_transect_point))
		# Initialize lists
		lats, lons, norm_lats, norm_lons, m_x, m_y = ([] for i in range(6))
		# Convert the UTM coords to LAT/LONG
		for i in range (0, len(UTM_x)):
			# Transform returns coords in long lat format
			lon, lat = transform(proj_in, proj_out, UTM_x[i], UTM_y[i])
			nlon, nlat = transform(proj_in, proj_out, UTM_normx[i], UTM_normy[i])
			lats.append(lat)
			lons.append(lon)
			norm_lats.append(nlat)
			norm_lons.append(nlon)

		for i in range (0, len(UTM_m_x)):
			# Transform returns coords in long lat format
			lon_m_x, lat_m_y = transform(proj_in, proj_out, UTM_m_x[i], UTM_m_y[i])
			m_x.append(lon_m_x)
			m_y.append(lat_m_y)

		# Replace the mean calculations with the means of the reaches
		lat0 = np.mean(lats)
		lon0 = np.mean(lons)

		if('RDN' in self.reaches_file):
			print("\nCreating a projection centered on the Regional District of Nanaimo ...\n")
			# Create Basemap instance centered around the area of interest (Currently the Georgia Straight for RDN)
			projection = Basemap(projection='stere', lon_0 = lon0, lat_0 = lat0, lat_ts = lat0,
						llcrnrlon = -125, llcrnrlat = 49.15,
						urcrnrlon = -123.5, urcrnrlat = 49.75,
						rsphere = 6378137, resolution = 'f',
						area_thresh = 0.01, epsg = 4326)
		elif('RDT' in self.reaches_file):
			print("\nCreating a projection centered on the Regional District of Tofino ...\n")
			# Create Basemap instance centered around the area of interest (Currently Tofino)
			projection = Basemap(projection='stere', lon_0 = lon0, lat_0 = lat0, lat_ts = lat0,
						llcrnrlon = -126, llcrnrlat = 49,
						urcrnrlon = -125.75, urcrnrlat = 49.25,
						rsphere = 6378137, resolution = 'f',
						area_thresh = 0.01, epsg = 4326)
		print("\nDone")
		return (projection, lats, lons, m_x, m_y, norm_lons, norm_lats, transect_points)

	# Taking self.data as input, write the contents to a text file
	# UTM coords must be translated to Long Lat
	# Writes each line of the text file as: "LONG, LAT, DEPTH, DIST_ALONG_NORM"
	def export_transects(self, basename):
		print("\nExporting each reach to a seperate .txt file...")
		path = os.getcwd()
		relative_path = path + "/transect_files"
		proj_in = Proj(init = "epsg:32610")
		proj_out = Proj(init = "epsg:4326")

		for key, value in self.data.items():
			# Skip first dict entry which contains column names
			if key == 0:
				continue
			# create a new text file to store output for every transect
			filename = (relative_path + '/' + basename + '%02d.txt' % key)
			if(len(value[4]) > 0):
				with open(filename, 'w') as f:
					# For UTM coord in in transect
					for tup in value[4]:
						depth = tup[2]
						dist_along_norm = tup[3]
						# Convert the UTM_x and UTM_y values to equvalent Long/Lat
						lon, lat = transform(proj_in, proj_out, tup[0], tup[1])
						info = str(lon) + "," + str(lat) + "," + str(depth) + "," + str(dist_along_norm)
						f.write(info)
						f.write('\n')
			else:
				print("Reach #%d is empty. Not writing a .txt for it\n" % key)
				continue
		print("Done!\n")

	# Write the points of interest to .kml format to put them into google earth
	def export_transects_to_kml(self, basename):
		print("Exporting the generated transect .txt files to Google Earth readable .kml files...\n")
		# simplekml takes coords in long/lat format
		path = os.getcwd()
		relative_path = path + "/transect_files"
		for x in range (1, len(self.data)):
			point_coords = []
			kml = simplekml.Kml()
			file = (relative_path + '/' + basename + '%02d.txt' % (x))
			# Check if the file exists
			if os.path.isfile(file):
				with open(file, 'r') as f:
					# Check if the file was empty before making a line object
					if(os.stat(file).st_size > 0):
						for index, line in enumerate(f):
							line = line.split(',')
							line[-1] = line[-1].strip()
							depth = float(line[2])
							point = kml.newpoint(name = ("%.2fm" % depth), coords = [(float(line[0]),float(line[1]), depth)])
							point.style.iconstyle.icon.href = None
							point_coords.append((float(line[0]),float(line[1]), depth))
						line = kml.newlinestring(name = ("Transect %d" % (x)), coords = point_coords)
						line.altitudemode = simplekml.GxAltitudeMode.relativetoseafloor
						kml.save(relative_path + '/' + basename + '%02d.kml' % (x))
					else:
						print("Reach #%d is empty. Not writing a .kml for it\n" % x)
						continue
		print("Done!\n")


# The main function to run the class functions from. Currently contains a usage example
def main(argv):

	# Create new class object for Nanaimo -------------------------------------------------------------------------

	Nanaimo_Transects = Transect('Salish_DEM_20180615.npz', 'RDN_Reaches_UTM.txt')
	# Get the reach segment info from a file provided by the user
	Nanaimo_Reaches = Nanaimo_Transects.read_reach_file()
	# Inputs for Nanaimo:
	Nanaimo_Transects.get_transects(Nanaimo_Reaches, 1200, 0.25, 500, 5, -10, False)

	# Plot the self.data dict using matplotlib and basemap
	#projection, lats, lons, m_x, m_y, norm_lons, norm_lats, transect_points = Nanaimo_Transects.create_projection()
	#Nanaimo_Transects.plot_data(projection, lats, lons, m_x, m_y, norm_lons, norm_lats, transect_points)

	# Write the contents of the Nanaimo_Transects.reaches dict to text files for Nanaimo
	Nanaimo_Transects.export_transects('RDN_Transect')
	# Write text file transects to kml format for use on google earth for Nanaimo
	Nanaimo_Transects.export_transects_to_kml('RDN_Transect')

	# Create new class object for Tofino -------------------------------------------------------------------------

	Tofino_Transects = Transect('RDT_DEM.npz', 'RDT_Reaches_UTM.txt')
	# Get the reach segment info from a file provided by the user
	Tofino_Reaches = Tofino_Transects.read_reach_file()
	# Inputs for Tofino:
	Tofino_Transects.get_transects(Tofino_Reaches, 600, 0.25, 500, 5, -10, False)

	# Plot the self.data dict using matplotlib and basemap
	#projection, lats, lons, m_x, m_y, norm_lons, norm_lats, transect_points = Tofino_Transects.create_projection()
	#Tofino_Transects.plot_data(projection, lats, lons, m_x, m_y, norm_lons, norm_lats, transect_points)

	# Write the contents of the Tofino_Transects.reaches dict to text files for Tofino
	Tofino_Transects.export_transects('RDT_Transect')
	# Write text file transects to kml format for use on google earth for Tofino
	Tofino_Transects.export_transects_to_kml('RDT_Transect')

	sys.exit(0)

if __name__ == "__main__":
	main(sys.argv)
