README for shoreline_transects.py
Author: Henri De Boever
2018/07/23

A shoreline transect generator written by Henri De Boever for use by Cascadia Coast Research.

---------------------- Table of Contents ----------------------

    1. Core Dependencies
    2. Getting the Reach information into the program
      2.1. Converted reach file example
    3. Instantiating a new Transect Class Object
    4. Getting the reaches from the reach file into the object
    5. Defining Transects
      5.1. Pseudo code for get_transects
    6. Plotting Reach and Transect Data
    7. Exporting Reach and Transect Data
      7.1. Export usage example
    8. Function by Function Breakdown
    9. Usage Example
    10. Potential Improvements/Current Issues

This module serves to generate perpendicular transects given an input list of reaches that define an area of interest.

1. Core Dependencies --------------------------------------------------

  This code was run and tested on Anaconda for Python 3.6.5 (conda 4.5.2)

  from scipy.interpolate import LinearNDInterpolator   # Used to interpolate depths and altitudes from the bathymetric data
  from mpl_toolkits.Basemap import Basemap             # Used for basic plotting
  from collections import OrderedDict                  # Ordered Dict to store data etc.
  from pyproj import Proj, transform                   # Used to convert coordinates in different projections (in this case EPSG:32610 to EPSG:4326)
  import matplotlib.pyplot as plt                      # Used to generate transect and Basemap plots
  import numpy as np, simplekml                        # Used for the mathematical operations on lines, used to generate .kml files
  import sys, os, urllib                               # Used for error handling, file and directory structure manipulation, HTTP error checking

2. Getting the Reach information into the program ---------------------------------------------------

  At the higher level, the tool requires a .kml or .kmz file containing the reach endpoints. The typical way to generate
  reach files is to use the path making tool on Google earth. Exporting the file from Google Earth will give the user a .kml file.

  The next step is to convert the raw .kml file to .txt format. For the purpose, this website performs the trick: http://zonums.com/online/kml2x/
  Convert the .kml file to a UTM meters projection. In the current case, since we were working in the regional districts of Tofino and Nanaimo,
  the corresponding UTM Zone to select is Zone 10.

  Usage:

  I. Generate a path on Google Earth to define reaches.
  II. Export the .kml file and save it to the transect/reach_kmz_files directory: e.g. Example_Reaches.kmz
  III. Go to http://zonums.com/online/kml2x/, export a .txt format file in the corresponding UTM meters projection needed: e.g. Example_Reaches_UTM.txt
  IV. Ensure that the text file does not have an extra new line character at the end of it!

      2.1 An example of a converted reach file should look like this:
      ----------------------------------------------------------
      KML2x Export:
      Units: UTM(meters), Zone 10
      Datum: WGS84

      id XCoo YCoo Elevation
      RDT_Reaches vertex: 6
      288981.54464152345 5441936.811227419 0
      289093.15019439533 5442279.848374281 0
      290069.1394485142 5442241.687932079 0
      290056.7059098609 5443199.614280624 0
      289725.7736959002 5443321.594058154 0
      289655.8516967584 5443591.262702237 0
      ----------------------------------------------------------

  The above file contains the clockwise definition of 6 points, and 5 reach segments.

  IMPORTANT: To ensure correct functionality of this program, the reaches must be traced on Google Earth in a clockwise fashion.
  ie. if you were walking along each reach in the order they were drawn, the shore will always be to your right, and the ocean to your left.

3. Instantiating a new Transect Class Object --------------------------------------------------

  When creating a instance of a Transect object, the user must pass in two arguments:

  The first argument is the bathymetric data file with a .nzp file extension. The second parameter is the name of the .txt file generated in the previous
  step which contains the UTM reaches data. Ensure that the .npz file provided overlaps with the region in which the reaches are declared.

  The following example is creating a Transect instance for the Regional District of Nanaimo:

  # Creating new class object
	transect_object = Transect('Salish_DEM_20180615.npz', 'RDN_Reaches_UTM.txt')

4. Getting the reaches from the reach file into the object --------------------------------------------------

  In the __init__ function of the Transects class, we declare 2 file references, as well as an Ordered Dictionary to contain
  the attributes of the reaches and the transects. The reach attributes are known from the reaches file, and since all the attributes of the transects
  must be derived from the reach info, the next step is to read the reaches file. read_reach_file ignores the header of the text file by
  skipping the first 6 rows before reading the points.

  # Get the reach segment info from a file provided by the user
	reach_coords = transect_object.read_reach_file()

  The read_reach_file(self) function only takes self in as a parameter. It proceeds to read the UTM X, Y coordinates and appends them to a list.
  The list contains tuples defining a line in 2 dimensional space. Since each row in the reach file represents a point, a reach can be drawn
  between 2 points and stored into the list as a reach.

  A reach in the list takes the form ((x1, y1),(x2, y2)). The list is composed of reaches in the form described just previously. As such, a reach file must
  contain an even number of rows (points) and therefore the number of reaches defined in the file is equal to (number of rows - 1).

  The list generated by read_reach_file is then passed onto the get_transects function, where the transects definitions and calculations take place.

5. Defining Transects --------------------------------------------------

  We are able to establish perpendicular transects at the midpoint of each reach since we know that the midpoint is a point on both lines,
  we can determine the slope intercept formula for the reach, and since we can calculate the slope of each line. Given these known values,
  we can calculate the endpoints of the perpendicular normal at the midpoint of each reach out to a desired distance on either side of the reach.

  get_transects(reaches list, length of the normal, epsilon factor for RDP, number of points to generate along each normal, max elevation, min min_elevation, plot flag)

  Usage example

  # Populate the self.data dict
  # Inputs for Tofino: reach_coords, 600, 0.25, 500, 5, -10, False
  # Inputs for Nanaimo: reach_coords, 1200, 0.25, 500, 5, -10, False
  transect_object.get_transects(reach_coords, 1200, 0.25, 500, 5, -10, False)

  -  reach_coords is the list returned by the read_reach_file() function.
  -  1200 is the desired distance in metres desired for the maximum transect length.
  -  0.25 is the smoothing factor that is taken by the RDP algorithm for line simplification(retain elevation differences of 0.25 between significant points)
  -  500 is the number equidistant points along each transect normal that will be interpolated
  -  5 is the maximum elevation upper bound in metres along a transect.
  -  -10 is the minimum elevation lower bound for a transect (not actually used, but passed in case we find use for it in the future)
  -  False is the state of the plotting flag. If turned on to True, the user will be able to see the shape of each transect as it is
     calculated while the program is running. The plotting functionality is mostly implemented using matplotlib

  get_transects makes usage of a Linear Interpolator from scipy.interpolate to determine the depths and elevations given a UTM coordinate.

  5.1 Pseudo code for get_transects --------------------------------------------------

  instantiate the interpolator by loading the .npz file data
  write into the the first index of the self.data ordered Dictionary the names of the columns
  # Append the column names to the dict
  self.data[0] = ("reach", "normal", "midpoint_coord", "reach_length", "transect_UTM_coords")

  for each reach in reaches list:

    # compute the reach length, reach midpoint, and normal endpoints

    # Generate a list of points along the normal
    # For each point on the normal, get the elevation at that point using the interpolator and the .npz file
    # Combine each elevation value with its distance along the normal starting from the seaward terminus of the normal
      while the elevation is less than the upper elevation cutoff.

    # Perform the Ramer Douglas Peucker algorithm on the list of points that define each transect. (This preserves the general shape of the
    transect, but greatly reduces the number of points that need to be stored in memory.) In the case that RDP fails, it is because the reach and
    the normal for the current reach was defined in an area where the points and their elevations did not meet the criteria (too high, too low,
    too far away, etc.). In that case, the transect for that reach will be empty.

    # If the plot flag is on, plot the transect for the user

    # For every point that has been kept after RDP along the transect, get the associated UTM_X and UTM_Y coordinates.
    # Finally, store reach, normal, midpoint_coord, reach_length, transect_UTM_coords into the index of the OrderedDict self.data

6. Plotting Reach and Transect Data --------------------------------------------------

  There are 2 ways to interact with the ouput of shoreline_transects.py so far. This first way is to use the built in Basemap plotting functions that
  are written. The second way involves exporting the transects from self.data to a .kml file format. This allows the user to see the output of the
  program in Google Earth! A much more pleasant, and interactive experience! Basemap was integrated in the case that a future user couldn't get Google
  Earth working on their machine, but Basemap offers little compared to Goolge Earth. Plotting in Basemap can also induce projection based distortions
  and aspect ratio problems, which I have tried my best to address in the current versioning on the code. Still, perpendiculars will not always look
  orthogonal on a Basemap instance. As the Basemap imagery is polled from server, calls to these functions may take up to a couple minutes to run.
  The area of interest is determined on the names of the input reaches text files at the moment.
  Basemap projecting should be used as a last resort when all else fails.

  In either case, this is how the natively coded Basemap functions are invoked:

  # Plot the self.data dict using matplotlib and Basemap
	projection, lats, lons, m_x, m_y, norm_lons, norm_lats, transect_points = transect_object.create_projection()
	transect_object.plot_data(projection, lats, lons, m_x, m_y, norm_lons, norm_lats, transect_points)

  The first call returns attributes about the projection which are then plotted in the second call. A mix of scattering and plotting
  was used to get the desired effects on the Basemap layers.

7. Exporting Reach and Transect Data --------------------------------------------------

  The exportation of the data into .kmz format takes place in 2 steps. Firstly, the self.data dictionary's transect points are
  written to text files in transect/transect_files. Then, each file in transect/transect_files in read and a .kml copy is made
  in the same directory. For this, I used the python module simplekml.

    7.1 Export usage ---------------

    # Write the contents of the transect_object.reaches dict to text files
  	transect_object.export_transects()

  	# Write text file transects to kml format for use on google earth
  	transect_object.export_transects_to_kml()

  After the above calls, you should have .kml files that you can know import into Google earth. The elevation profiles of each transect
  can also be viewed in Google Earth since the elevations have been initialized relative to sea level with simplekml.
  If the files already exist in the target directory, they will be overwritten.

8. Function by Function Breakdown --------------------------------------------------------

  def __init__(self, bathymetric_data_file, reaches_file):

    Input parameters: relevant .npz file, and relevant .txt reach file

    This function gets called every time an instance of the class is declared. 3 instance variables are declared in this function.
    - self.data is an ordered dict which will contain the information about the reach, its midpoint, and the relevant normal extending to a desired
      distance perpendicular to the midpoint. The last item in the dict in a list of the points of interest along each normal with their depths.
    - self.npz_file is a reference to the relevant bathymetric data file for the current instance.
    - self.reaches_file is a reference to the relevant converted .kml text file for the current instance

    After declaring these variables, the class instance checks for the existence of 2 sub directories: transect_files, and reach_files. These sub
    directories are required to store and use the relevant data for the running of this program. If the reach_files folder is not created, the
    program will create it, and then prompt the user to put reach .txt data into it before exiting;
    the program cannot run without a reach file to read from!


  def read_reach_file(self):

    Input parameters: self
    This function is an instance function which reads the reaches text file from self.reaches_file.
    Before the function reads from the file, it checks for the file's existence in the reach_files directory. If that is false, then the user is
    told that a file of that name wasn't was not found before the program exits.

    If the file is found, then the function opens it in read mode ('r').  We skip the first 6 lines where the header of the file is located, then
    we read the coordinates line by line. Since a reach is defined by 2 coordinates, we zip consecutive coordinates together before appending them to a
    list.

    The function returns the list of reaches. The list contains tuples of 2 nodes each, which each define a reach.


  def get_transects(self, reach_coords, norm_length, epsilon, num_points, upper_elevation_limit, lower_elevation_limit, transect_plot_flag):

    Input parameters: Desired length of normal, epsilon scaling facotr for RDP, number of points to generate along norm, upper elevation limit
                      of a transect, lower elevation limit of a transect, the optional transect plotting flag.

    The bulk of the logic for this program is encapsulated within this function. Most of the defining factors of a transect are passed from the
    user to the program as inputs to this function. Since we are able to derive the transect normals only given the reach coordinates, this function
    is centred around populating the self.data ordered dictionary for each reach/transect one by one by incrementally calculating desired values
    before appending the relevant data to the dict.

    The function declares a local instance of the interpolator by invoking the reference held at self.npz_file. Then the "column names" of the
    data are appended to self.data at the 0th index.

    Then for each reach element in the reach coords list, we compute the reach length, the reach midpoint, and the end points of the normal out the desired
    distance. Then equidistant points are generated along the normal from the seaward endpoint towards land. Then the elevation at each one of these points
    is calculated using the interpolator which takes UTM_x and UTM_y locations on the .npz file to return the associated elevation value.

    Then in order for us to be able to establish a transect, we must determine the horizontal distance of each point of interest along the normal
    in relation to the seaward normal endpoint. The next function zips the z values and the step out distance from the seaward normal endpoint
    into a list. The resulting list is reduced in complexity by the RDP algorithm. The function then optionally plots the transect to the user if
    the plot flag was set. It is again important to state that this program will has defined functionality only for clockwise defined reaches.
    Otherwise there is no way to tell between a seaward endpoint and a landward endpoint of a transect.

    The next step is to determine the UTM_x and UTM_y coordinates given step out distance along each normal of each point along the line.

    Once the function has completed the aforementioned steps, the dictionary self.data at the current (index+1) of the loop is updated to contain
    the reach coords, the normal coords, the midpoint coord, the reach length, and the list of the coordinates and depths of the points along each transect

  def get_reach_length(self, reach):

    Input parameters: a tuple containing the 2 coordinates which define a reach

    This function is called within get_transects to get the length of the reach. Given the 4 values that define a line ((x1, y1),(x2, y2))
    the length of the reach is determined using the formula for Euclidean distance. Uses np.sqrt.

    Return values: The length of the reach.

  def get_reach_midpoint(self, reach):

    Input parameters: a tuple containing the 2 coordinates which define a reach

    This function is called within get_transects to get the midpoint coordinates of the reach. Given the 4 values that define a line ((x1, y1),(x2, y2))
    the midpoint of the reach is determined using the midpoint formula.

    Return values: A tuple containing (midpoint_x, midpoint_y)

  def get_reach_normal(self, reach, midpoint, distance, direction):

    Input parameters: a tuple containing the 2 coordinates which define a reach, tuple containing (midpoint_x, midpoint_y), desired length of normal/2,
    direction which simply interchanges the formula used to calculate the coordinates of a point a certain distance away on either side of the reach.

    This function is called by get_transects to get the normal of each reach.
    Computes the reach slope and the slope of the normal. Given these values and a point which is common to both lines, the equation of
    the transect normal can be computed.

    Return values: A tuple which contains the endpoint of a transect normal, either the seaward one, or the landward one, depending on the direction.
    Clarification on the math can be found on the wolfram alpha link below
    Nota Benne: http://www.wolframalpha.com/input/?i=solve+for+x+500+%3D+(((x-g)%5E2)+%2B+((mx+%2B+b+-+h)%5E2))%5E(1%2F2)

  def get_points_along_norm(self, norm, num_points, reach):

    Input parameters: A tuple containing the 2 endpoints of the transect normal, the desired number of points along the normal, a tuple containing
    the reach coordinates.

    This function is called by get_transects to get a list of coordinates at equal intervals along each transect normal.
    The reach is used to tell the orientation with respect to the shore of the transect. Number of points is used to
    calculate the distance between each point along the line. The norm is used to solve for the y value of the line for a
    given x a at the number of steps away from the seaward endpoint of the normal.

    Return values: A list of UTM_x and UTM_y tuples of coordinates along the line, a boolean to denote the orientation of the transect with respect to the
    shore based on a clockwise definition of reaches.

  def get_z_values(self, points, file, interpolator):

    Input parameters: list of points generated by get_points_along_norm, reference to the .npz file in self.npz_file, the interpolator declared earlier
    in get_transects.

    This function is called in get_transects to get a list of corresponding elevation values at each point in the points list. Uses interpolator from
    scipy.interpolate LinearNDInterpolator.

    Return values: 1 dimensional list of z values


  def create_points_list(self, norm_length, points, step, max_elevation, min_elevation, increasing):

    Input parameters: length of the normal, points from get_points_along_norm, step size, min and max elevations, boolean for shore orientation

    This function is called in get_transects to zip the elevation values foundin get_z_values with the points generated in get_points_along_norm.
    Depending on the orientation of the shore, the distances are zipped either from the landward end of the normal, or from the seaward end of the normal.

    Return values: A list containing the elevation of each point along the transect, and it's distance from the seaward endpoint of the normal.

  def ramer_douglas_peucker(self, points_distances, epsilon):

    Input Parameters: the 2 dimensional list containing tuples with distances and elevations, the resolution scaling factor epsilon.
    This function is an implementation of the Ramer-Douglas-Peuker Algorithm for 2 dimensional curve simplification. It is called within get_transects
    in order to simplify the number of points in the final list to be stored. The call of RDP conserves general curve shape while removing points
    that don't meet the epsilon criteria fro being kept in the curve. You lose a bit of curve resolution. but the advantage is that the number of points
    is diminished greatly.

    Return vluaes: A shortened list of points and elevation along each transect.
    More info here: https://github.com/sebleier/RDP/blob/master/__init__.py

  def distance(self, a, b):

    Input parameters: points a and b

    A helper function called in ramer_douglas_peucker to determine the 2D distance between 2 points.

    Returns the distance between a and b.

  def perpendicular_distance(self, point, start, end):

    Input parameters: points a and b

    A helper function called in ramer_douglas_peucker to point to line distance.

    Returns the distance between a point, and the line defined by start and end.

  def get_UTM_coords(self, list, norm, midpoint, facing_shore, norm_length):

    Input parameters: list of step distances and elevations along the transect, the 2 coordinates that define the transect normal
    the midpoint coordinate,the orientation boolean, the length of the normal

    This function is called in get_transects. Using the know start and end points, the step size distances of the points along the normal and
    their depths, we can obtain a list of the UTM_x and UTM y coordinates associated with the depths at those locations.

    Return values: A list continaing tuples of (UTM_x, UTM_y, elevation in metres, step out distance of that point from the seaward normal endpoint)

  def instantiate_interpolator(self, file):

    Input parameters: reference to self.npz_file
    Return values: Returns an interpolator object associated with the bathymetric data file provided.

  def plot_list(self, list, transect_number):

    Interacts with matplotlib to plot a transect to the user. It is called by get_transects when the plot flag is turned on. Displays
    a cross section of the transect. Can be good for debugging purposes, or simply curiosity.

    Return values : The input parameters for plot_data


  def plot_data(self, projection, lats, lons, m_x, m_y, norm_lons, norm_lats, transect_points):

    Interacts with matplotlib to plot a a geographical projection of the reaches, and the points along the transects.
    It can be used as a bare-bones alternative to google earth when the former is unavailable.


  def create_projection(self):

    Creates a projection onto the appropriate longitudes and latitudes for Nanaimo or for Tofino and passes it to plot_data.
    Currently the projection areas are hard coded to work only for Tofino or Nanaimo, based on the name of the reach file
    provided to the program.


  def export_transects(self, basename):

    Exports the self.data dictionary to a text file containing the coordinates and the depth of each point of interest of a
    transect on each line. The basename is derived from the reach file name. Outputs .txt files for each reach in the transect_files
    directory. One file per transect. If the file does not exist it will be created. If the file exists, it will be overwritten.


  def export_transects_to_kml(self, basename):

    Reads the transect files in transect_files and generates google earth readable .kml files. Deposits them into the transect_files folder.
    One file per transect. If the file does not exist it will be created. If the file exists, it will be overwritten. Another stipulation
    to writing the elevation data to kml is to ensure that the values for the layer are specified relative to the seafloor. This feature
    allows the elevation profiles to also be viewed within google earth on the "Show Elevation Profile" option for each .kml transect path.

9. Usage example ------------------------------------------------------------

    The following example shows how the class can be used and interacted with. It contains 2 examples of usage for RDN and RDT.
    This assumes that you have the python script and the files in question in the same directory as shoreline_transects.py. I still
    need to figure out a way to install it globally.


    from shoreline_transects import Transect
    import sys

    # The main function to run the class functions from
    def main(argv):

    	# Create new class object for Nanaimo

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

    	# Create new class object for Tofino

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

10. Potential Improvements/Current Issues --------------------------------------------------

  - Cannot adequately handle a transect that straddles an island, thereby creating a kind of bowl shape. Either we search the transect for
    the entire distances out of the midpoint of the reach, or we simply define the transects to be shorter to avoid that problem.
    Also redrawing reaches in an optimal way to avoid those cases as much as possible.

  - Relies on reaches to drawn well to function well.

  - Lack of adequate warning when a transect is empty. This can cause issues when reading the text files in to other programs.

  - It wouldn't be a bad idea to re-write the read_reach_file function to accept formats different from the ones provided
    by: http://zonums.com/online/kml2x/

  - THIS HAS SINCED BEEN FIXED: At times, when the transect is empty, after the export to .kml format,
    you may get locations placed at 0 longitude and 0 latitude. Simply manually edit the text file
    and/or the .kml file, or ignore the bug.
