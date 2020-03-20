# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview



import numpy as np
from paraview.simple import *
from os.path import isfile, join
from os import listdir

data_dir = 'D:/Jones_Lab/hnn_params/input_time_sweep/input_time_data/'
save_dir = 'D:/Jones_Lab/hnn_params/input_time_sweep/input_time_skeleton/'

file_list = [f for f in listdir(data_dir) if isfile(join(data_dir, f))]


for csd_file in file_list:
    csd_file_name = csd_file.strip('.csv')
    csd_file_path = data_dir + csd_file

    #Create path + file name to save skeleton node and arc data
    csd_node_string = save_dir + csd_file_name + '_nodes.csv'
    csd_arc_string = save_dir + csd_file_name + '_arcs.csv'


    #-------- Paraview Scripting Code (start)--------------

    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'CSV Reader'
    surface_csv = CSVReader(FileName=[csd_file_path])

    # Properties modified on surface_csv
    surface_csv.HaveHeaders = 0

    # Create a new 'SpreadSheet View'
    spreadSheetView1 = CreateView('SpreadSheetView')
    spreadSheetView1.ColumnToSort = ''
    spreadSheetView1.BlockSize = 1024
    # uncomment following to set a specific view size
    # spreadSheetView1.ViewSize = [400, 400]

    # get layout
    layout1 = GetLayout()

    # place view in the layout
    layout1.AssignView(2, spreadSheetView1)

    # show data in view
    test_sweept_evprox_1_0_t_evdist_1_0csvDisplay = Show(surface_csv, spreadSheetView1)

    # update the view to ensure updated data information
    spreadSheetView1.Update()

    # create a new 'Table To Points'
    tableToPoints1 = TableToPoints(Input=surface_csv)

    # Properties modified on tableToPoints1
    tableToPoints1.XColumn = 'Field 0'
    tableToPoints1.YColumn = 'Field 1'
    tableToPoints1.ZColumn = 'Field 2'

    # show data in view
    tableToPoints1Display = Show(tableToPoints1, spreadSheetView1)

    # hide data in view
    Hide(surface_csv, spreadSheetView1)

    # update the view to ensure updated data information
    spreadSheetView1.Update()


    # Fetch points from server and store z column bounds for elevation
    data_fetch = servermanager.Fetch(tableToPoints1) # Look into limiting fetch to just the array bounds
    point_bounds = data_fetch.GetPoints().GetBounds()
    elevation_min, elevation_max = point_bounds[4], point_bounds[5]

    # create a new 'Delaunay 2D'
    delaunay2D1 = Delaunay2D(Input=tableToPoints1)

    # show data in view
    delaunay2D1Display = Show(delaunay2D1, spreadSheetView1)

    # hide data in view
    Hide(tableToPoints1, spreadSheetView1)

    # update the view to ensure updated data information
    spreadSheetView1.Update()

    # create a new 'Tetrahedralize'
    tetrahedralize1 = Tetrahedralize(Input=delaunay2D1)

    # show data in view
    tetrahedralize1Display = Show(tetrahedralize1, spreadSheetView1)

    # hide data in view
    Hide(delaunay2D1, spreadSheetView1)

    # update the view to ensure updated data information
    spreadSheetView1.Update()

    # create a new 'Elevation'
    elevation1 = Elevation(Input=tetrahedralize1)

    # Properties modified on elevation1
    elevation1.LowPoint = [0, 0, elevation_min]
    elevation1.HighPoint = [0, 0, elevation_max]

    # show data in view
    elevation1Display = Show(elevation1, spreadSheetView1)

    # hide data in view
    Hide(tetrahedralize1, spreadSheetView1)

    # update the view to ensure updated data information
    spreadSheetView1.Update()

    # create a new 'TTK PersistenceDiagram'
    tTKPersistenceDiagram1 = TTKPersistenceDiagram(Input=elevation1)

    # show data in view
    tTKPersistenceDiagram1Display = Show(tTKPersistenceDiagram1, spreadSheetView1)

    # hide data in view
    Hide(elevation1, spreadSheetView1)

    # update the view to ensure updated data information
    spreadSheetView1.Update()

    # create a new 'TTK TopologicalSimplification'
    tTKTopologicalSimplification1 = TTKTopologicalSimplification(Domain=elevation1,
        Constraints=tTKPersistenceDiagram1)

    # Properties modified on tTKTopologicalSimplification1
    tTKTopologicalSimplification1.OutputOffsetScalarField = ''

    # show data in view
    tTKTopologicalSimplification1Display = Show(tTKTopologicalSimplification1, spreadSheetView1)

    # hide data in view
    Hide(tTKPersistenceDiagram1, spreadSheetView1)

    # hide data in view
    Hide(elevation1, spreadSheetView1)

    # update the view to ensure updated data information
    spreadSheetView1.Update()

    # create a new 'TTK Merge and Contour Tree (FTM)'
    tTKMergeandContourTreeFTM1 = TTKMergeandContourTreeFTM(Input=tTKTopologicalSimplification1)

    # show data in view
    tTKMergeandContourTreeFTM1Display = Show(tTKMergeandContourTreeFTM1, spreadSheetView1)

    # hide data in view
    Hide(tTKTopologicalSimplification1, spreadSheetView1)

    # show data in view
    tTKMergeandContourTreeFTM1Display_1 = Show(OutputPort(tTKMergeandContourTreeFTM1, 1), spreadSheetView1)

    # hide data in view
    Hide(tTKTopologicalSimplification1, spreadSheetView1)

    # show data in view
    tTKMergeandContourTreeFTM1Display_2 = Show(OutputPort(tTKMergeandContourTreeFTM1, 2), spreadSheetView1)

    # hide data in view
    Hide(tTKTopologicalSimplification1, spreadSheetView1)

    # update the view to ensure updated data information
    spreadSheetView1.Update()

    # show data in view
    tTKMergeandContourTreeFTM1Display_1 = Show(OutputPort(tTKMergeandContourTreeFTM1, 0), spreadSheetView1)

    # export view
    ExportView(csd_node_string, view=spreadSheetView1)

    # show data in view
    tTKMergeandContourTreeFTM1Display_1 = Show(OutputPort(tTKMergeandContourTreeFTM1, 1), spreadSheetView1)

    # Properties modified on spreadSheetView1
    spreadSheetView1.FieldAssociation = 'Cell Data'

    # export view
    ExportView(csd_arc_string, view=spreadSheetView1)

    ResetSession()

    #-------- Paraview Scripting Code (End)--------------
