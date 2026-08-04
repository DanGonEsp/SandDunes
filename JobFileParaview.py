# trace generated using paraview version 6.0.1
#import paraview
#paraview.compatibility.major = 6
#paraview.compatibility.minor = 0
import os
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# === 1. Input file ===
input_file = 'master_Solution.pvd'  # change this
input_folder = '/Users/gonzald/Desktop/Stabilization/TableTest_euler/Table3Debug_0000-tri-lev3-Stokes-RelVel-NoForce-AveDiff-NoSlip_copy/master_Solution.pvd'

# === 2. Get input directory ===
input_dir = os.path.dirname(input_folder)

# === 3. Create output folder next to input file ===
output_folder = os.path.join(input_dir, "paraview_images")
os.makedirs(output_folder, exist_ok=True)

# create a new 'PVD Reader'
master_Solutionpvd = PVDReader(registrationName=input_file, FileName=input_folder)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# Adjust camera

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# show data in view
master_SolutionpvdDisplay = Show(master_Solutionpvd, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
master_SolutionpvdDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera(False, 0.9)

#changing interaction mode based on data extents
renderView1.Set(
    InteractionMode='2D',
    CameraPosition=[39.5, 7.092800140380859, 264.65000000000003],
    CameraFocalPoint=[39.5, 7.092800140380859, 0.0],
)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(master_SolutionpvdDisplay, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
master_SolutionpvdDisplay.SetScalarBarVisibility(renderView1, True)
# Adjust camera

# current camera placement for renderView1
renderView1.Set(
    InteractionMode='2D',
    CameraPosition=[39.5, 7.092800140380859, 264.65000000000003],
    CameraFocalPoint=[39.5, 7.092800140380859, 0.0],
    CameraParallelScale=40.13175567840743,
)

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')
vtkBlockColorsLUT.Set(
    InterpretValuesAsCategories=1,
    AnnotationsInitialized=1,
    Annotations=['0', '0', '1', '1', '2', '2', '3', '3', '4', '4', '5', '5', '6', '6', '7', '7', '8', '8', '9', '9', '10', '10', '11', '11'],
    ActiveAnnotatedValues=['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11'],
    IndexedColors=[1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.63, 0.63, 1.0, 0.67, 0.5, 0.33, 1.0, 0.5, 0.75, 0.53, 0.35, 0.7, 1.0, 0.75, 0.5],
)

# get opacity transfer function/opacity map for 'vtkBlockColors'
vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

# get 2D transfer function for 'vtkBlockColors'
vtkBlockColorsTF2D = GetTransferFunction2D('vtkBlockColors')
# Adjust camera

# current camera placement for renderView1
renderView1.Set(
    InteractionMode='2D',
    CameraPosition=[39.5, 7.092800140380859, 264.65000000000003],
    CameraFocalPoint=[39.5, 7.092800140380859, 0.0],
    CameraParallelScale=40.13175567840743,
)

# set scalar coloring
ColorBy(master_SolutionpvdDisplay, ('POINTS', 'c'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(vtkBlockColorsLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
master_SolutionpvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
master_SolutionpvdDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'c'
cLUT = GetColorTransferFunction('c')
cLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'c'
cPWF = GetOpacityTransferFunction('c')
cPWF.ScalarRangeInitialized = 1

# get 2D transfer function for 'c'
cTF2D = GetTransferFunction2D('c')


# get color legend/bar for cLUT in view renderView1
cLUTColorBar = GetScalarBar(cLUT, renderView1)
cLUTColorBar.Set(
    Title='c',
    ComponentTitle='',
)




# change scalar bar placement
cLUTColorBar.Set(
    Orientation='Horizontal',
    WindowLocation='Any Location',
    Position=[0.37092975206611567, 0.84657060518732],
    ScalarBarLength=0.32999999999999996,
)

# Adjust camera

# current camera placement for renderView1
renderView1.Set(
    InteractionMode='2D',
    CameraPosition=[18.45437148293592, 5.667747943838653, 264.65000000000003],
    CameraFocalPoint=[18.45437148293592, 5.667747943838653, 0.0],
    CameraParallelScale=5.965329763017001,
)

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(2178, 1388)

# Adjust camera

# current camera placement for renderView1
renderView1.Set(
    InteractionMode='2D',
    CameraPosition=[18.45437148293592, 5.667747943838653, 264.65000000000003],
    CameraFocalPoint=[18.45437148293592, 5.667747943838653, 0.0],
    CameraParallelScale=5.965329763017001,
)

# === 8. Save animation as JPG sequence ===
output_pattern = os.path.join(output_folder, "frame.png")

# save animation
SaveAnimation(filename=output_pattern, viewOrLayout=renderView1, location=16, ImageResolution=[2178, 1388],
    TransparentBackground=1,
    FrameWindow=[0, 116])




#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(2178, 1388)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.Set(
    InteractionMode='2D',
    CameraPosition=[18.45437148293592, 5.667747943838653, 264.65000000000003],
    CameraFocalPoint=[18.45437148293592, 5.667747943838653, 0.0],
    CameraParallelScale=5.965329763017001,
)


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://www.paraview.org/paraview-docs/nightly/python/
##--------------------------------------------


print("Saved image sequence to:", output_folder)


# Borrar
#renderView1 = GetActiveViewOrCreate('RenderView')

#renderView1.Set(
    #InteractionMode='2D',
    #CameraPosition=[18.45437148293592, 5.667747943838653, 264.65000000000003],
    #CameraFocalPoint=[18.45437148293592, 5.667747943838653, 0.0],
    #CameraParallelScale=5.965329763017001,
#)
#Render()
