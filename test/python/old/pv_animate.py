#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1546, 587]

# get camera animation track for the view
cameraAnimationCue1 = GetCameraTrack(view=renderView1)

# create keyframes for this animation track

# create a key frame
keyFrame3955 = CameraKeyFrame()
keyFrame3955.Position = [1057.4099999999994, 547.5129999999999, 2082.96]
keyFrame3955.FocalPoint = [399.531, 399.553, 4.48686]
keyFrame3955.ViewUp = [-0.07843177085312021, 0.9958545678032245, -0.04606665937852193]
keyFrame3955.ParallelScale = 565.551151870903
keyFrame3955.PositionPathPoints = [1057.41, 547.513, 2082.96, 2427.425710630599, 595.8037327437048, 794.3452184646411, 2294.043003614791, 498.60193436520905, -1079.838326788099, 756.1463783072436, 327.96905408889955, -2150.1473884479906, -1046.13032594249, 210.40535058450092, -1623.098477651225, -1776.6526612318517, 233.06775913549984, 110.57583729868887, -893.841160412001, 379.1552942627901, 1765.600122176601]
keyFrame3955.FocalPathPoints = [399.531, 399.553, 4.48686]
keyFrame3955.ClosedPositionPath = 1

# create a key frame
keyFrame3956 = CameraKeyFrame()
keyFrame3956.KeyTime = 1.0
keyFrame3956.Position = [1057.4099999999994, 547.5129999999999, 2082.96]
keyFrame3956.FocalPoint = [399.531, 399.553, 4.48686]
keyFrame3956.ViewUp = [-0.07843177085312021, 0.9958545678032245, -0.04606665937852193]
keyFrame3956.ParallelScale = 565.551151870903

# initialize the animation track
cameraAnimationCue1.Mode = 'Path-based'
cameraAnimationCue1.KeyFrames = [keyFrame3955, keyFrame3956]

# get animation scene
animationScene1 = GetAnimationScene()

animationScene1.Play()

# current camera placement for renderView1
renderView1.CameraPosition = [1057.4099999999994, 547.5129999999999, 2082.96]
renderView1.CameraFocalPoint = [399.531, 399.553, 4.48686]
renderView1.CameraViewUp = [-0.07843177085312021, 0.9958545678032245, -0.04606665937852193]
renderView1.CameraParallelScale = 565.551151870903

# current camera placement for renderView1
renderView1.CameraPosition = [1057.4099999999994, 547.5129999999999, 2082.96]
renderView1.CameraFocalPoint = [399.531, 399.553, 4.48686]
renderView1.CameraViewUp = [-0.07843177085312021, 0.9958545678032245, -0.04606665937852193]
renderView1.CameraParallelScale = 565.551151870903

# save animation images/movie
WriteAnimation('/home/grogan/test.ogv', Magnification=1, FrameRate=15.0, Compression=True)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [1057.4099999999994, 547.5129999999999, 2082.96]
renderView1.CameraFocalPoint = [399.531, 399.553, 4.48686]
renderView1.CameraViewUp = [-0.07843177085312021, 0.9958545678032245, -0.04606665937852193]
renderView1.CameraParallelScale = 565.551151870903

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).