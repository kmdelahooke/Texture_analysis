#Run on copy of output from "grid_processing_hi.py"

import geomagic.app.v3
from geomagic.app.v3.imports import *

#delete mesh
models = geoapp.getModels()
for model in models:
	if model.name == "Mesh 1":
		geoapp.setActiveModel(model)

activeModel = geoapp.getActiveModel()
mesh = geoapp.getMesh(activeModel)
geoapp.deleteModel(mesh)

#downsample
for model in models:
   if model.name != "World":
       geoapp.setActiveModel(model)

activeModel = geoapp.getActiveModel()
pts = geoapp.getPoints(activeModel)

sel = PointSelection(pts)
sel.selectAll()

mod = SamplePointsRandomly()
target = int(pts.numPoints*0.005) # down sample to 0.5%
print(target)
mod.targetNumPoints = target
mod.selection = sel
mod.run()

geo.clear_all()

#clean
geo.repair_point_normals(0)

#mesh
activeModel = geoapp.getActiveModel()
pts = geoapp.getPoints(activeModel)
wrapper = WrapPoints()
wrapper.points = pts
wrapper.run()
mesh = wrapper.mesh

repair = RepairMesh(mesh)
repair.run()

#remove spikes
removeSpikes = RemoveSpikes(mesh)
removeSpikes.maxRounds = 50
removeSpikes.deficiency = 0.001
removeSpikes.run()

#relax
relaxor = Relax(mesh)
relaxor.strength = 1.0
relaxor.iterations = 50
relaxor.fixBoundaries = True
relaxor.run()

geoapp.addModel(mesh, u"Mesh 2")
geoapp.redraw(False)

geo.save()