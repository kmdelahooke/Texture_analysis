import geomagic.app.v3
from geomagic.app.v3.imports import *

#downsample
geo.select_all_objects()
geo.sample(0.0005, 7, False, 0)
geo.save()
geo.clear_all()

#clean
geo.repair_point_normals(0)
geo.save()

#merge scans
geo.global_registration(0, 100, 2000, False, 20, True, True, False)
geo.merge_point_objects(u'Combined Points 1', False, False)

#clean
geo.reduce_noise(0, 1, True, 0.000597154, 0, False, 1, 0.00249936, False, u'')
geo.delete()
geo.select_outliers(100)
geo.delete()

# Get the active model.
activeModel = geoapp.getActiveModel()

# Get the points from the active model.
pts = geoapp.getPoints(activeModel)

ptSelection = PointSelection(pts)
ptSelection.selectAll()

#wrap
wrapper = WrapPoints()
wrapper.points = pts
wrapper.run()
mesh = wrapper.mesh

repair = RepairMesh(mesh)
repair.run()

geoapp.addModel(mesh, u"Mesh 1")
geoapp.redraw(False)