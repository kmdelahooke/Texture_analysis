import geomagic.app.v3
from geomagic.app.v3.imports import *

models = geoapp.getModels()
for model in models:
	if model.name == "Mesh 2":
		geoapp.setActiveModel(model)

activeModel = geoapp.getActiveModel()
mesh = geoapp.getMesh(activeModel)

#write
path = "D:/lo_res_grid_obj"
saveOptions = FileSaveOptions()
saveOptions.units = Length.Millimeters

writer = WriteFile()
writer.filename = path + "/x35y25.obj"
writer.options = saveOptions
writer.mesh = mesh
writer.run()
