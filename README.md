**This repository contains the code needed to reproduce the analyses of:**

<em>Delahooke, K. M., Liu, A.G, Stephenson, N. P., Kenchington, C.G., and Mitchell, E.G, 2024. Quantification of MISS to untangle early animal - matground interactions (in prep).</em>

All data and results will be made public upon publication.

PART 1: Spatial analysis of texture images of the entire surface MISS topography map, in combination with fossil location data
PART 2: Classification of topography maps of Ivesheadiomorphs using Persistent Homology and Surface Metrology

![image](https://github.com/user-attachments/assets/cddf22b4-7a29-4d9e-b336-e2b543d67d8a)

**Geomagic processing scripts:**

These scripts automate scan processing steps using the Geomagic Wrap 2021 Python API  
These scripts were applied to aligned point clouds, cropped to the area of interest (ivesheadiomorph/grid square).   
These scripts automated the process of merging, cleaning and meshing each cropped point cloud, and can be applied via batch processing.  

<ul>
  <li><em>grid_processing_hi.py</em> automates downsampling, merging cleaning and meshing to produce hi-resolution meshes</li>
  <li><em>grid_processing_lo.py</em> was applied to theoutput of 'grid_processing_hi.py' and automates downsampling (to 0.5%), and smoothes the mesh by removing spikes and relaxing the mesh over 50 iterations.</li>
  <li><em>save_obj_file.py</em> exports meshes as .obj files. Currently the filenames needs changing manually.</li>
</ul>








