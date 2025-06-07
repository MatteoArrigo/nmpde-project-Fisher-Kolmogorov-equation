// Set mesh size factor (smaller values = finer mesh)
Mesh.MeshSizeFactor = 1.0;

// Import STL
Merge "../brain-h3.0.stl";

// Create a surface loop from all surfaces
Surface Loop(1) = Surface{:};

// Create a volume from the surface loop
Volume(1) = {1};

// Set mesh algorithm (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, etc.)
Mesh.Algorithm = 6;

// Set 3D mesh algorithm (1=Delaunay, 4=Frontal, 7=MMG3D, etc.)
Mesh.Algorithm3D = 1;

// Set mesh order (1=linear elements, 2=quadratic elements)
Mesh.ElementOrder = 1;

// Generate 3D mesh
Mesh 3;

// Save as MSH format
Save "../data/volume.msh";