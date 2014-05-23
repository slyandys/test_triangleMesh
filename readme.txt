Tracking code by Stefanie Wuhrer, work in collaboration with Jochen Lang and Chang Shu. Published in:
S. Wuhrer, J. Lang, C. Shu.
Tracking Complete Deformable Objects with Finite Elements.
In Conference on 3D Imaging Modeling Processing Visualization and Transmission, 2012. 

Template-based tracking approach of data (input needs to be either a point cloud with normals or a triangle mesh). Approach deforms the template to the data using an energy optimization scheme, then uses the deformed areas as input to a
linear FEM in order to predict the unpbserved side, and finally, performs another optimization on the unobserved side.

Library dependencies:
- GeomLib (NRC library to store and access triangle meshes and other geometric structures)
- TriangeMeshLib (NRC library to process triangle meshes, e.g. compute geodesics, Spin images, curvatures, etc.)
- tetgen (used to tetrahedralize mesh before FEM step; version 1.4.3, see http://tetgen.berlios.de/)
- CLAPACK (used for matrix operations in FEM code; version I use is from 2000, see http://www.netlib.org/clapack/faq.html)
- Insight toolkit (used for energy optimization; version 3.16.0, see http://www.itk.org/)
- ANN (compute nearest neighbors; version 1.1.1, see http://www.cs.umd.edu/~mount/ANN/)



Some command line arguments for tracking WITHOUT FEM:

"-track" "Data\\EarModel\\ear_aligned_fixed_20K.wrl" "Data\\EarModel\\Bumblebee\\1\\frames_with_normals.txt" "n" "ear_1" "5.0"

"-track" "Data\\neck\\template_closed.wrl" "Data\\neck\\frames_with_normals_every_2_to_end.txt" "n" "neck" "5.0" "18" "C:\\Stefanie\\Projects\\Tracking\\Data\\neck\\landmark_file.txt"
"-track" "Data\\dinoKinectAT\\template_closed.wrl" "Data\\dinoKinectAT\\frames_with_normals_every_2_to_end.txt" "n" "neckKinect" "5.0"
"-track" "Data\\flap\\template_closed_def.wrl" "Data\\flap\\frames_with_normals.txt" "n" "flap" "5.0"

"-track" "Data\\GroundTruthDataSets\\Gipshand\\input_param\\template.wrl" "Data\\GroundTruthDataSets\\Gipshand\\gipshand_linear\\frames_enhanced.txt" "y" "handLinear" "1.5"
"-track" "Data\\GroundTruthDataSets\\Gipshand\\input_param\\template.wrl" "Data\\GroundTruthDataSets\\Gipshand\\gipshand_NonLinear1\\frames_enhanced.txt" "y" "handNonLinear1" "1.5"

"-track" "Data\\GroundTruthDataSets\\Buste\\input_param\\template.wrl" "Data\\GroundTruthDataSets\\Buste\\buste_linear\\frames.txt" "y" "busteLinear" "1.5"
"-track" "Data\\GroundTruthDataSets\\Buste\\input_param\\template.wrl" "Data\\GroundTruthDataSets\\Buste\\buste_linear\\frames_noise.txt" "y" "busteLinearNoise" "1.5"
"-track" "Data\\GroundTruthDataSets\\Buste\\input_param\\template.wrl" "Data\\GroundTruthDataSets\\Buste\\buste_linear\\frames_outliers.txt" "y" "busteLinearOutliers" "1.5"
"-track" "Data\\GroundTruthDataSets\\Buste\\input_param\\template.wrl" "Data\\GroundTruthDataSets\\Buste\\buste_linear\\frames_res.txt" "y" "busteLinearRes" "1.5"
"-track" "Data\\GroundTruthDataSets\\Buste\\input_param\\template.wrl" "Data\\GroundTruthDataSets\\Buste\\buste_NonLinear1\\frames.txt" "y" "busteNonLinear1" "1.5"

"-track" "Data\\GroundTruthDataSets\\Duck\\duckMediumCoarse_input_param\\template.wrl" "Data\\GroundTruthDataSets\\Duck\\duckMediumCoarse_nonlinear\\frames.txt" "y" "duck" "1.5"



Some command line arguments for tracking WITH FEM:

"-femtrack" "Data\\EarModel\\ear_aligned_fixed_20K.wrl" "Data\\EarModel\\Bumblebee\\1\\frames_with_normals.txt" "n" "fem_ear_1" "Data\\EarModel\\Bumblebee\\1\\contact_point_in_tetrahedral.txt" "5.0" "Data\\EarModel\\ear_aligned_fixed_3K.wrl"

"-femtrack" "Data\\neck\\template_closed.wrl" "Data\\neck\\frames_with_normals_every_2_to_end.txt" "n" "fem_neck" "Data\\neck\\contact_point_in_tetrahedral.txt" "5.0" "Data\\neck\\template_simplified_3000.wrl" "18" "C:\\Stefanie\\Projects\\Tracking\\Data\\neck\\landmark_file.txt"
"-femtrack" "Data\\dinoKinectAT\\template_closed.wrl" "Data\\dinoKinectAT\\frames_with_normals_every_2_to_end.txt" "n" "fem_neckKinect" "Data\\dinoKinectAT\\contact_point_in_tetrahedral.txt" "5.0" "Data\\dinoKinectAT\\template_simplified_3000.wrl"
"-femtrack" "Data\\flap\\template_closed_def.wrl" "Data\\flap\\frames_with_normals.txt" "n" "fem_flap" "Data\\flap\\contact_point_in_tetrahedral.txt" "5.0" "Data\\flap\\template_closed_def_simplified_3000.wrl"

"-femtrack" "Data\\GroundTruthDataSets\\Gipshand\\input_param\\template.wrl" "Data\\GroundTruthDataSets\\Gipshand\\gipshand_linear\\frames_enhanced.txt" "y" "fem_handLinear" "Data\\GroundTruthDataSets\\Gipshand\\input_param\\contact_point_in_tetrahedral.txt" "1.5" "Data\\GroundTruthDataSets\\Gipshand\\input_param\\template_simplified_1000.wrl"
"-femtrack" "Data\\GroundTruthDataSets\\Gipshand\\input_param\template.wrl" "Data\\GroundTruthDataSets\\Gipshand\\gipshand_NonLinear1\\frames_enhanced.txt" "y" "fem_handNonLinear1" "Data\\GroundTruthDataSets\\Gipshand\\input_param\\contact_point_in_tetrahedral.txt" "1.5" "Data\\GroundTruthDataSets\\Gipshand\\input_param\\template_simplified_1000.wrl"

"-femtrack" "Data\\GroundTruthDataSets\\Buste\\input_param\\template.wrl" "Data\\GroundTruthDataSets\\Buste\\buste_linear\\frames.txt" "y" "fem_busteLinear" "Data\\GroundTruthDataSets\\Buste\\input_param\\contact_point_in_tetrahedral.txt" "1.5" "Data\\GroundTruthDataSets\\Buste\\input_param\\template_simplified.wrl"
"-femtrack" "Data\\GroundTruthDataSets\\Buste\\input_param\\template.wrl" "Data\\GroundTruthDataSets\\Buste\\buste_linear\\frames_noise.txt" "y" "fem_busteLinearNoise" "Data\\GroundTruthDataSets\\Buste\\input_param\\contact_point_in_tetrahedral.txt" "1.5" "Data\\GroundTruthDataSets\\Buste\\input_param\\template_simplified.wrl"
"-femtrack" "Data\\GroundTruthDataSets\\Buste\\input_param\\template.wrl" "Data\\GroundTruthDataSets\\Buste\\buste_linear\\frames_outliers.txt" "y" "fem_busteLinearOutliers" "Data\\GroundTruthDataSets\\Buste\\input_param\\contact_point_in_tetrahedral.txt" "1.5" "Data\\GroundTruthDataSets\\Buste\\input_param\\template_simplified.wrl"
"-femtrack" "Data\\GroundTruthDataSets\\Buste\\input_param\\template.wrl" "Data\\GroundTruthDataSets\\Buste\\buste_linear\\frames_res.txt" "y" "fem_busteLinearRes" "Data\\GroundTruthDataSets\\Buste\\input_param\\contact_point_in_tetrahedral.txt" "1.5" "Data\\GroundTruthDataSets\\Buste\\input_param\\template_simplified.wrl"
"-femtrack" "Data\\GroundTruthDataSets\\Buste\\input_param\\template.wrl" "Data\\GroundTruthDataSets\\Buste\\buste_NonLinear1\\frames.txt" "y" "fem_busteNonLinear1" "Data\\GroundTruthDataSets\\Buste\\input_param\\contact_point_in_tetrahedral.txt" "1.5" "Data\\GroundTruthDataSets\\Buste\\input_param\\template_simplified.wrl"

"-femtrack" "Data\\GroundTruthDataSets\\Duck\\duckMediumCoarse_input_param\\template.wrl" "Data\\GroundTruthDataSets\\Duck\\duckMediumCoarse_nonlinear\\frames.txt" "y" "fem_duck" "Data\\GroundTruthDataSets\\Duck\\duckMediumCoarse_input_param\\contact_point.txt" "1.5" "Data\\GroundTruthDataSets\\Duck\\duckMediumCoarse_input_param\\template_low_res_for_tet.wrl"
