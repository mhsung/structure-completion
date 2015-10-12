### Data-driven Structural Priors for Shape Completion<br>
Minhyuk Sung, Vladimir G. Kim, Roland Angst, and Leonidas Guibas<br>
SIGGRAPH Asia 2015<br>
<p>

**Citation:**<br>
Please cite our paper if you use this code:<br>
>Minhyuk Sung, Vladimir G. Kim, Roland Angst, and Leonidas Guibas<br>
Data-driven Structural Priors for Shape Completion,<br>
ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2015)<br>
<p>

=================

### Installation
See 'lib/install_libraries.sh' for required libraries.


### Prepare a new dataset

0. Prepare mesh files and label files (*.off, *.gt)<br>
Make sure that all mesh files have *UNIT LENGTH* from the bounding box center to the farthest point.

1. Create a mesh and label file directory for the new dataset<br>
`($shape2pose)/data/1_input/($dataset_name)`<br>
This directory should have both `off` and `gt` directories, which are mesh(.off) and label(.gt) file directories, respectively.

2. Create a information file directory for the new dataset<br>
`($shape2pose)/data/0_body/($dataset_name)`<br>
This directory should have following files:<br>

    * `regions.txt`: This file is for both *shape2pose* code and *cuboid-prediction* code.<br>
Each line shows part `(part_name) pnts 1`<br>
The first line part corresponds to label number 0, the next line corresponds to label number 1, and so on.<br>

    * `regions_symmetry.txt`: This file is for both *shape2pose* code and *cuboid-prediction* code.<br>
Each line shows a set of symmetric parts `(part_name_1) (part_name_2) ... (part_name_k)`<br>
All part names should appear in `regions.txt` file.<br>
If a part has no symmetric parts, it should also be written in a line without any other part name.<br>

    * `symmetry_groups.txt`: This file is only for *cuboid-prediction* code.<br>
In contrast to `regions_symmetry.txt` file which has symmetric part information for learning local classifiers,<br>
`symmetry_groups.txt` file has information of symmetric parts in terms of the part structure.<br>
For example, all legs of chairs are considered as symmetric each other when training local classifiers, but the front legs and rear legs are considered as separate symmetric groups in the part structure. <br>

       Each symmetry group should be recorded in the following format:<br>
       ```
       symmetry_group (rotation/reflection) (axis_index:[0,1,2])
       single_label_indices (label_number_0 label_number_1 ... label_number_k)
       pair_label_indices (label_number_pair_a_0 label_number_pair_b_0 ... label_number_pair_a_k label_number_pair_b_k)
       ```
       A symmetry group can be either rotation group or reflection group.<br>
       Axis index means part local axis index ([x, y, z] &rarr; [0, 1, 2]) of reflection plane normal or rotation axis.<br>
       Single label indices indicate single parts which are symmetric in terms of the symmetry axis.<br>
       Pair label indices indicate pair of parts which are symmetric in terms of the symmetry axis.<br>
       Each pair of successive label numbers `label_number_pair_a_i label_number_pair_b_i` shows a (i-th) pair.<br>
       A rotation symmetry group *MUST NOT* has pair label indices (currently not supported),<br>
       and also a reflection symmetry group *MAY NOT* have pair label indices (optional).<br>

       Ex)
       ```
       symmetry_group reflection 0
       single_label_indices 0
       pair_label_indices 1 2 3 4
       ```
       Let [+, +, +] indicates that a cuboid corner of part, which has position<br>
       (center_x + 0.5 * size_x, center_y + 0.5 * size_y, center_z + 0.5 * size_z)<br>
       The [-, +, +] corner of part with label 1 is symmetric with the [+, +, +] corner of part with label 2.<br>
<br>


### Train/test local point classifiers

Make sure that all files are prepared as mentioned above.

1. **[IMPORTANT]**<br>
Copy `regions.txt` and `regions_symmetry.txt` files in `($shape2pose)/data/0_body/($dataset_name)` to `($shape2pose)/data/0_body`.<br>
Double check these files.

2. Make a `($dataset_name)_all.txt` file in `($shape2pose)/script/examples`.<br>
This file should have all mesh file names (without extension).<br>
If you train a subset of files, make the name list file with the subset.

3. **[IMPORTANT]**<br>
If you do cross-validation, copy `trainRegions_cv.py` and `predictRegions_cv.py` files to `trainRegions.py` and `predictRegions.py` in `($shape2pose)/script/scriptlibs`, respectively. <br> 
If you run for the subset of meshes and don't do cross-validation, copy `trainRegions_origin.py` and `predictRegions_origin.py` files to `trainRegions.py` and `predictRegions.py` in `($shape2pose)/script/scriptlibs`, respectively. <br>

4. For training, run the following command in `($shape2pose)/script`:<br>
`./train.py ($dataset_name) exp1_($dataset_name) examples/($dataset_name)_all.txt`<br><br>
The classifier files (`train_(part_name).arff`, `weka_(part_name).model`) are generated in<br>
`($shape2pose)/data/3_trained/classifier/exp1_($dataset_name)`<br>
If you do cross-validation, the classifier files are generated in each mesh name directory.<br>
Make sure that all mesh name directories have the same number of files (classifier files for all parts).

5. For testing, run the following command in `($shape2pose)/script`:<br>
`./test.py ($dataset_name) exp1_($dataset_name) examples/($dataset_name)_all.txt`<br><br>
The prediction files are generated in<br>
`($shape2pose)/data/4_experiments/classifier/exp1_($dataset_name)`<br>
<br>


### Run experiments

0. Compile code<br>
In `../../build/OSMesaViewer/build`, `make`.<br>

1. Make an experiment directory<br>
`($cuboid-prediction)/experiments/exp1_($dataset_name)`<br>
This directory should have following files:<br>
    * `arguments.txt`: The following is the example of arguments.<br>

         ```
         --data_root_path=($shape2pose)
         --label_info_path=data/0_body/($dataset_name)/
         --mesh_path=data/1_input/($dataset_name)/off/
         --sample_path=data/2_analysis/($dataset_name)/points/even1000/
         --dense_sample_path=data/2_analysis/($dataset_name)/points/random100000/
         --mesh_label_path=data/1_input/($dataset_name)/gt/
         --sample_label_path=data/4_experiments/exp1_($dataset_name)/1_prediction/
         ```
    **[IMPORTANT]**<br>
    Make sure that `data_root_path` is set correctly.
    * `pose.txt`: Camera pose file for rendering.<br>

2. Run experiments<br>
In `($cuboid-prediction)/python`, run the following command:<br>
`./batch_exec.py ($exp_type) ($shape2pose)/data/1_input/($dataset_name)/off/ ../experiments/($dataset_name)/`<br><br>
Run the command in the following `($exp_type)` order:<br>

    1) **ground_truth_cuboids**: Create ground truth cuboids in `../experiments/($dataset_name)/training`<br>
        **[IMPORTANT]**<br>
        After this, run the following command in `../experiments/($dataset_name)` for generating part relation statistics files:<br>
        `../../build/OSMesaViewer/build/Build/bin/OSMesaViewer --run_training --flagfile=arguments.txt`<br>

    2) **prediction**: Run our method. Files are generated in `../experiments/($dataset_name)/output`<br>
    3) **part_assembly**: Run part assembly. Files are generated in `../experiments/($dataset_name)/part_assembly`<br>
    4) **symmetry_detection**: Run symmetry detection. Files are generated in `../experiments/($dataset_name)/symmetry_detection`<br>
        **[IMPORTANT]**<br>
        *BEFORE* this, run the following command in `($cuboid-prediction)/python`:<br>
        `./batch_symmetry_detection.py ($shape2pose)/data/1_input/($dataset_name)/off/ 
../experiments/($dataset_name)/`<br>
        Make sure that `binDir` variable in `./batch_symmetry_detection.py` is correctly set.

    5) **baseline**: Compute baseline. Files are generated in `../experiments/($dataset_name)/baseline`<br>
    6) **render_assembly**: Render part assembly cuboids. Should be executed after part assembly.<br>
    7) **render_evaluation** [optional]: Re-render all experimental result images (including part assembly and symmetry detection). Used when rendering with new parameters.<br>
    8) **extract_symmetry_info** [optional]: Used when extracting symmetry axes information of our method results.<br>

3. Generate HTML result pages<br>
In `($cuboid-prediction)/python`, run the following command:<br>
`./generate_all.sh  ../experiments/($dataset_name)/`<br>
<br>
The resulting files are created in the following directories:<br>
`($cuboid-prediction)/output/exp1_($dataset_name)`<br>
`($cuboid-prediction)/output/assembly_exp1_($dataset_name)`<br>
`($cuboid-prediction)/output/symm_detection_exp1_($dataset_name)`<br>
`($cuboid-prediction)/output/baseline_exp1_($dataset_name)`<br>
<br>
For generating paper figures, run the following command in `($cuboid-prediction)/python/figures`:<br>
`./fig_N.py fig_N.txt`<br>
Select examples and record in `fig_N.txt` files.<br>
Tex files and relates image files are generated in `($cuboid-prediction)/report` and `($cuboid-prediction)/report/images`.<br>
<br>


### Parameters

The followings are remarkable parameters (can be set by adding in the `arguments.txt` file):<br>

* occlusion_pose_filename:<br>
**[IMPORTANT]** If this is set to "", random occlusion pose is generated based on `random_view_seed`.<br>
* random_view_seed:<br>
Seed number of random occlusion view generation.<br>
* param_min_num_symmetric_point_pairs:<br>
If the number of symmetric point pairs is less than this value, the symmetric point pairs are not considered in optimization.<br>
* param_min_sample_point_confidence:<br>
In the initial step, only points with confidence greater than this value are clustered. A lower value can be better when there are noise in the input points.<br>
* param_sparse_neighbor_distance:<br>
Point neighbor distance (in most cases).<br>
* param_cuboid_split_neighbor_distance:<br>
Point neighbor distance for splitting initial cuboids.<br>
* param_occlusion_test_neighbor_distance:<br>
Point neighbor distance for occlusion test.<br>
* param_fusion_grid_size:<br>
Voxel size for fusion.<br>
* param_min_cuboid_overall_visibility:<br>
If the overall visibility of the missing cuboid is greater than this value, it is considered as created in the visible area, and ignored.<br>
* param_fusion_visibility_smoothing_prior:<br>
MRF smoothing prior value for fusion. <br>
* param_eval_min_neighbor_distance:<br>
Minimum error value for accuracy/completeness rendering.<br>
Run `render_evaluation` for rendering with new parameter values.<br>
* param_eval_max_neighbor_distance:<br>
Maximum error value for accuracy/completeness rendering.<br>
Run `render_evaluation` for rendering with new parameter values.<br>
* use_view_plane_mask:<br>
**[IMPORTANT]** Set true if one uses view plane 2D occlusion mask.<br>
* param_view_plane_mask_proportion:<br>
The view plane 2D occlusion mask is created so that this proportion of points are occluded more *AFTER* self-occlusion.<br>
<br>


### Experiment result files
`($shape2pose)/data/0_body/`, `($shape2pose)/data/1_input/`<br>
assembly_airplanes<br>
assembly_bicycles<br>
assembly_chairs<br>
coseg_chairs<br>
shapenet_tables<br>
scan_chairs<br>

`($shape2pose)/scripts/examples`<br>
assembly_airplanes_all.txt<br>
assembly_bicycles_all.txt<br>
assembly_chairs_all.txt<br>
coseg_chairs_all.txt<br>
shapenet_tables_all.txt<br>
scan_chairs_all.txt<br>
coseg_chairs_all_N_train.txt, coseg_chairs_all_N_test.txt<br>
('N' is proportion of training mesh files).

`($shape2pose)/data/3_trained/classifier`, `($shape2pose)/data/4_experiments`<br>
exp1_cv_assembly_airplanes_all.txt<br>
exp1_cv_assembly_bicycles_all.txt<br>
exp1_cv_assembly_chairs_all.txt<br>
exp1_cv_coseg_chairs_all.txt<br>
exp1_cv_shapenet_tables_all.txt<br>
exp1_scan_chairs_all.txt<br>
exp_coseg_chairs_N<br>

`($cuboid-prediction)/experiments`<br>
final3_assembly_airplanes_all<br>
final3_assembly_bicycles_all<br>
final3_assembly_chairs_all<br>
final3_coseg_chairs_all<br>
final3_shapenet_tables_all<br>
test_scan_chairs<br>
partial_coseg_chairs_N<br>

Results with view plane mask with 30% proportion 2D occlusion:<br>
view_mask_assembly_airplanes_all<br>
view_mask_assembly_bicycles_all<br>
view_mask_assembly_chairs_all<br>
view_mask_coseg_chairs_all<br>
view_mask_shapenet_tables_all<br>
