# GLS_last_version
# GLS
descriptors GLS
This repository has been created to compute the GLS (growing least square) descriptor on point clouds/meshes. 
It is the implementation in C++ of the paper "Relative Scale Estimation and 3D Registration of Multi-Modal Geometry Using Growing Least Squares" written by Nicolas Melldo, Matteo Dellepiane and Roberto Scopigno. (https://www.researchgate.net/publication/286510263_Relative_Scale_Estimation_and_3D_Registration_of_Multi-Modal_Geometry_Using_Growing_Least_Squares)
Our code relies on some functions from the Patate library (avalaible on inria gitlab "https://gitlab.inria.fr/patate/patate"). 
This library have then been replaced by the PONCA library (https://www.irit.fr/recherches/STORM/MelladoNicolas/patate-library/)
We didn't test our code with the new version of the library. That's why we provided in zip file the code of the former library developped by Mellado and it's team. 
To make the GLS code works, the user must set the path to the Patate library in his C++ project. Please, refer to the README included in the Patate zip file to add external dependency (Eigen)

This project comes with data provided for testing purpose. 
Our main code is located in the file "/gls/compute_gls.cpp".
4 applications can be generated with our GLS code : 
   - 1st application : computing GLS descriptors. 
     Uncomment from line 138 to line 254 and Comment from line 256 to end
     Then, to test the code, use the data contained in the folder "data_example_1" and run the following command in a terminal: 
     bunny.ply  bunny_source_interest_point.ply  0.1292 4.9115 500 bunny_profiles.txt 1  <br />
     The first parameter is the original mesh/point cloud.<br />
     The second parameter is a subset of keypoints from the original point cloud. GLS will be computed on these keypoints. <br />
     The third and fourth parameters are respectively the minimum and maximum scale (see "https://www.researchgate.net/publication/286510263_Relative_Scale_Estimation_and_3D_Registration_of_Multi-Modal_Geometry_Using_Growing_Least_Squares" p6)
     A good choice for the minimum scale is to take the minimum of the average edge length (for mesh) (or point distance (for point cloud)) between the two objects to register. 
     A good choice for the maximum scale is to take the bigest bounding box diagonal between the two objects to register. Another valid choice is to take the maximum of the average edge length between the two objects, and multiply it by a big consant (let's say 5). <br />
     The fifth parameter is the number of scale samples (see "https://www.researchgate.net/publication/286510263_Relative_Scale_Estimation_and_3D_Registration_of_Multi-Modal_Geometry_Using_Growing_Least_Squares" p6)<br />
     The sixth parameter is the output filename where the algorithm will write the point profiles (i.e. the descriptor)
     and the optional last parameter is a boolean to set to 1 if we want the algorithm to compute and write the geometric variations. <br />
     
     - 2nd application : Estimate the relative scale between the two object to register & Estimate the matching points between the two objects
      Uncomment from line 256 to line 331 and Comment from line 138 to 254 and from line 334 to end. 
      This application is the continuity of the first one. To make it works, the user must first compute the descriptors on the two objects A and B he wants to register. 
      So, he ill have to run 1rst application on object A, 1rst application on object B, and then use both the generated descriptor files in application 2. 
      The folder "data_example_2" contains the results that one would obtained if he ran application 1 on both source and target meshes contained "data_example_1"
      To test the code, use the data contained in the folder "data_example_2" and run the following command in a terminal: 
      bunny_vsa_complete_profiles.txt  bunny_downsampled_vsa_complete_profiles.txt  bunny_vsa_bunny_downsampled_vsa_matching.txt 1
      The first parameter is the descriptor file of object A (in our case source_bunny). 
      The second parameter is the descriptor file of object B (in our case target_bunny).
      The third parameter is the output filename where the algorithm will write the corresponding pairs of points and the estimated scale factor. 
      The last parameter is optional and must be consistant with what you choose for last parameter in application 1. If you didnt add the optional parameter in application 1, don't add it here. Otherwise, do so. 
      
      - 3rd application : Compute the matching pair after having rescaled upstream the object
        Comment from line 138 to 331 and from line 372 to the end. Uncomment from line 334 to line 369. 
        This third application is an alternative to the second application, in the case where the 2 objects share the same scale (or have been rescaled upstream). 
        Then, this application computes the matching pairs from the generated descriptor files of the two same-scale objects. The folder "data_example_3" contains two subfolders. The subfolder "meshes" contains source and target meshes at the same scale that the user can use to generate profiles files with application 1. The subfolder "profiles" contains the file profiles generated with the data from "meshes". 
        To test the code, use the data contained in the folder "data_example_3/profiles" and run the following command in a terminal:
        profiles\bunny_vsa_complete_profiles.txt  profiles\bunny_downsampled_vsa_complete_profiles.txt bunny_vsa_bunny_downsampled_vsa_rescaled_matching.txt
        The first parameter is the descriptor file of object A (in our case source_bunny).
        The second parameter is the descriptor file of object B (in our case target_bunny).
        The third parameter is the output filename where the algorithm will write the matching pair of points. 
        
      - 4th application : Apply RANSAC algorithm (like in the original paper) to estimate the transform
        Comment from line 138 to line 370. Uncomment from line 372 to the end. 
        This application takes the profiles files (computed in app 1 - with geometric variations -) of objects A and B  in order to estimate the matching points. Then, this subset of matching point is given to the RANSAC algorithm that will output the transform T between the two objects to register. 
        From 2 objects to register, the following steps to get to the estimated transform are : run app 1 (with geometric variation) and run app 4. 
        The folder "data_example_4" contains the matching points files obtained from app 1 with data from "data_example_1"
        To test the code, use the data contained in the folder "data_example_4" and run the following command in a terminal:
        bunny_vsa_1000_profiles_and_geom_var.txt  bunny_downsampled_vsa_1000_profiles_and_geom_var.txt transform.txt
        The first parameter is the descriptor file of object A (in our case source_bunny).
        The second parameter is the descriptor file of object B (in our case target_bunny).
        The third parameter is the output filename where the algorithm will write the estimated transform. 
        
        NB 1 : The 4th application contains other parameters directly set inside the function. 
             The user is asked to change them if he wants. 
             "nb_iterations" : the number of iterations for the RANSAC algorithm
             "max_err_scale" : the maximal scale error authorized for a set of k pairs of matching points. Each pair of points comes with it related scale. So scale can not be coherent for all the matching points. This threshold enables to keep a set of matching points with similar scale. 
             "max_err_reg" : parameter of RANSAC algorithm. At iteration k, the computed transformed is applied only if the registration error is below this threshold. 
             "max_err_norm" : parameter of RANSAC algorithm. In order to find a transform, RANSAC algorithm needs at least 3 matching pairs (minimum number of pairs).
             The threshold "max_error_norm" makes sure that at iteration k, the normals of the 3 points in object A are similar to the normal of their 3 corresponding in object B. If so, the transform is applied; otherwise it isn't. 
             "bbox_diag_ratio" : the ratio of the diagonal of the bounding box of object A (source object) over the diagonal of object B (target object)
             "bbox_diag" : diagonal of the bounding box of object B (target object)
             
         NB 2 :  The 4th application contains a debug mode that the user can comment. This debug mode consists in writing intermediate results in text files.
         The debug mode regards : 
            - the filename "kpairs_filename" (line 422 to 427) associated to the function "write_kpairs" (line 435). "write_kpairs" produces a k-lines and 8 columns file containing, for each row : a new matching pair, for the 3 first columns : the 3D coordinates of the target point , for the 3 following columns : the 3D coordinates of the matching source point, for the 7th column, the estimated scale for the pair of points, and for the 8th column an indication of good (1) or bad (0) correspondance based on the correlation value. 
            - the filename "debug_file" (line 437-440) used in the function "write_ply_file_for_debug" (see file "ransac_scheme.cpp" line 835 and file "IO" line 516-562. This function enables to visualize by color in one of the point cloud (A or B) the set of points that have been seleted as half-pair to compute the RANSAC algorithm (in green color). If one of the point belongs to the pair with highest priority, then this point is colored in red. Otherwise, the rest of the points are in black. 
            - the filename "filename_transform" (line 452 to 459) and used in the function "write_matrix_transform" (line 459) enables to save the iterative results of the transform for different values of lambda, in order to see which lambda value in the metric gives the better results for RANSAC. The parameter lambda comes in a metric for pair selection for RANSAC algorithm. The metric is used to chose k pairs of corresponding points for the RANSAC algorithm, based on both the geometric variation and the distance between the selected points in target point cloud. For a given target point p belonging to 1 out of K corresponding pairs, the metric is the following:  _sum_K(sqrt(distL2(p, q_k) / bbox_diag) + lambda*((max_cost_ - diss_cost(p)*diss_cost(q_k))/(max_cost_ - min_cost_))_. 
            
            With :  - sum_k : the sum over the K considered matching pairs
                    - q_k : the target point of the kth pair (k goes from 1 to K)
                    - bbox_diag : diagonal of the bounding box for the target point cloud
                    - lambda : parameter to evaluate in the debug mode (nevertheless, the user can set it and remove the debug mode)
                    - diss_cost(p) : dissimilarity cost (based on the geometric variation) of the pair (p, p'). With p' the correponding point of p in the source point cloud. 
                    - diss_cost(q_k) : dissimilarity cost (based on the geometric variation) of the pair (q_k, q_k'). With q_k' the correponding point of q_k in the  source point cloud. 
                    - max_cost and min_cost : the maximum and minimum dissimilarity cost for all the pairs of matching points. 

         For simplicity, the user can commment the code lines used for debug. Otherwise, he will have to replace our paths and filnemanes for debug files by its own. 
         
        
        NB 3 : In the 4th application, we used the function "pop_k_farthest_pairs" (in "ransac_scheme.cpp") to compute the k pairs of mathing points to give to the RANSAC algorithm. This function is based on the lamdab-metric using both the distance and the geometric variation and described just above. If the user doesn't want to use a simpler metric, he can call "pop_triplet" or "pop_3_farthest_pairs" in replace of "popk_farthest_pairs" in the function "ransac_algorithm". "pop_triplet" will pop a triplet of pairs of matching point according to their ascending dissimilarity cost (i.e. the pair with the smallest dissimilarity cost will be poped first, and so on). "pop_3_farthest_pairs" will first pop the pair of points (p, p') that have the smallest dissimilarity cost (i.e. the one in the top of the priority queue). Then, two other pairs of points (q, q') and (r, r') will complete it by maximizing the distance |p - q| + |p - r| + |q - r| (with p, q and r belonging to the target point cloud). There is an other mode to "pop_3_farthest_pairs" where the 3 pairs of points will be selected according to both the dissimilarity cost (computed with geometric variation) and the distance to each other. This mode needs the user to set a lambda parameter. 
          
              
         
         
