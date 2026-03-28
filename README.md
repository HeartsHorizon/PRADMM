# PRADMM
This repository contains the implementation for the AAAI 2026 paper:

"Parallelizable Riemannian Alternating Direction Method of Multipliers for Non-convex Pose Graph Optimization".

Pose graph optimization (PGO) is fundamental to robot perception and navigation systems, serving as the mathematical backbone for solving simultaneous localization and mapping (SLAM). Existing solvers suffer from polynomial growth in computational complexity with graph size, hindering real-time deployment in large-scale scenarios. In this paper, by duplicating variables and introducing equality constraints, we reformulate the problem and propose a Parallelizable Riemannian Alternating Direction Method of Multipliers (PRADMM) to solve it efficiently. Compared with the state-of-the-art methods that usually exhibit polynomial time complexity growth with graph size, PRADMM enables efficient parallel computation across vertices regardless of graph size. Crucially, all subproblems admit closed-form solutions, ensuring PRADMM maintains exceptionally stable performance. Furthermore, by carefully exploiting the structures of the coefficient matrices in the constraints, we establish the global convergence of PRADMM under mild conditions, enabling larger relaxation step sizes within the interval (0,2). Extensive empirical validation on two synthetic datasets and multiple real-world 3D SLAM benchmarks confirms the superior computational performance of PRADMM.

If you find our work useful in your research, please consider citing:

Chen, X., Cui, C., Han, D. and Qi, L. 2026. Parallelizable Riemannian Alternating Direction Method of Multipliers for Non-convex Pose Graph Optimization. Proceedings of the AAAI Conference on Artificial Intelligence. 40, 43 (Mar. 2026), 36811-36819. DOI:https://doi.org/10.1609/aaai.v40i43.41007.
