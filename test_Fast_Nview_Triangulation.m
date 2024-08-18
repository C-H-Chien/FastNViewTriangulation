%> Test certifiable fast multiple view triangulation as a mexFunction in
%  the integration of C++ and MATLAB code. This test example code does a
%  simple 3-view triangulation.

%> 1) Preparing data
K = [512, 0, 512; 0, 512, 512; 0, 0, 1];
%> absolute GT rotations
R1 = [0.896288, -0.375944, 0.235231; 0.437885, 0.834166, -0.335297; -0.070169, 0.403527, 0.912273];
R2 = [0.872666, 0.0423745, 0.486475; -0.0995012, 0.990757, 0.0921907; -0.478072, -0.128857, 0.868817];
R3 = [0.980714, 0.156866, -0.116591; -0.182924, 0.946782, -0.264846; 0.0688409, 0.281065, 0.957216];
%> absolute GT translations
t1 = [0.181751; -0.0476752; -0.158032];
t2 = [0.01608; 0.018602; 0.0462539];
t3 = [1.26283; 1.14514; 0.0474746];
%> Noisy point correspondences in metrics and their conversions to pixels
gamma1 = [0.538634, 0.382500, 1.000000]';
gamma2 = [1.974536, 1.198738, 1.000000]';
gamma3 = [0.707729, 0.379656, 1.000000]';
p1 = K * gamma1;
p2 = K * gamma2;
p3 = K * gamma3;
%> Construct a feature track
Feature_Track = [p1(1:2,:), p2(1:2,:), p3(1:2,:)];

%> 2) Reshape the data for inputs of mex function
K_ = reshape(K', [1,9]);
Rs = [reshape(R1', [1,9]), reshape(R2', [1,9]), reshape(R3', [1,9])];
Ts = [t1; t2; t3];

%> 3) Run multiview triangulation from mexFunction in CPP code
debug = 0;
[corrected_features, reproj_errs, is_sol_global_optimal, ~] = ...
    fast_multiview_triangulation_mex(Feature_Track, K_, Rs, Ts, debug);

N = size(Feature_Track, 2);
%> Check if the returned reprojection erros are correct and making senses
for i = 1:N
    assert(abs(norm(Feature_Track(:,i) - corrected_features(:,i)) - reproj_errs(i)) < 1e-8);
    assert(abs(reproj_errs(i)) < 5);
end