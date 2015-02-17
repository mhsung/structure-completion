#include "MeshCuboidSolver.h"

#include "MeshCuboidParameters.h"
#include "Utilities.h"

#include <cstdint>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCore>
#include <MRFEnergy.h>
#include <EigenQP.h>


bool write_eigen_matrix_binary(const char* filename, const Eigen::MatrixXd& _matrix){
	std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
	if (!out.good()) return false;
	int16_t rows = static_cast<int16_t>(_matrix.rows());
	int16_t cols = static_cast<int16_t>(_matrix.cols());
	out.write((char*)(&rows), sizeof(int16_t));
	out.write((char*)(&cols), sizeof(int16_t));
	out.write((char*)_matrix.data(), rows*cols*sizeof(Eigen::MatrixXd::Scalar));
	out.close();
	return true;
}

bool read_eigen_matrix_binary(const char* filename, Eigen::MatrixXd& matrix){
	std::ifstream in(filename, std::ios::in | std::ios::binary);
	if (!in.good()) return false;
	int16_t rows = 0, cols = 0;
	in.read((char*)(&rows), sizeof(int16_t));
	in.read((char*)(&cols), sizeof(int16_t));
	matrix.resize(rows, cols);
	in.read((char *)matrix.data(), rows*cols*sizeof(Eigen::MatrixXd::Scalar));
	in.close();
	return true;
}

bool write_eigen_vector_binary(const char* filename, const Eigen::VectorXd& _vector){
	std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
	if (!out.good()) return false;
	int16_t rows = static_cast<int16_t>(_vector.rows());
	int16_t cols = static_cast<int16_t>(_vector.cols());
	out.write((char*)(&rows), sizeof(int16_t));
	out.write((char*)(&cols), sizeof(int16_t));
	out.write((char*)_vector.data(), rows*cols*sizeof(Eigen::VectorXd::Scalar));
	out.close();
	return true;
}

bool read_eigen_vector_binary(const char* filename, Eigen::VectorXd& matrix){
	std::ifstream in(filename, std::ios::in | std::ios::binary);
	if (!in.good()) return false;
	int16_t rows = 0, cols = 0;
	in.read((char*)(&rows), sizeof(int16_t));
	in.read((char*)(&cols), sizeof(int16_t));
	matrix.resize(rows, cols);
	in.read((char *)matrix.data(), rows*cols*sizeof(Eigen::VectorXd::Scalar));
	in.close();
	return true;
}

std::vector<int> solve_markov_random_field(const unsigned int _num_nodes, const unsigned int _num_labels,
	const Eigen::MatrixXd& _energy_mat)
{
	assert(_energy_mat.rows() == _num_nodes * _num_labels);
	assert(_energy_mat.cols() == _num_nodes * _num_labels);

	MRFEnergy<TypeGeneral>* mrf;
	MRFEnergy<TypeGeneral>::NodeId* nodes;
	MRFEnergy<TypeGeneral>::Options options;
	TypeGeneral::REAL energy, lower_bound;

	std::list<TypeGeneral::REAL *> energy_term_list;
	mrf = new MRFEnergy<TypeGeneral>(TypeGeneral::GlobalSize());
	nodes = new MRFEnergy<TypeGeneral>::NodeId[_num_nodes];


	// Data term.
	for (unsigned int node_index = 0; node_index < _num_nodes; ++node_index)
	{
		TypeGeneral::REAL *single_energy = new TypeGeneral::REAL[_num_labels];
		energy_term_list.push_back(single_energy);

		for (unsigned int label_index = 0; label_index < _num_labels; ++label_index)
		{
			unsigned int mat_index = node_index * _num_labels + label_index;
			Real energy = _energy_mat(mat_index, mat_index);
			single_energy[label_index] = energy;
		}

		nodes[node_index] = mrf->AddNode(TypeGeneral::LocalSize(_num_labels),
			TypeGeneral::NodeData(single_energy));
	}


	// Smoothness term.
	for (unsigned int node_index_1 = 0; node_index_1 < _num_nodes - 1; ++node_index_1)
	{
		for (unsigned int node_index_2 = node_index_1 + 1; node_index_2 < _num_nodes; ++node_index_2)
		{
			TypeGeneral::REAL *pair_energy = new TypeGeneral::REAL[_num_labels * _num_labels];
			energy_term_list.push_back(pair_energy);
			memset(pair_energy, 0, _num_labels * _num_labels * sizeof(TypeGeneral::REAL));

			for (unsigned int label_index_1 = 0; label_index_1 < _num_labels; label_index_1++)
			{
				for (unsigned int label_index_2 = 0; label_index_2 < _num_labels; label_index_2++)
				{
					unsigned int mat_index_1 = node_index_1 * _num_labels + label_index_1;
					unsigned int mat_index_2 = node_index_2 * _num_labels + label_index_2;

					Real energy = _energy_mat(mat_index_1, mat_index_2);
					// NOTE:
					// Check symmetry.
					assert(energy == _energy_mat(mat_index_2, mat_index_1));

					// NOTE:
					// Check out the index. It should be the following:
					// label_index_1 + label_index_2 * num_labels.
					pair_energy[label_index_1 + label_index_2 * _num_labels] = energy;
				}
			}

			mrf->AddEdge(nodes[node_index_1], nodes[node_index_2],
				TypeGeneral::EdgeData(TypeGeneral::GENERAL, pair_energy));
		}
	}

	std::vector<int> output_labels(_num_nodes);


	// Function below is optional - it may help if, for example, nodes are added in a random order
	//mrf->SetAutomaticOrdering();
	options.m_iterMax = 100; // maximum number of iterations
	options.m_printIter = 10;
	options.m_printMinIter = 0;

	//////////////////////// BP algorithm ////////////////////////
	//mrf->ZeroMessages();
	//mrf->AddRandomMessages(0, 0.0, 1.0);
	//mrf->Minimize_BP(options, energy);
	//std::cout << "Energy = " << energy << std::endl;

	/////////////////////// TRW-S algorithm //////////////////////
	mrf->ZeroMessages();
	mrf->AddRandomMessages(0, 0.0, 1.0);
	mrf->Minimize_TRW_S(options, lower_bound, energy);
	std::cout << "Energy = " << energy << std::endl;


	Eigen::VectorXd solution_vec = Eigen::VectorXd::Zero(_num_nodes * _num_labels);
	for (unsigned int node_index = 0; node_index < _num_nodes; ++node_index)
	{
		int label_index = mrf->GetSolution(nodes[node_index]);
		output_labels[node_index] = label_index;
		solution_vec(node_index * _num_labels + label_index) = 1.0;
	}

	double energy_verified = solution_vec.transpose() * _energy_mat * solution_vec;
	std::cout << "Energy [Verified] = " << energy_verified << std::endl;

	for (std::list<TypeGeneral::REAL *>::iterator it = energy_term_list.end();
		it != energy_term_list.end(); ++it)
		delete[](*it);
	delete[] nodes;
	delete mrf;

	return output_labels;
}

/*
std::vector<int> solve_markov_random_field(const unsigned int _num_nodes, const unsigned int _num_labels,
	const Eigen::MatrixXd& _energy_mat)
{
	const unsigned int mat_size = _num_nodes * _num_labels;
	assert(_energy_mat.rows() == mat_size);
	assert(_energy_mat.cols() == mat_size);

	// Find the minimum/maximum potential.
	double min_potential = k_max_potential, max_potential = 0.0;

	for (unsigned int i = 0; i < mat_size; ++i)
	{
		for (unsigned int j = 0; j < mat_size; ++j)
		{
			assert(_energy_mat(i, j) >= 0);
			if (_energy_mat(i, j) < k_max_potential)
			{
				min_potential = std::min(min_potential, _energy_mat(i, j));
				max_potential = std::max(max_potential, _energy_mat(i, j));
			}
		}
	}

	assert(max_potential >= min_potential);
	double potential_range = max_potential - min_potential;
	//double sigma = -(potential_range * potential_range) / std::log(1.0E-6);
	double sigma = -0.25 / std::log(0.01);
	Eigen::MatrixXd exp_energy_mat = Eigen::MatrixXd::Zero(mat_size, mat_size);

	for (unsigned int i = 0; i < mat_size; ++i)
	{
		for (unsigned int j = 0; j < mat_size; ++j)
		{
			assert(_energy_mat(i, j) >= 0);
			if (_energy_mat(i, j) >= k_max_potential)
			{
				exp_energy_mat(i, j) = 0;
			}
			else
			{
				assert(std::abs(_energy_mat(i, j) - _energy_mat(j, i)) <= 1.0E-6);
				double potential = 0.5 * (_energy_mat(i, j) + _energy_mat(j, i));
				potential -= min_potential;
				exp_energy_mat(i, j) = std::exp(-potential * potential / sigma);
				//exp_energy_mat(i, j) = 1.0 / min_potential;
			}
		}
	}

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
	es.compute(exp_energy_mat);

	Eigen::VectorXd solution_vec = es.eigenvectors().col(mat_size - 1);

	std::vector<int> output_labels(_num_nodes);
	for (unsigned int node_index = 0; node_index < _num_nodes; ++node_index)
	{
		Eigen::VectorXd node_solution_vec = solution_vec.segment(
			node_index * _num_labels, _num_labels);
		int max_label = 0;
		node_solution_vec.maxCoeff(&max_label);
		output_labels[node_index] = max_label;
	}

	double energy_verified = solution_vec.transpose() * exp_energy_mat * solution_vec;
	std::cout << "Energy [Verified] = " << energy_verified << std::endl;

	return output_labels;
}
*/

Eigen::VectorXd solve_quadratic_programming(
	const Eigen::MatrixXd& _quadratic_term,
	const Eigen::VectorXd& _linear_term,
	const double _constant_term,
	Eigen::VectorXd* _init_values_vec = NULL,
	Eigen::MatrixXd* _init_values_mask = NULL,
	double _quadprog_ratio = 1.0)
{
	const unsigned int n = _quadratic_term.cols();
	assert(_quadratic_term.rows() == n);
	assert(_linear_term.rows() == n);

	Eigen::VectorXd x = Eigen::VectorXd::Zero(n);

	Eigen::MatrixXd G = _quadratic_term;
	Eigen::VectorXd g0 = _linear_term;
	Eigen::MatrixXd CE(n, 0);
	Eigen::VectorXd ce0(0.0);
	Eigen::MatrixXd CI(n, 0);
	Eigen::VectorXd ci0(0.0);

	std::cout << "Energy: (";

	// Optional.
	if (_init_values_vec)
	{
		assert((*_init_values_vec).rows() == n);

		double init_error = 0;
		Eigen::VectorXd x0 = (*_init_values_vec);
		init_error += x0.transpose() * _quadratic_term * x0;
		init_error += 2 * _linear_term.transpose() * x0;
		init_error += _constant_term;
		std::cout << "initial = " << init_error << ", ";

		//if (_init_values_mask)
		//{
		//	assert((*_init_values_mask).cols() == n);
		//	const unsigned int p = (*_init_values_mask).rows();

		//	CE.resize(n, p);
		//	ce0.resize(p);
		//	CE = (*_init_values_mask).transpose();
		//	ce0 = -(*_init_values_mask) * x0;
		//}
		//else
		//{
		//	// Initial Prior.
		//	if (_quadprog_ratio > 0)
		//	{
		//		std::cout << "quadprog_ratio = " << _quadprog_ratio << std::endl;
		//		G = (1.0 / init_error) * G + _quadprog_ratio * Eigen::MatrixXd::Identity(dimension, dimension);
		//		g0 = (1.0 / init_error) * g0 - _quadprog_ratio * (*_init_values_vec);
		//	}
		//}
	}
	

	// Make the quadratic term symmetric.
	G = 0.5 * (G + G.transpose());

	// Regularizer using initial values.
	//G = G + 1.0E-12 * Eigen::MatrixXd::Identity(dimension, dimension);
	//g0 = g0 - 1.0E-12 * (*_init_values_vec);

	// Direct solution.
	//x = -G.inverse() * g0;
	
	double energy = QP::solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
	//std::cout << "QuadProg++ Error = " << energy << std::endl;

	double final_error = 0;
	final_error = x.transpose() * _quadratic_term * x;
	final_error += 2 * _linear_term.transpose() * x;
	final_error += _constant_term;

	std::cout << "final = " << final_error << ")" << std::endl;

	return x;
}

/*
Eigen::VectorXd solve_quadratic_programming(
	const Eigen::MatrixXd& _quadratic_term,
	const Eigen::VectorXd& _linear_term,
	const double _constant_term,
	Eigen::VectorXd* _init_values_vec = NULL,
	Eigen::MatrixXd* _init_values_mask = NULL,
	double _quadprog_ratio = 1.0)
{
	const unsigned int dimension = _quadratic_term.cols();
	assert(_quadratic_term.rows() == dimension);
	assert(_linear_term.rows() == dimension);

	Eigen::VectorXd x = Eigen::VectorXd::Zero(dimension);

	// Use an external application.
	bool ret;
	ret = write_eigen_matrix_binary("quadratic_term.dat", _quadratic_term);
	assert(ret);

	ret = write_eigen_vector_binary("linear_term.dat", _linear_term);
	assert(ret);

	if (_init_values_vec)
	{
		ret = write_eigen_vector_binary("init_values_vec.dat", *_init_values_vec);
		assert(ret);
	}

	std::cout << "Waiting for the result from quadprog_solver... " << std::endl;

	QProcess process;
	process.setStandardOutputFile("log_quadprog_solver.txt");
	process.start("quadprog_solver.exe");
	process.waitForFinished();
	qDebug() << process.readAllStandardOutput();
	process.close();

	while (!read_eigen_vector_binary("quadprog_solver.dat", x));
	assert(x.rows() == dimension);

	QFile::remove("quadratic_term.dat");
	QFile::remove("linear_term.dat");
	QFile::remove("init_values_vec.dat");
	QFile::remove("quadprog_solver.dat");
	std::cout << "Done." << std::endl;

	std::cout << "Energy: (";

	// Optional.
	if (_init_values_vec)
	{
		double init_error = 0;
		Eigen::VectorXd x0 = (*_init_values_vec);
		init_error += x0.transpose() * _quadratic_term * x0;
		init_error += 2 * _linear_term.transpose() * x0;
		init_error += _constant_term;
		std::cout << "initial = " << init_error << ", ";
	}

	double final_error = 0;
	final_error = x.transpose() * _quadratic_term * x;
	final_error += 2 * _linear_term.transpose() * x;
	final_error += _constant_term;

	std::cout << "final = " << final_error << ")" << std::endl;

	return x;
}
*/

void update_cuboid_surface_points(
	MeshCuboidStructure &_cuboid_structure,
	const Real _modelview_matrix[16])
{
	const Real radius = FLAGS_param_observed_point_radius * _cuboid_structure.mesh_->get_object_diameter();

	std::vector<MeshCuboid *> all_cuboids = _cuboid_structure.get_all_cuboids();
	for (std::vector<MeshCuboid *>::iterator it = all_cuboids.begin(); it != all_cuboids.end(); ++it)
	{
		MeshCuboid *cuboid = (*it);
		cuboid->create_grid_points_on_cuboid_surface(
			FLAGS_param_num_cuboid_surface_points);

		cuboid->compute_cuboid_surface_point_visibility(
			_modelview_matrix, radius, _cuboid_structure.sample_points_);
	}
}

void segment_sample_points(
	MeshCuboidStructure &_cuboid_structure)
{
	// Parameter.
	const double null_cuboid_probability = 0.1;
	
	assert(_cuboid_structure.mesh_);
	double neighbor_distance = FLAGS_param_sample_point_neighbor_distance *
		_cuboid_structure.mesh_->get_object_diameter();
	double lambda = -(neighbor_distance * neighbor_distance)
		/ std::log(null_cuboid_probability);
	double min_probability = 1 / FLAGS_param_max_potential;

	unsigned int num_sample_points = _cuboid_structure.num_sample_points();
	const int num_neighbors = std::min(FLAGS_param_num_sample_point_neighbors,
		static_cast<int>(num_sample_points));

	std::vector<MeshCuboid *> all_cuboids = _cuboid_structure.get_all_cuboids();
	unsigned int num_cuboids = all_cuboids.size();



	//
	std::vector<ANNpointArray> cuboid_ann_points(num_cuboids);
	std::vector<ANNkd_tree*> cuboid_ann_kd_tree(num_cuboids);

	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
	{
		MeshCuboid *cuboid = all_cuboids[cuboid_index];

		unsigned int num_cuboid_surface_points = cuboid->num_cuboid_surface_points();
		Eigen::MatrixXd cuboid_surface_points(3, num_cuboid_surface_points);
		assert(num_cuboid_surface_points > 0);

		for (unsigned int point_index = 0; point_index < num_cuboid_surface_points; ++point_index)
		{
			for (unsigned int i = 0; i < 3; ++i)
				cuboid_surface_points.col(point_index)(i) =
				cuboid->get_cuboid_surface_point(point_index)->point_[i];
		}

		cuboid_ann_kd_tree[cuboid_index] = ICP::create_kd_tree(cuboid_surface_points,
			cuboid_ann_points[cuboid_index]);
		assert(cuboid_ann_points[cuboid_index]);
		assert(cuboid_ann_kd_tree[cuboid_index]);
	}

	ANNpoint q = annAllocPt(3);
	ANNidxArray nn_idx = new ANNidx[1];
	ANNdistArray dd = new ANNdist[1];
	//



	// Single potential.
	// NOTE: The last column is for the null cuboid.
	Eigen::MatrixXd single_potentials(num_sample_points, num_cuboids + 1);

	for (unsigned int point_index = 0; point_index < num_sample_points; ++point_index)
	{
		MeshSamplePoint *sample_point = _cuboid_structure.sample_points_[point_index];
		for (unsigned int i = 0; i < 3; ++i)
			q[i] = sample_point->point_[i];

		for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
		{
			MeshCuboid *cuboid = all_cuboids[cuboid_index];
			unsigned int label_index = cuboid->get_label_index();

			assert(cuboid_ann_kd_tree[cuboid_index]);
			cuboid_ann_kd_tree[cuboid_index]->annkSearch(q, 1, nn_idx, dd);
			double distance = dd[0];
			assert(distance >= 0);

			double label_probability = sample_point->label_index_confidence_[label_index];
			double energy = distance * distance - lambda * std::log(label_probability);

			if (cuboid->is_group_cuboid())
				energy = FLAGS_param_max_potential;

			single_potentials(point_index, cuboid_index) = energy;
		}

		// For null cuboid.
		double energy = FLAGS_param_max_potential;
		single_potentials(point_index, num_cuboids) = energy;
	}

	// Deallocate ANN.
	annDeallocPt(q);
	delete[] nn_idx;
	delete[] dd;

	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
	{
		annDeallocPts(cuboid_ann_points[cuboid_index]);
		delete cuboid_ann_kd_tree[cuboid_index];
	}


	// Construct a KD-tree.
	Eigen::MatrixXd sample_points(3, num_sample_points);
	for (SamplePointIndex point_index = 0; point_index < num_sample_points; ++point_index)
	{
		for (unsigned int i = 0; i < 3; ++i)
			sample_points.col(point_index)(i) =
			_cuboid_structure.sample_points_[point_index]->point_[i];
	}

	const int dim = 3;
	ANNpointArray sample_ann_points = annAllocPts(num_sample_points, dim);	// allocate data points

	for (SamplePointIndex point_index = 0; point_index < num_sample_points; ++point_index)
	{
		for (int i = 0; i < dim; i++)
			sample_ann_points[point_index][i] = sample_points.col(point_index)[i];
	}

	ANNkd_tree *sample_kd_tree = new ANNkd_tree(sample_ann_points, num_sample_points, dim);
	q = annAllocPt(dim);
	nn_idx = new ANNidx[num_neighbors];
	dd = new ANNdist[num_neighbors];


	// Pair potentials.
	std::vector< Eigen::Triplet<double> > pair_potentials;
	pair_potentials.reserve(num_sample_points * num_neighbors);

	for (unsigned int point_index = 0; point_index < num_sample_points; ++point_index)
	{
		for (unsigned int i = 0; i < 3; ++i)
			q[i] = sample_points.col(point_index)[i];

		int num_searched_neighbors = sample_kd_tree->annkFRSearch(q,
			neighbor_distance, num_neighbors, nn_idx, dd);

		for (int i = 0; i < std::min(num_neighbors, num_searched_neighbors); i++)
		{
			unsigned int n_point_index = (int)nn_idx[i];

			// NOTE: Avoid symmetric pairs.
			if (n_point_index <= point_index)
				continue;

			//
			double distance = (neighbor_distance - dd[i]);
			assert(distance >= 0);
			//

			double energy = distance * distance;
			pair_potentials.push_back(Eigen::Triplet<double>(point_index, n_point_index, energy));
		}
	}

	delete[] nn_idx;
	delete[] dd;
	annDeallocPt(q);
	annDeallocPts(sample_ann_points);
	delete sample_kd_tree;
	annClose();


	// MRF.
	MRFEnergy<TypePotts>* mrf;
	MRFEnergy<TypePotts>::NodeId* nodes;
	MRFEnergy<TypePotts>::Options options;
	TypePotts::REAL energy, lower_bound;

	const int num_nodes = single_potentials.rows();
	const int num_labels = single_potentials.cols();

	std::list<TypeGeneral::REAL *> energy_term_list;
	mrf = new MRFEnergy<TypePotts>(TypePotts::GlobalSize(num_labels));
	nodes = new MRFEnergy<TypePotts>::NodeId[num_nodes];

	// Data term.
	for (unsigned int node_index = 0; node_index < num_nodes; ++node_index)
	{
		TypeGeneral::REAL *D = new TypeGeneral::REAL[num_labels];
		energy_term_list.push_back(D);

		for (unsigned int label_index = 0; label_index < num_labels; ++label_index)
			D[label_index] = static_cast<TypeGeneral::REAL>(
				single_potentials(node_index, label_index));
		nodes[node_index] = mrf->AddNode(TypePotts::LocalSize(), TypePotts::NodeData(D));
	}

	// Smoothness term.
	for (std::vector< Eigen::Triplet<double> >::iterator it = pair_potentials.begin();
		it != pair_potentials.end(); ++it)
	{
		unsigned int node_index_i = (*it).row();
		assert(node_index_i < num_nodes);
		unsigned int node_index_j = (*it).col();
		assert(node_index_j < num_nodes);
		double potential = (*it).value();
		mrf->AddEdge(nodes[node_index_i], nodes[node_index_j], TypePotts::EdgeData(potential));
	}


	// Function below is optional - it may help if, for example, nodes are added in a random order
	//mrf->SetAutomaticOrdering();
	options.m_iterMax = 100; // maximum number of iterations
	options.m_printIter = 10;
	options.m_printMinIter = 0;

	//////////////////////// BP algorithm ////////////////////////
	//mrf->ZeroMessages();
	//mrf->AddRandomMessages(0, 0.0, 1.0);
	//mrf->Minimize_BP(options, energy);
	//std::cout << "Energy = " << energy << std::endl;

	/////////////////////// TRW-S algorithm //////////////////////
	mrf->ZeroMessages();
	mrf->AddRandomMessages(0, 0.0, 1.0);
	mrf->Minimize_TRW_S(options, lower_bound, energy);
	std::cout << "Energy = " << energy << std::endl;


	std::vector<int> output_labels(num_nodes);
	for (unsigned int node_index = 0; node_index < num_nodes; ++node_index)
		output_labels[node_index] = mrf->GetSolution(nodes[node_index]);

	for (std::list<TypeGeneral::REAL *>::iterator it = energy_term_list.end();
		it != energy_term_list.end(); ++it)
		delete[](*it);
	delete[] nodes;
	delete mrf;


	// Reassign sample points to cuboids.
	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
	{
		MeshCuboid *cuboid = all_cuboids[cuboid_index];
		cuboid->clear_sample_points();
	}

	for (unsigned int point_index = 0; point_index < num_sample_points; ++point_index)
	{
		MeshSamplePoint *sample_point = _cuboid_structure.sample_points_[point_index];
		int cuboid_index = output_labels[point_index];
		//int cuboid_index;
		//single_potentials.row(point_index).minCoeff(&cuboid_index);

		if (cuboid_index < num_cuboids)
			all_cuboids[cuboid_index]->add_sample_point(sample_point);
	}

	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
	{
		MeshCuboid *cuboid = all_cuboids[cuboid_index];
		std::cout << "[" << cuboid_index << "]: " << cuboid->num_sample_points() << std::endl;
		cuboid->update_point_correspondences();
	}
}

void compute_labels_and_axes_configuration_potentials(
	const std::vector<Label>& _labels,
	const std::vector<MeshCuboid *>& _cuboids,
	const MeshCuboidPredictor &_predictor,
	Eigen::MatrixXd &_potential_mat,
	const std::vector< std::list<LabelIndex> > *_label_symmetries,
	bool _add_dummy_label)
{
	unsigned int num_labels = _labels.size();
	unsigned int num_cuboids = _cuboids.size();
	unsigned int num_axis_configurations = MeshCuboid::num_axis_configurations();
	assert(num_axis_configurations > 0);
	if (num_labels == 0 || num_cuboids == 0) return;

	unsigned int num_cases = num_labels * num_axis_configurations;
	
	if (_add_dummy_label)
	{
		num_cases = num_cases + 1;
	}
	
	unsigned int mat_size = num_cuboids * num_cases;
	_potential_mat = Eigen::MatrixXd(mat_size, mat_size);
	_potential_mat.setZero();


	// NOTE:
	// The potential matrix should be symmetric.

	// Cuboid with various labels and axes configurations.
	std::vector< std::vector<MeshCuboid *> > cuboids_with_label_and_axes(num_cuboids);
	std::vector< std::vector<MeshCuboidAttributes> > attributes_with_label_and_axes(num_cuboids);
	std::vector< std::vector<MeshCuboidTransformation> > transformations_with_label_and_axes(num_cuboids);

	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
	{
		cuboids_with_label_and_axes[cuboid_index].resize(num_cases);
		attributes_with_label_and_axes[cuboid_index].resize(num_cases);
		transformations_with_label_and_axes[cuboid_index].resize(num_cases);

		for (int case_index = 0; case_index < num_cases; ++case_index)
		{
			if (_add_dummy_label && case_index >= (num_cases - 1))
			{
				continue;
			}

			unsigned int label_index = (case_index / num_axis_configurations);
			unsigned int axis_configuration_index =
				(case_index % num_axis_configurations);
			assert(label_index < num_labels);
			assert(axis_configuration_index < num_axis_configurations);

			// FIXME:
			// The cuboid should not deep copy all cuboid surface points.
			MeshCuboid *cuboid = new MeshCuboid(*_cuboids[cuboid_index]);	// Copy constructor.
			cuboid->set_label_index(label_index);
			cuboid->set_axis_configuration(axis_configuration_index);

			cuboids_with_label_and_axes[cuboid_index][case_index] = cuboid;
			attributes_with_label_and_axes[cuboid_index][case_index].compute_attributes(cuboid);
			transformations_with_label_and_axes[cuboid_index][case_index].compute_transformation(cuboid);
		}
	}


	// Single potentials.
	std::cout << "Unary potentials..." << std::endl;

	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
	{
		for (unsigned int case_index = 0; case_index < num_cases; ++case_index)
		{
			if (_add_dummy_label && case_index >= (num_cases - 1))
			{
				continue;
			}

			unsigned int mat_index = cuboid_index * num_cases + case_index;
			unsigned int label_index = (case_index / num_axis_configurations);
			unsigned int axis_configuration_index =
				(case_index % num_axis_configurations);
			assert(label_index < num_labels);
			assert(axis_configuration_index < num_axis_configurations);
			Label label = _labels[label_index];

			const MeshCuboid *cuboid = cuboids_with_label_and_axes[cuboid_index][case_index];
			const MeshCuboidAttributes &attributes = attributes_with_label_and_axes[cuboid_index][case_index];
			const MeshCuboidTransformation &transformation = transformations_with_label_and_axes[cuboid_index][case_index];

			//Real potential = _predictor.get_single_potential(cuboid, &attributes, &transformation, label_index);
			//assert(potential >= 0.0);

			Real potential = 0.0;
			if (_cuboids[cuboid_index]->get_label_index() != label_index)
				potential = FLAGS_param_max_potential;

			// NOTE:
			// If label symmetry information is given, symmetric labels have zero potential value.
			if (_label_symmetries)
			{
				assert((*_label_symmetries).size() == num_labels);
				const std::list<LabelIndex> &label_symmetry = (*_label_symmetries)[label_index];
				for (std::list<LabelIndex>::const_iterator it = label_symmetry.begin(); it != label_symmetry.end(); ++it)
				{
					if (_cuboids[cuboid_index]->get_label_index() == (*it))
					{
						potential = 0.0;
						break;
					}
				}
			}

			_potential_mat(mat_index, mat_index) = potential;
		}
	}

	// Pairwise potentials.
	std::cout << "Binary potentials..." << std::endl;

	for (unsigned int cuboid_index_1 = 0; cuboid_index_1 < num_cuboids - 1; ++cuboid_index_1)
	{
		for (unsigned int case_index_1 = 0; case_index_1 < num_cases; ++case_index_1)
		{
			if (_add_dummy_label && case_index_1 >= (num_cases - 1))
			{
				continue;
			}

			unsigned int mat_index_1 = cuboid_index_1 * num_cases + case_index_1;
			unsigned int label_index_1 = (case_index_1 / num_axis_configurations);
			unsigned int axis_configuration_index_1 =
				(case_index_1 % num_axis_configurations);
			assert(label_index_1 < num_labels);
			assert(axis_configuration_index_1 < num_axis_configurations);
			Label label_1 = _labels[label_index_1];

			const MeshCuboid *cuboid_1 = cuboids_with_label_and_axes[cuboid_index_1][case_index_1];
			const MeshCuboidAttributes &attributes_1 = attributes_with_label_and_axes[cuboid_index_1][case_index_1];
			const MeshCuboidTransformation &transformation_1 = transformations_with_label_and_axes[cuboid_index_1][case_index_1];

			for (unsigned int cuboid_index_2 = cuboid_index_1 + 1; cuboid_index_2 < num_cuboids; ++cuboid_index_2)
			{
				assert(cuboid_index_1 != cuboid_index_2);

#pragma omp parallel for
				for (int case_index_2 = 0; case_index_2 < num_cases; ++case_index_2)
				{
					if (_add_dummy_label && case_index_2 >= (num_cases - 1))
					{
						continue;
					}

					unsigned int mat_index_2 = cuboid_index_2 * num_cases + case_index_2;

					// Symmetric matrix. Diagonals are for single potentials.
					if (mat_index_1 >= mat_index_2)
						continue;

					unsigned int label_index_2 = (case_index_2 / num_axis_configurations);
					unsigned int axis_configuration_index_2 =
						(case_index_2 % num_axis_configurations);
					assert(label_index_2 < num_labels);
					assert(axis_configuration_index_2 < num_axis_configurations);
					Label label_2 = _labels[label_index_2];

					MeshCuboid *cuboid_2 = cuboids_with_label_and_axes[cuboid_index_2][case_index_2];
					const MeshCuboidAttributes &attributes_2 = attributes_with_label_and_axes[cuboid_index_2][case_index_2];
					const MeshCuboidTransformation &transformation_2 = transformations_with_label_and_axes[cuboid_index_2][case_index_2];

					if (label_index_1 == label_index_2)
					{
						// NOTE:
						// Currently, it is NOT allowed that multiple parts have the same label.
						_potential_mat(mat_index_1, mat_index_2) = FLAGS_param_max_potential;
						_potential_mat(mat_index_2, mat_index_1) = FLAGS_param_max_potential;
						continue;
					}
					assert(label_1 != label_2);

					Real potential = _predictor.get_pair_potential(
						cuboid_1, cuboid_2, &attributes_1, &attributes_2,
						&transformation_1, &transformation_2, label_index_1, label_index_2);
					assert(potential >= 0.0);

					// Symmetric matrix.
					_potential_mat(mat_index_1, mat_index_2) = potential;
					_potential_mat(mat_index_2, mat_index_1) = potential;
				}
			}
		}
	}


	// Set dummy label pair potentials.
	if (_add_dummy_label)
	{
		const unsigned int dummy_case_index = num_cases - 1;
		const Real dummy_potential = FLAGS_param_dummy_potential;

		for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
		{
			for (unsigned int case_index = 0; case_index < num_cases; ++case_index)
			{
				unsigned int mat_index = cuboid_index * num_cases + case_index;

				for (unsigned int dummy_cuboid_index = 0; dummy_cuboid_index < num_cuboids; ++dummy_cuboid_index)
				{
					if (cuboid_index == dummy_cuboid_index)
						continue;

					unsigned int dummy_mat_index = dummy_cuboid_index * num_cases + dummy_case_index;
					_potential_mat(mat_index, dummy_mat_index) = dummy_potential;
					_potential_mat(dummy_mat_index, mat_index) = dummy_potential;
				}
			}
		}
	}


	// Deallocate.
	for (std::vector< std::vector<MeshCuboid *> >::iterator it = cuboids_with_label_and_axes.begin();
		it != cuboids_with_label_and_axes.end(); ++it)
	{
		for (std::vector<MeshCuboid *>::iterator jt = (*it).begin(); jt != (*it).end(); ++jt)
		{
			delete (*jt);
		}
	}
}

void recognize_labels_and_axes_configurations(
	MeshCuboidStructure &_cuboid_structure,
	const MeshCuboidPredictor &_predictor,
	const std::string _log_filename,
	bool _use_symmetry_info,
	bool _add_dummy_label)
{
	std::ofstream log_file(_log_filename, std::ofstream::out | std::ofstream::app);
	assert(log_file);

	const std::vector<Label>& labels = _cuboid_structure.labels_;
	std::vector<MeshCuboid *> all_cuboids = _cuboid_structure.get_all_cuboids();

	unsigned int num_labels = labels.size();
	unsigned int num_cuboids = all_cuboids.size();
	unsigned int num_axis_configurations = MeshCuboid::num_axis_configurations();
	assert(num_axis_configurations > 0);
	if (num_labels == 0 || num_cuboids == 0) return;

	unsigned int num_cases = num_labels * num_axis_configurations;

	if (_add_dummy_label)
	{
		num_cases = num_cases + 1;
	}

	// NOTE:
	// This matrix should be symmetric.
	Eigen::MatrixXd potential_mat;

	if (_use_symmetry_info)
	{
		compute_labels_and_axes_configuration_potentials(
			labels, all_cuboids, _predictor, potential_mat, &_cuboid_structure.label_symmetries_, _add_dummy_label);
	}
	else
	{
		compute_labels_and_axes_configuration_potentials(
			labels, all_cuboids, _predictor, potential_mat, NULL, _add_dummy_label);
	}


	// Solve MRF.
	std::vector<int> output = solve_markov_random_field(num_cuboids, num_cases, potential_mat);
	assert(output.size() == num_cuboids);


	// Print results.
	Eigen::VectorXd part_case_pair_vec;

	part_case_pair_vec = Eigen::VectorXd::Zero(num_cuboids * num_cases);
	log_file << " -- Input -- " << std::endl;

	// TEST.
	bool *is_label_visited = new bool[num_labels];
	memset(is_label_visited, false, num_labels*sizeof(bool));

	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
	{
		MeshCuboid *cuboid = all_cuboids[cuboid_index];
		LabelIndex label_index = cuboid->get_label_index();

		// TEST.
		if (is_label_visited[label_index])
		{
			log_file << "[" << cuboid_index << "]: Dummy" << std::endl;
			unsigned int mat_index = cuboid_index * num_cases + (num_cases - 1);
			part_case_pair_vec(mat_index) = 1.0;
			continue;
		}
		is_label_visited[label_index] = true;

		Label label = labels[label_index];
		unsigned int axis_configuration_index = 0;
		unsigned int case_index = label_index * num_axis_configurations + axis_configuration_index;
		unsigned int mat_index = cuboid_index * num_cases + case_index;
		part_case_pair_vec(mat_index) = 1.0;

		log_file << "[" << cuboid_index << "]: " << label << ", " << axis_configuration_index << std::endl;
	}

	// TEST.
	delete [] is_label_visited;

	log_file << "Energy = " << part_case_pair_vec.transpose() * potential_mat * part_case_pair_vec;
	log_file << std::endl;


	part_case_pair_vec = Eigen::VectorXd::Zero(num_cuboids * num_cases);
	log_file << " -- Output -- " << std::endl;

	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
	{
		MeshCuboid *cuboid = all_cuboids[cuboid_index];
		assert(output[cuboid_index] >= 0);

		unsigned int case_index = static_cast<unsigned int>(output[cuboid_index]);

		// Dummy label.
		if (_add_dummy_label && case_index >= (num_cases - 1))
		{
			log_file << "[" << cuboid_index << "]: Dummy" << std::endl;
			unsigned int mat_index = cuboid_index * num_cases + (num_cases - 1);
			part_case_pair_vec(mat_index) = 1.0;

			delete cuboid;
			all_cuboids[cuboid_index] = NULL;
			continue;
		}

		unsigned int mat_index = cuboid_index * num_cases + case_index;
		part_case_pair_vec(mat_index) = 1.0;

		LabelIndex label_index = (case_index / num_axis_configurations);
		unsigned int axis_configuration_index =
			(case_index % num_axis_configurations);
		Label label = labels[label_index];

		cuboid->set_label_index(label_index);
		cuboid->set_axis_configuration(axis_configuration_index);

		log_file << "[" << cuboid_index << "]: " << label << ", " << axis_configuration_index << std::endl;
	}
	log_file << "Energy = " << part_case_pair_vec.transpose() * potential_mat * part_case_pair_vec;
	log_file << std::endl;

	log_file.close();


	// Put recognized cuboids.
	_cuboid_structure.label_cuboids_.clear();
	_cuboid_structure.label_cuboids_.resize(num_labels);
	for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
		_cuboid_structure.label_cuboids_[label_index].clear();

	for (std::vector<MeshCuboid *>::const_iterator it = all_cuboids.begin();
		it != all_cuboids.end(); ++it)
	{
		MeshCuboid *cuboid = (*it);

		// Dummy label.
		if (!cuboid) continue;

		LabelIndex label_index = cuboid->get_label_index();
		assert(label_index < _cuboid_structure.num_labels());
		_cuboid_structure.label_cuboids_[label_index].push_back(*it);
	}
}

void test_recognize_labels_and_axes_configurations(
	const std::vector<Label>& _labels,
	const std::vector<MeshCuboid *>& _cuboids,
	const MeshCuboidPredictor &_predictor,
	const std::string _log_filename)
{
	std::ofstream log_file(_log_filename, std::ofstream::out | std::ofstream::app);
	assert(log_file);

	unsigned int num_labels = _labels.size();
	unsigned int num_cuboids = _cuboids.size();
	if (num_labels == 0 || num_cuboids == 0) return;

	// NOTE:
	// This matrix should be symmetric.
	Eigen::MatrixXd potential_mat(num_labels, num_labels);
	potential_mat.setZero();


	// Single potentials.
	std::cout << "Unary potentials..." << std::endl;

	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
	{
		MeshCuboid *cuboid = _cuboids[cuboid_index];
		LabelIndex label_index = cuboid->get_label_index();
		MeshCuboidAttributes attributes;
		MeshCuboidTransformation transformation;
		attributes.compute_attributes(cuboid);
		transformation.compute_transformation(cuboid);

		Real potential = _predictor.get_single_potential(
			cuboid, &attributes, &transformation, label_index);
		assert(potential >= 0.0);

		assert(potential_mat(label_index, label_index) == 0);
		potential_mat(label_index, label_index) = potential;
	}

	// Pairwise potentials.
	std::cout << "Binary potentials..." << std::endl;

	for (unsigned int cuboid_index_1 = 0; cuboid_index_1 < num_cuboids - 1; ++cuboid_index_1)
	{
		MeshCuboid *cuboid_1 = _cuboids[cuboid_index_1];
		LabelIndex label_index_1 = cuboid_1->get_label_index();
		MeshCuboidAttributes attributes_1;
		MeshCuboidTransformation transformation_1;
		attributes_1.compute_attributes(cuboid_1);
		transformation_1.compute_transformation(cuboid_1);

		for (unsigned int cuboid_index_2 = 0; cuboid_index_2 < num_cuboids; ++cuboid_index_2)
		{
			if (cuboid_index_1 == cuboid_index_2)
				continue;

			MeshCuboid *cuboid_2 = _cuboids[cuboid_index_2];
			LabelIndex label_index_2 = cuboid_2->get_label_index();
			MeshCuboidAttributes attributes_2;
			MeshCuboidTransformation transformation_2;
			attributes_2.compute_attributes(cuboid_2);
			transformation_2.compute_transformation(cuboid_2);


			Real potential = _predictor.get_pair_potential(
				cuboid_1, cuboid_2, &attributes_1, &attributes_2,
				&transformation_1, &transformation_2, label_index_1, label_index_2);
			assert(potential >= 0.0);

			// Symmetric matrix.
			assert(potential_mat(label_index_1, label_index_2) == 0);
			assert(potential_mat(label_index_2, label_index_1) == 0);
			potential_mat(label_index_1, label_index_2) = potential;
			potential_mat(label_index_2, label_index_1) = potential;


#ifdef DEBUG_TEST
			const unsigned int num_attributes = MeshCuboidAttributes::k_num_attributes;
			const unsigned int mat_size = num_cuboids * num_attributes;
			Eigen::MatrixXd pair_quadratic_term(mat_size, mat_size);
			Eigen::VectorXd pair_linear_term(mat_size);
			double pair_constant_term;

			Real same_potential = 0;
			unsigned int start_index_1 = (cuboid_index_1 * num_attributes);
			unsigned int start_index_2 = (cuboid_index_2 * num_attributes);
			Eigen::VectorXd x1 = attributes_1.get_attributes();
			Eigen::VectorXd x2 = attributes_2.get_attributes();
			Eigen::VectorXd x = Eigen::VectorXd::Zero(mat_size);
			x.block<num_attributes, 1>(start_index_1, 0) = x1;
			x.block<num_attributes, 1>(start_index_2, 0) = x2;

			_predictor.get_pair_quadratic_form(cuboid_1, cuboid_2,
				cuboid_index_1, cuboid_index_2,
				label_index_1, label_index_2,
				pair_quadratic_term, pair_linear_term, pair_constant_term);
			same_potential += (x.transpose() * pair_quadratic_term * x);
			same_potential += (2 * pair_linear_term.transpose() * x);
			same_potential += pair_constant_term;

			CHECK_NUMERICAL_ERROR(__FUNCTION__, potential, same_potential);
#endif
		}
	}


	// Print results.
	Real sum_potential = 0;
	for (unsigned int label_index_1 = 0; label_index_1 < num_labels - 1; ++label_index_1)
	{
		for (unsigned int label_index_2 = label_index_1 + 1; label_index_2 < num_labels; ++label_index_2)
		{
			// NOTE:
			// Symmetric potential matrix.
			assert(potential_mat(label_index_1, label_index_2) == potential_mat(label_index_2, label_index_1));

			log_file << " - (" << label_index_1 << ", " << label_index_2 << "): "
				<< potential_mat(label_index_1, label_index_2) << std::endl;

			sum_potential += potential_mat(label_index_1, label_index_2);
		}
	}

	log_file << "Total: " << sum_potential << std::endl;
	log_file.close();
	//
}

void get_optimization_formulation(
	const std::vector<MeshCuboid *>& _cuboids,
	const MeshCuboidPredictor &_predictor,
	Eigen::VectorXd &_init_values,
	Eigen::MatrixXd &_single_quadratic_term, Eigen::MatrixXd &_pair_quadratic_term,
	Eigen::VectorXd &_single_linear_term, Eigen::VectorXd &_pair_linear_term,
	double &_single_constant_term, double &_pair_constant_term,
	double &_single_total_energy, double &_pair_total_energy)
{
	const unsigned int num_attributes = MeshCuboidAttributes::k_num_attributes;
	unsigned int num_cuboids = _cuboids.size();
	unsigned int mat_size = num_cuboids * num_attributes;

	_single_quadratic_term = Eigen::MatrixXd(mat_size, mat_size);
	_pair_quadratic_term = Eigen::MatrixXd(mat_size, mat_size);

	_single_linear_term = Eigen::VectorXd(mat_size);
	_pair_linear_term = Eigen::VectorXd(mat_size);

	_single_quadratic_term.setZero(); _pair_quadratic_term.setZero();
	_single_linear_term.setZero(); _pair_linear_term.setZero();
	_single_constant_term = 0; _pair_constant_term = 0;

	_init_values = Eigen::VectorXd(mat_size);

	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
	{
		MeshCuboid *cuboid = _cuboids[cuboid_index];
		MeshCuboidAttributes attributes;
		attributes.compute_attributes(cuboid);
		_init_values.segment<num_attributes>(num_attributes * cuboid_index)
			= attributes.get_attributes();
	}


	// Single energy (ICP prior energy).
	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
	{
		MeshCuboid *cuboid = _cuboids[cuboid_index];
		//LabelIndex label_index = cuboid->get_label_index();

		Eigen::MatrixXd each_single_quadratic_term(mat_size, mat_size);
		Eigen::VectorXd each_single_linear_term(mat_size);
		double each_single_constant_term;

		_predictor.get_single_quadratic_form(cuboid, cuboid_index,
			each_single_quadratic_term, each_single_linear_term, each_single_constant_term);

		_single_quadratic_term = _single_quadratic_term + each_single_quadratic_term;
		_single_linear_term = _single_linear_term + each_single_linear_term;
		_single_constant_term = _single_constant_term + each_single_constant_term;
	}

	_single_total_energy = 0;
	_single_total_energy += _init_values.transpose() * _single_quadratic_term * _init_values;
	_single_total_energy += 2 * _single_linear_term.transpose() * _init_values;
	_single_total_energy += _single_constant_term;


	// Pairwise energy.
	double same_pair_total_energy = 0;

	for (unsigned int cuboid_index_1 = 0; cuboid_index_1 < num_cuboids; ++cuboid_index_1)
	{
		MeshCuboid *cuboid_1 = _cuboids[cuboid_index_1];
		LabelIndex label_index_1 = cuboid_1->get_label_index();

		//for (unsigned int cuboid_index_2 = 0; cuboid_index_2 < num_cuboids; ++cuboid_index_2)
		//{
		//	if (cuboid_index_1 == cuboid_index_2)
		//		continue;

		// Using bilateral relations.
		// (A, B) and (B, A) pairs are the same.
		for (unsigned int cuboid_index_2 = cuboid_index_1 + 1; cuboid_index_2 < num_cuboids; ++cuboid_index_2)
		{
			MeshCuboid *cuboid_2 = _cuboids[cuboid_index_2];
			LabelIndex label_index_2 = cuboid_2->get_label_index();

			Eigen::MatrixXd each_pair_quadratic_term(mat_size, mat_size);
			Eigen::VectorXd each_pair_linear_term(mat_size);
			double each_pair_constant_term;

			Real energy = _predictor.get_pair_quadratic_form(cuboid_1, cuboid_2,
				cuboid_index_1, cuboid_index_2,
				label_index_1, label_index_2,
				each_pair_quadratic_term, each_pair_linear_term, each_pair_constant_term);

			_pair_quadratic_term = _pair_quadratic_term + each_pair_quadratic_term;
			_pair_linear_term = _pair_linear_term + each_pair_linear_term;
			_pair_constant_term = _pair_constant_term + each_pair_constant_term;

			same_pair_total_energy += energy;
		}
	}

	_pair_total_energy = 0;
	_pair_total_energy += _init_values.transpose() * _pair_quadratic_term * _init_values;
	_pair_total_energy += 2 * _pair_linear_term.transpose() * _init_values;
	_pair_total_energy += _pair_constant_term;

#ifdef DEBUG_TEST
	CHECK_NUMERICAL_ERROR(__FUNCTION__, _pair_total_energy, same_pair_total_energy);
#endif
}

void get_optimization_error(
	const std::vector<MeshCuboid *>& _cuboids,
	const MeshCuboidPredictor &_predictor,
	double &_single_total_energy, double &_pair_total_energy)
{
	Eigen::VectorXd init_values;
	Eigen::MatrixXd single_quadratic_term, pair_quadratic_term;
	Eigen::VectorXd single_linear_term, pair_linear_term;
	double single_constant_term, pair_constant_term;

	get_optimization_formulation(_cuboids, _predictor, init_values,
		single_quadratic_term, pair_quadratic_term,
		single_linear_term, pair_linear_term,
		single_constant_term, pair_constant_term,
		_single_total_energy, _pair_total_energy);
}

void optimize_attributes_once(
	const std::vector<MeshCuboid *>& _cuboids,
	const MeshCuboidPredictor& _predictor,
	const double _quadprog_ratio)
{
	const unsigned int num_attributes = MeshCuboidAttributes::k_num_attributes;
	unsigned int num_cuboids = _cuboids.size();
	unsigned int mat_size = num_cuboids * num_attributes;

	Eigen::VectorXd init_values;
	Eigen::MatrixXd single_quadratic_term, pair_quadratic_term;
	Eigen::VectorXd single_linear_term, pair_linear_term;
	double single_constant_term, pair_constant_term;
	double single_total_energy, pair_total_energy;

	get_optimization_formulation(_cuboids, _predictor, init_values,
		single_quadratic_term, pair_quadratic_term,
		single_linear_term, pair_linear_term,
		single_constant_term, pair_constant_term,
		single_total_energy, pair_total_energy);


	Eigen::MatrixXd quadratic_term = pair_quadratic_term + _quadprog_ratio * single_quadratic_term;
	Eigen::VectorXd linear_term = pair_linear_term + _quadprog_ratio * single_linear_term;
	double constant_term = pair_constant_term + _quadprog_ratio * single_constant_term;


	// Solve quadratic programming.
	Eigen::VectorXd new_values = solve_quadratic_programming(quadratic_term, linear_term, constant_term,
		&init_values, NULL, _quadprog_ratio);

	assert(new_values.rows() == mat_size);


	// Update cuboid.
	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
	{
		MeshCuboid *cuboid = _cuboids[cuboid_index];
		Eigen::VectorXd new_attributes_vec = new_values.segment(
			num_attributes * cuboid_index, num_attributes);

		MyMesh::Point new_bbox_center(0.0);
		std::array<MyMesh::Point, MeshCuboid::k_num_corners> new_bbox_corners;

		//for (unsigned int i = 0; i < 3; ++i)
		//	new_bbox_center[i] = new_attributes_vec[MeshCuboidAttributes::k_center_index + i];

		for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
		{
			for (unsigned int i = 0; i < 3; ++i)
			{
				new_bbox_corners[corner_index][i] = new_attributes_vec[
					MeshCuboidAttributes::k_corner_index + 3 * corner_index + i];
			}
			new_bbox_center += new_bbox_corners[corner_index];
		}
		new_bbox_center = new_bbox_center / MeshCuboid::k_num_corners;

		cuboid->set_bbox_center(new_bbox_center);
		cuboid->set_bbox_corners(new_bbox_corners);

		// Update cuboid surface points.
		for (unsigned int point_index = 0; point_index < cuboid->num_cuboid_surface_points();
			++point_index)
		{
			MeshCuboidSurfacePoint *cuboid_surface_point = cuboid->get_cuboid_surface_point(point_index);

			MyMesh::Point new_point(0.0);
			for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
				new_point += cuboid_surface_point->corner_weights_[corner_index] * new_bbox_corners[corner_index];

			cuboid_surface_point->point_ = new_point;
		}
	}
}

void optimize_attributes(
	MeshCuboidStructure &_cuboid_structure,
	const Real _modelview_matrix[16],
	const MeshCuboidPredictor &_predictor,
	const double _quadprog_ratio,
	const std::string _log_filename,
	const unsigned int _max_num_iterations,
	GLViewerCore *_viewer)
{
	std::ofstream log_file(_log_filename, std::ofstream::out | std::ofstream::app);
	assert(log_file);

	unsigned int num_labels = _cuboid_structure.num_labels();

	std::vector<MeshCuboid *> all_cuboids = _cuboid_structure.get_all_cuboids();
	unsigned int num_cuboids = all_cuboids.size();

	std::stringstream sstr;
	double single_total_energy, pair_total_energy, total_energy;


	// NOTE: Keep the lowest energy cuboids.
	std::vector<MeshCuboid *> final_all_cuboids;
	double final_total_energy = std::numeric_limits<double>::max();
	unsigned int final_cuboid_iteration = 0;

	final_all_cuboids.resize(num_cuboids);
	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
		final_all_cuboids[cuboid_index] = new MeshCuboid(*all_cuboids[cuboid_index]);

	update_cuboid_surface_points(_cuboid_structure, _modelview_matrix);
	if (_viewer) _viewer->updateGL();

	//
	get_optimization_error(all_cuboids, _predictor,
		single_total_energy, pair_total_energy);
	total_energy = pair_total_energy + _quadprog_ratio * single_total_energy;
	sstr.str(std::string());
	sstr << "Energy: (pair = " << pair_total_energy
		<< ", single = " << _quadprog_ratio * single_total_energy
		<< ", total = " << total_energy << ")" << std::endl;
	std::cout << sstr.str(); log_file << sstr.str();
	//

	std::cout << std::endl; log_file << std::endl;


	unsigned int iteration = 1;
	for (; iteration <= _max_num_iterations; ++iteration)
	{
		sstr.str(std::string());
		sstr << "iteration [" << iteration << "]" << std::endl;
		std::cout << sstr.str(); log_file << sstr.str();

		optimize_attributes_once(all_cuboids, _predictor, _quadprog_ratio);
		if (_viewer) _viewer->updateGL();

		//
		get_optimization_error(all_cuboids, _predictor,
			single_total_energy, pair_total_energy);
		total_energy = pair_total_energy + _quadprog_ratio * single_total_energy;
		sstr.str(std::string());
		sstr << " - After optimization" << std::endl;
		sstr << "Energy: (pair = " << pair_total_energy
			<< ", single = " << _quadprog_ratio * single_total_energy
			<< ", total = " << total_energy << ")" << std::endl;
		std::cout << sstr.str(); log_file << sstr.str();
		//

		// Cuboidize.
		for (std::vector<MeshCuboid *>::const_iterator it = all_cuboids.begin();
			it != all_cuboids.end(); ++it)
			(*it)->cuboidize();

		update_cuboid_surface_points(_cuboid_structure, _modelview_matrix);
		if (_viewer) _viewer->updateGL();

		//
		get_optimization_error(all_cuboids, _predictor,
			single_total_energy, pair_total_energy);
		total_energy = pair_total_energy + _quadprog_ratio * single_total_energy;
		sstr.str(std::string());
		sstr << " - After cuboidization" << std::endl;
		sstr << "Energy: (pair = " << pair_total_energy
			<< ", single = " << _quadprog_ratio * single_total_energy
			<< ", total = " << total_energy << ")";
		if (total_energy < final_total_energy)
			sstr << " - Minimum.";
		sstr << std::endl;

		std::cout << sstr.str(); log_file << sstr.str();
		//

		std::cout << std::endl; log_file << std::endl;
				
		if (total_energy < final_total_energy)
		{
			//
			final_total_energy = total_energy;
			final_cuboid_iteration = iteration;

			for (std::vector<MeshCuboid *>::iterator it = final_all_cuboids.begin();
				it != final_all_cuboids.end(); ++it)
				delete (*it);

			final_all_cuboids.resize(num_cuboids);
			for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
				final_all_cuboids[cuboid_index] = new MeshCuboid(*all_cuboids[cuboid_index]);
			//
		}
		else if (total_energy > 1.5 * final_total_energy)
		{
			sstr.str(std::string());
			sstr << "Last iteration is worse than the previous one..." << std::endl;
			sstr << "Recover to the previous iteration... Stop." << std::endl;
			std::cout << sstr.str(); log_file << sstr.str();

			break;
		}
	}

	if (iteration >= _max_num_iterations)
	{
		sstr.str(std::string());
		sstr << "# of iteration exceeds maximum number of iterations ("
			<< _max_num_iterations << ") ... Stop." << std::endl;
		std::cout << sstr.str(); log_file << sstr.str();
	}

	// Copy final cuboids.
	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
		(*all_cuboids[cuboid_index]) = (*final_all_cuboids[cuboid_index]);
	
	for (std::vector<MeshCuboid *>::iterator it = final_all_cuboids.begin();
		it != final_all_cuboids.end(); ++it)
		delete (*it);

	update_cuboid_surface_points(_cuboid_structure, _modelview_matrix);
	if (_viewer) _viewer->updateGL();

	//
	sstr.str(std::string());
	sstr << "Final cuboid iteration: " << final_cuboid_iteration << std::endl;
	std::cout << sstr.str(); log_file << sstr.str();

	get_optimization_error(all_cuboids, _predictor,
		single_total_energy, pair_total_energy);
	total_energy = pair_total_energy + _quadprog_ratio * single_total_energy;
	sstr.str(std::string());
	sstr << "Energy: (pair = " << pair_total_energy
		<< ", single = " << _quadprog_ratio * single_total_energy
		<< ", total = " << total_energy << ")" << std::endl;
	std::cout << sstr.str(); log_file << sstr.str();
	//

	log_file.close();
}

void add_missing_cuboids_once(
	const std::vector<MeshCuboid *>& _given_cuboids,
	const std::list<LabelIndex> &_missing_label_indices,
	const MeshCuboidPredictor &_predictor,
	std::vector<MeshCuboid *>& _new_cuboids)
{
	if (_given_cuboids.empty() || _missing_label_indices.empty())
		return;


	unsigned int num_given_cuboids = _given_cuboids.size();
	std::vector<MeshCuboid *> all_cuboids;
	all_cuboids.insert(all_cuboids.end(), _given_cuboids.begin(), _given_cuboids.end());

	for (std::list<LabelIndex>::const_iterator it = _missing_label_indices.begin();
		it != _missing_label_indices.end(); ++it)
	{
		LabelIndex label_index = (*it);
		MeshCuboid *missing_cuboid = new MeshCuboid(label_index);
		all_cuboids.push_back(missing_cuboid);
	}

	unsigned int num_cuboids = all_cuboids.size();
	unsigned int num_new_cuboids = num_cuboids - num_given_cuboids;
	assert(num_new_cuboids > 0);

	const unsigned int num_attributes = MeshCuboidAttributes::k_num_attributes;
	unsigned int mat_size = num_cuboids * num_attributes;

	Eigen::MatrixXd quadratic_term(mat_size, mat_size);
	Eigen::VectorXd linear_term(mat_size);
	double constant_term = 0;
	quadratic_term.setZero();
	linear_term.setZero();

	for (unsigned int cuboid_index_1 = 0; cuboid_index_1 < num_cuboids; ++cuboid_index_1)
	{
		MeshCuboid *cuboid_1 = all_cuboids[cuboid_index_1];
		assert(cuboid_1);
		LabelIndex label_index_1 = cuboid_1->get_label_index();

		//for (unsigned int cuboid_index_2 = 0; cuboid_index_2 < num_cuboids; ++cuboid_index_2)
		//{
		//	if (cuboid_index_1 == cuboid_index_2) continue;

		// Using bilateral relations.
		// (A, B) and (B, A) pairs are the same.
		for (unsigned int cuboid_index_2 = cuboid_index_1 + 1; cuboid_index_2 < num_cuboids; ++cuboid_index_2)
		{
			MeshCuboid *cuboid_2 = all_cuboids[cuboid_index_2];
			assert(cuboid_2);
			LabelIndex label_index_2 = cuboid_2->get_label_index();

			Eigen::MatrixXd pair_quadratic_term(mat_size, mat_size);
			Eigen::VectorXd pair_linear_term(mat_size);
			double pair_constant_term;

			_predictor.get_pair_quadratic_form(cuboid_1, cuboid_2,
				cuboid_index_1, cuboid_index_2,
				label_index_1, label_index_2,
				pair_quadratic_term, pair_linear_term, pair_constant_term);

			quadratic_term = quadratic_term + pair_quadratic_term;
			linear_term = linear_term + pair_linear_term;
			constant_term = constant_term + pair_constant_term;
		}
	}


	// Attributes of given cuboids are fixed.
	unsigned int given_mat_size = num_given_cuboids * num_attributes;
	unsigned int new_mat_size = num_new_cuboids * num_attributes;
	assert(given_mat_size + new_mat_size == mat_size);

	Eigen::MatrixXd new_quadratic_term(new_mat_size, new_mat_size);
	Eigen::VectorXd new_linear_term(new_mat_size);
	double new_constant_term = 0;
	new_quadratic_term.setZero();
	new_linear_term.setZero();

	Eigen::VectorXd given_init_values = Eigen::VectorXd::Zero(given_mat_size);
	for (unsigned int cuboid_index = 0; cuboid_index < num_given_cuboids; ++cuboid_index)
	{
		MeshCuboid *cuboid = all_cuboids[cuboid_index];
		assert(cuboid);
		MeshCuboidAttributes attributes;
		attributes.compute_attributes(cuboid);
		given_init_values.segment<num_attributes>(num_attributes * cuboid_index) = attributes.get_attributes();
	}

	// NOTE:
	// First 'num_given_cuboids' cuboids are given cuboids,
	// and last 'num_new_cuboids' cuboids are missing cuboids.
	Eigen::MatrixXd quadratic_term_11 = quadratic_term.block(0, 0, given_mat_size, given_mat_size);
	Eigen::MatrixXd quadratic_term_12 = quadratic_term.block(0, given_mat_size, given_mat_size, new_mat_size);
	Eigen::MatrixXd quadratic_term_21 = quadratic_term.block(given_mat_size, 0, new_mat_size, given_mat_size);
	Eigen::MatrixXd quadratic_term_22 = quadratic_term.block(given_mat_size, given_mat_size, new_mat_size, new_mat_size);

	Eigen::VectorXd linear_tern_1 = linear_term.segment(0, given_mat_size);
	Eigen::VectorXd linear_tern_2 = linear_term.segment(given_mat_size, new_mat_size);

	new_quadratic_term = quadratic_term_22;

	new_linear_term += (0.5 * given_init_values.transpose() * quadratic_term_12).transpose();
	new_linear_term += (0.5 * given_init_values.transpose() * quadratic_term_21.transpose()).transpose();
	new_linear_term += linear_tern_2;

	new_constant_term += (given_init_values.transpose() * quadratic_term_11 * given_init_values);
	new_constant_term += (2 * linear_tern_1.transpose() * given_init_values);
	new_constant_term += constant_term;


	// Solve quadratic programming.
	Eigen::VectorXd new_values = solve_quadratic_programming(new_quadratic_term, new_linear_term, new_constant_term);

	assert(new_values.rows() == new_mat_size);


	_new_cuboids.clear();

	for (unsigned int new_cuboid_index = 0; new_cuboid_index < num_new_cuboids; ++new_cuboid_index)
	{
		unsigned int cuboid_index = num_given_cuboids + new_cuboid_index;
		MeshCuboid *cuboid = all_cuboids[cuboid_index];
		// NOTE:
		// new_cuboid_index.
		Eigen::VectorXd new_attributes_vec = new_values.segment<num_attributes>(
			num_attributes * new_cuboid_index);

		MyMesh::Point new_bbox_center(0.0);
		std::array<MyMesh::Point, MeshCuboid::k_num_corners> new_bbox_corners;

		//for (unsigned int i = 0; i < 3; ++i)
		//	new_bbox_center[i] = new_attributes_vec[MeshCuboidAttributes::k_center_index + i];
		for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
		{
			for (unsigned int i = 0; i < 3; ++i)
			{
				new_bbox_corners[corner_index][i] = new_attributes_vec[
					MeshCuboidAttributes::k_corner_index + 3 * corner_index + i];
			}
			new_bbox_center += new_bbox_corners[corner_index];
		}
		new_bbox_center = new_bbox_center / MeshCuboid::k_num_corners;

		cuboid->set_bbox_center(new_bbox_center);
		cuboid->set_bbox_corners(new_bbox_corners);

		// Cuboidize.
		cuboid->cuboidize();

		_new_cuboids.push_back(cuboid);
	}
}

void add_missing_cuboids(
	MeshCuboidStructure &_cuboid_structure,
	const Real _modelview_matrix[16],
	const std::list<LabelIndex> &_missing_label_indices,
	const MeshCuboidPredictor &_predictor)
{
	if (_missing_label_indices.empty())
		return;

	const Real radius = FLAGS_param_observed_point_radius * _cuboid_structure.mesh_->get_object_diameter();

	// NOTE:
	// Each existing cuboid creates candidates of missing cuboids.
	// The best missing cuboid for the same label will be found in the next iteration.
	std::vector<MeshCuboid *> all_cuboids = _cuboid_structure.get_all_cuboids();
	unsigned int num_all_cuboids = all_cuboids.size();
	for (unsigned int cuboid_index = 0; cuboid_index < num_all_cuboids; ++cuboid_index)
	{
		std::vector<MeshCuboid *> given_cuboids;
		given_cuboids.push_back(all_cuboids[cuboid_index]);
		std::vector<MeshCuboid *> new_cuboids;

		add_missing_cuboids_once(given_cuboids, _missing_label_indices, _predictor, new_cuboids);

		for (std::vector<MeshCuboid *>::iterator it = new_cuboids.begin(); it != new_cuboids.end(); ++it)
		{
			MeshCuboid *cuboid = (*it);
			assert(cuboid);

			cuboid->create_grid_points_on_cuboid_surface(
				FLAGS_param_num_cuboid_surface_points);

			cuboid->compute_cuboid_surface_point_visibility(
				_modelview_matrix, radius, _cuboid_structure.sample_points_);

			//if (cuboid->get_cuboid_overvall_visibility() > FLAGS_param_min_cuboid_overall_visibility)
			//{
			//	// The cuboid is placed in the visible area.
			//	delete cuboid;
			//}
			//else
			//{
				LabelIndex label_index = cuboid->get_label_index();
				_cuboid_structure.label_cuboids_[label_index].push_back(cuboid);
			//}
		}
	}
}

void symmetrize_cuboids(MeshCuboidStructure &_cuboid_structure)
{
	unsigned int num_given_labels = _cuboid_structure.num_labels();

	for (LabelIndex label_index_1 = 0; label_index_1 < num_given_labels; ++label_index_1)
	{
		//if (_cuboid_structure.label_symmetries_[label_index_1].empty())
		//	continue;

		Label label_1 = _cuboid_structure.labels_[label_index_1];
		MeshCuboid *cuboid_1 = NULL;

		// NOTE:
		// The current implementation assumes that there is only one part for each label.
		assert(_cuboid_structure.label_cuboids_[label_index_1].size() <= 1);
		if (!_cuboid_structure.label_cuboids_[label_index_1].empty())
			cuboid_1 = _cuboid_structure.label_cuboids_[label_index_1].front();


		// FIXME.
		// Assume that all cuboids have either intra-cuboid or inter-cuboid symmetry.
		if (cuboid_1 && _cuboid_structure.label_symmetries_[label_index_1].empty())
		{
			const int flip_axis_index = FLAGS_param_intra_cuboid_symmetry_axis;

			Eigen::Vector3d reflection_plane_point(0.0, 0.0, 0.0);
			Eigen::Vector3d reflection_plane_normal(0.0, 0.0, 0.0);

			for (unsigned int i = 0; i < 3; ++i)
			{
				reflection_plane_point[i] = cuboid_1->get_bbox_center()[i];
				reflection_plane_normal[i] = cuboid_1->get_bbox_axis(flip_axis_index)[i];
			}

			Eigen::Matrix3d R = Eigen::Matrix3d::Identity()
				- 2 * (reflection_plane_normal * reflection_plane_normal.transpose());
			Eigen::Vector3d t = 2 * (reflection_plane_normal * reflection_plane_normal.transpose())
				* reflection_plane_point;

			CHECK_NUMERICAL_ERROR(__FUNCTION__, (R * R - Eigen::Matrix3d::Identity()).norm());
			CHECK_NUMERICAL_ERROR(__FUNCTION__, (t + R * t).norm());


			std::vector<MeshSamplePoint *>cuboid_sample_points_1;
			cuboid_sample_points_1.insert(cuboid_sample_points_1.end(),
				cuboid_1->get_sample_points().begin(), cuboid_1->get_sample_points().end());

			for (std::vector<MeshSamplePoint *>::const_iterator it = cuboid_sample_points_1.begin();
				it != cuboid_sample_points_1.end(); ++it)
			{
				MeshSamplePoint *new_sample_point = new MeshSamplePoint(**it);

				Eigen::Vector3d p(0.0, 0.0, 0.0);
				for (unsigned int i = 0; i < 3; ++i)
					p[i] = new_sample_point->point_[i];
				Eigen::Vector3d transformed_p = R * p + t;
				for (unsigned int i = 0; i < 3; ++i)
					new_sample_point->point_[i] = transformed_p[i];

				cuboid_1->add_sample_point(new_sample_point);
				_cuboid_structure.sample_points_.push_back(new_sample_point);
			}
		}


		for (std::list<LabelIndex>::const_iterator it = _cuboid_structure.label_symmetries_[label_index_1].begin();
			it != _cuboid_structure.label_symmetries_[label_index_1].end(); ++it)
		{
			LabelIndex label_index_2 = (*it);
			assert(label_index_2 < num_given_labels);

			if (label_index_1 > label_index_2)
				continue;

			Label label_2 = _cuboid_structure.labels_[label_index_2];
			MeshCuboid *cuboid_2 = NULL;
				
			// NOTE:
			// The current implementation assumes that there is only one part for each label.
			assert(_cuboid_structure.label_cuboids_[label_index_2].size() <= 1);
			if (!_cuboid_structure.label_cuboids_[label_index_2].empty())
				cuboid_2 = _cuboid_structure.label_cuboids_[label_index_2].front();

			if (!cuboid_1 || !cuboid_2)
				continue;


			// NOTE:
			// Assume that all symmetric cases are reflectional symmetry.

			// Find the cuboid flipping axis.
			// Use heuristic. Just check the distance between corresponding cuboid corners.
			Eigen::Vector3d sum_axis_distances(0.0);
			for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
			{
				MyMesh::Normal diff = cuboid_1->get_bbox_corner(corner_index) - cuboid_2->get_bbox_corner(corner_index);
				for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
					sum_axis_distances[axis_index] += std::abs(diff[axis_index]);
			}

			int flip_axis_index = 0;
			sum_axis_distances.maxCoeff(&flip_axis_index);

			// Temporarily flip the axis of one cuboid.
			cuboid_2->flip_axis(flip_axis_index);

			
			// Find the accurate reflection plane,
			// which might not be perfectly perpendicular with one of x, y, or z axis.
			Eigen::Vector3d reflection_plane_point(0.0, 0.0, 0.0);
			Eigen::Vector3d reflection_plane_normal(0.0, 0.0, 0.0);

			MeshCuboid *group_cuboid = NULL;

			int group_label_index = _cuboid_structure.find_parent_label_index(label_index_1, label_index_2);
			if (group_label_index >= 0)
			{
				// NOTE:
				// The current implementation assumes that there is only one part for each label.
				assert(_cuboid_structure.label_cuboids_[group_label_index].size() <= 1);
				if (!_cuboid_structure.label_cuboids_[group_label_index].empty())
					group_cuboid = _cuboid_structure.label_cuboids_[group_label_index].front();
			}

			if (group_cuboid)
			{
				printf(">> Group Cuboid Index  = %d\n", group_label_index);
				for (unsigned int i = 0; i < 3; ++i)
				{
					reflection_plane_point[i] = group_cuboid->get_bbox_center()[i];
					reflection_plane_normal[i] = group_cuboid->get_bbox_axis(flip_axis_index)[i];
				}
			}
			else
			{
				for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
				{
					MyMesh::Point p = cuboid_1->get_bbox_corner(corner_index) + cuboid_2->get_bbox_corner(corner_index);
					for (unsigned int i = 0; i < 3; ++i)
						reflection_plane_point[i] += p[i];
				}
				reflection_plane_point /= static_cast<Real>(MeshCuboid::k_num_corners);

				// FIXME.
				//printf(">> reflection_plane_point[0] = %lf\n", reflection_plane_point[0]);
				//reflection_plane_point[0] = 0.0;

				Eigen::MatrixXd A(MeshCuboid::k_num_corners, 3);
				for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
				{
					MyMesh::Point p = cuboid_1->get_bbox_corner(corner_index) - cuboid_2->get_bbox_corner(corner_index);
					for (unsigned int i = 0; i < 3; ++i)
						A.row(corner_index)[i] = p[i];
				}

				// The first column of matrix V in SVD is the vector maximizing
				// cos(angle) with all columns in the given matrix A.
				Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
				reflection_plane_normal = svd.matrixV().col(0);

				// FIXME.
				//printf(">> reflection_plane_normal = %lf, %lf, %lf\n",
				//	reflection_plane_normal[0], reflection_plane_normal[1], reflection_plane_normal[2]);
				//reflection_plane_normal = Eigen::Vector3d::UnitX();
			}
			reflection_plane_normal.normalize();
			printf(">> reflection_plane_point[0] = %lf\n", reflection_plane_point[0]);
			printf(">> reflection_plane_normal = %lf, %lf, %lf\n",
				reflection_plane_normal[0], reflection_plane_normal[1], reflection_plane_normal[2]);


			// Define symmetry transformation matrix (' = transpose, n'n = 1).
			// p = reflection_plane_point
			// n = reflection_plane_normal
			// x - 2{n'(x - p)}n
			// = x - 2n{n'(x - p)}
			// = x - 2nn'x + 2nn'p
			// = (I - 2nn')x + 2nn'p
			// R = (I - 2nn'), t = 2nn'p.

			// The inverse transformation is the same.
			// RR = (I - 2nn')(I - 2nn') = I - 4nn' + 4nn'nn'
			// = I - 4nn' + 4nn' = I (since n'n = 1). R^(-1) = R.
			// -R^(-1)t = -Rt = -(I - 2nn')2nn'p = -2nn'p + 4nn'nn'p
			// = -2nn'p + 4nn'p = 2nn'p = t (since n'n = 1).

			Eigen::Matrix3d R = Eigen::Matrix3d::Identity()
				- 2 * (reflection_plane_normal * reflection_plane_normal.transpose());
			Eigen::Vector3d t = 2 * (reflection_plane_normal * reflection_plane_normal.transpose())
				* reflection_plane_point;

			CHECK_NUMERICAL_ERROR(__FUNCTION__, (R * R - Eigen::Matrix3d::Identity()).norm());
			CHECK_NUMERICAL_ERROR(__FUNCTION__, (t + R * t).norm());


			// Fit cuboid corners to become perfect symmetry.
			MyMesh::Point new_bbox_center_1(0.0);
			MyMesh::Point new_bbox_center_2(0.0);
			std::array<MyMesh::Point, MeshCuboid::k_num_corners> new_bbox_corners_1;
			std::array<MyMesh::Point, MeshCuboid::k_num_corners> new_bbox_corners_2;
			for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
			{
				Eigen::Vector3d p1(0.0), p2(0.0);
				for (unsigned int i = 0; i < 3; ++i)
				{
					p1[i] = cuboid_1->get_bbox_corner(corner_index)[i];
					p2[i] = cuboid_2->get_bbox_corner(corner_index)[i];
				}

				Eigen::Vector3d transformed_p2 = R * p2 + t;

				Eigen::Vector3d average_p1 = 0.5 * (p1 + transformed_p2);
				// The inverse transformation is the same.
				Eigen::Vector3d average_p2 = R * average_p1 + t;

				for (unsigned int i = 0; i < 3; ++i)
				{
					new_bbox_center_1[i] += average_p1[i];
					new_bbox_center_2[i] += average_p2[i];
					new_bbox_corners_1[corner_index][i] = average_p1[i];
					new_bbox_corners_2[corner_index][i] = average_p2[i];
				}
			}

			new_bbox_center_1 /= static_cast<Real>(MeshCuboid::k_num_corners);
			new_bbox_center_2 /= static_cast<Real>(MeshCuboid::k_num_corners);

			cuboid_1->set_bbox_center(new_bbox_center_1);
			cuboid_1->set_bbox_corners(new_bbox_corners_1);
			cuboid_1->cuboidize();

			cuboid_2->set_bbox_center(new_bbox_center_2);
			cuboid_2->set_bbox_corners(new_bbox_corners_2);
			cuboid_2->cuboidize();

			
			// Copy sample points.
			// 1 -> 2.
			std::vector<MeshSamplePoint *>cuboid_sample_points_1, cuboid_sample_points_2;
			cuboid_sample_points_1.insert(cuboid_sample_points_1.end(),
				cuboid_1->get_sample_points().begin(), cuboid_1->get_sample_points().end());
			cuboid_sample_points_2.insert(cuboid_sample_points_2.end(),
				cuboid_2->get_sample_points().begin(), cuboid_2->get_sample_points().end());

			for (std::vector<MeshSamplePoint *>::const_iterator it = cuboid_sample_points_1.begin();
				it != cuboid_sample_points_1.end(); ++it)
			{
				MeshSamplePoint *new_sample_point = new MeshSamplePoint(**it);

				Eigen::Vector3d p(0.0, 0.0, 0.0);
				for (unsigned int i = 0; i < 3; ++i)
					p[i] = new_sample_point->point_[i];
				Eigen::Vector3d transformed_p = R * p + t;
				for (unsigned int i = 0; i < 3; ++i)
					new_sample_point->point_[i] = transformed_p[i];

				cuboid_2->add_sample_point(new_sample_point);
				_cuboid_structure.sample_points_.push_back(new_sample_point);
			}

			// 2 -> 1.
			for (std::vector<MeshSamplePoint *>::const_iterator it = cuboid_sample_points_2.begin();
				it != cuboid_sample_points_2.end(); ++it)
			{
				MeshSamplePoint *new_sample_point = new MeshSamplePoint(**it);

				Eigen::Vector3d p(0.0, 0.0, 0.0);
				for (unsigned int i = 0; i < 3; ++i)
					p[i] = new_sample_point->point_[i];
				Eigen::Vector3d transformed_p = R * p + t;
				for (unsigned int i = 0; i < 3; ++i)
					new_sample_point->point_[i] = transformed_p[i];

				cuboid_1->add_sample_point(new_sample_point);
				_cuboid_structure.sample_points_.push_back(new_sample_point);
			}


			// Recover the flipped axis.
			cuboid_2->flip_axis(flip_axis_index);
		}
	}
}

/*
MeshCuboid *test_joint_normal_training(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const LabelIndex _label_index_1, const LabelIndex _label_index_2,
	const std::vector< std::vector<MeshCuboidJointNormalRelations> >& _relations)
{
	MeshCuboid *cuboid_2_copy = new MeshCuboid(*_cuboid_2);

	const MeshCuboidJointNormalRelations &relation_12 = _relations[_label_index_1][_label_index_2];
	Eigen::VectorXd mean = relation_12.get_mean();
	Eigen::VectorXd mean_21_vec = mean.segment(MeshCuboidFeatures::k_num_features, MeshCuboidFeatures::k_num_features);

	MeshCuboidTransformation transformation_1;
	transformation_1.compute_transformation(_cuboid_1);
	Eigen::Matrix3d inverse_rotation_1;
	Eigen::Vector3d inverse_translation_1;
	transformation_1.get_inverse_transformation(inverse_rotation_1, inverse_translation_1);

	Eigen::VectorXd inverse_transformed_mean_21_vec = mean_21_vec;
	for (unsigned int i = 0; i < MeshCuboidFeatures::k_num_local_points; ++i)
	{
		Eigen::VectorXd sub_values = mean_21_vec.block(3 * i, 0, 3, 1);
		sub_values = (inverse_rotation_1 * sub_values) + inverse_translation_1;
		inverse_transformed_mean_21_vec.block(3 * i, 0, 3, 1) = sub_values;
	}

	// FIXME:
	// Adhoc. Get attributes from features.
	Eigen::VectorXd new_attributes_vec_2 = inverse_transformed_mean_21_vec.segment(
		3, MeshCuboidAttributes::k_num_attributes);

	MyMesh::Point new_bbox_center(0.0);
	std::array<MyMesh::Point, MeshCuboid::k_num_corners> new_bbox_corners;

	//for (unsigned int i = 0; i < 3; ++i)
	//	new_bbox_center[i] = inverse_transformed_mean_21_vec[MeshCuboidAttributes::k_center_index + i];

	for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
	{
		for (unsigned int i = 0; i < 3; ++i)
		{
			new_bbox_corners[corner_index][i] = new_attributes_vec_2[
				MeshCuboidAttributes::k_corner_index + 3 * corner_index + i];
		}
		new_bbox_center += new_bbox_corners[corner_index];
	}
	new_bbox_center = new_bbox_center / MeshCuboid::k_num_corners;

	cuboid_2_copy->set_bbox_center_and_corners(new_bbox_center, new_bbox_corners);
	return cuboid_2_copy;
}
*/