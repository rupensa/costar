#ifndef MULTICOCLUST_H
#define MULTICOCLUST_H

#include <iostream>
#include <map>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <sstream>
#include <vector>
#include "Matrix.h"
#include "randomc.h"

class MultiViewCoClust {
private:
    int num_views_;

    // totals and computation related to to the dataset matrices
    double *n_tot_;
    double *two_divided_by_n_tot_;
    double *n_tot_square_;

    // constant to initialize in the algorithm
    // One I and C for each of the view
    double *Ic_;
    double *Cc_;

    double Ir_;
    double Cr_;

    // number of row cluster
    int n_row_clusters_;

    // number of column cluster x view
    int *n_col_clusters_;

    // assignment of the row points to each row cluster
    int *row_assignment_;

    // assignment of the column points to each column cluster x view
    // a matrix |# of view| * |column|
    int **col_assignment_;

    // actual distribution with row and column cluster x view
    // the first index is an index on a specific view
    // for each view we have a matrix that represents the actual distribution with
    // row and column cluster
    double ***T_;

    // original set of matrix with dimension
    std::vector<Matrix> dataset_;

    // the object to generate random number
    CRandomMersenne *random_generator_;

    // min number of row_cluster
    int min_n_row_cluster_;

    // output file name
    std::string output_filename_;
    
    void discreteInitialization();

    void randomInitialization();

    int randomInitCol(int view_index);

    int randomInitRow();

    void updateIandC();
    
    void rowPartitionCluster();
    
    void columnPartitionCluster(int view_index);

    double deltaTauR(double **sum_per_col, double *sum_lambdas, double **lambdas, 
                int source_cluster,
                int destination_cluster);
    
    double deltaTauC(int view_index, double *sum_col, double sum_lambdas,
                double *lambdas, int source_cluster, int destination_cluster);

    void modifyRowCluster(double **lambdas, int source_cluster, int destination_cluster);
    
    void modifyColumnCluster(int view_index, double *lambdas, int source_cluster, int destination_cluster);

    int checkEmptyRowCluster(int row_cluster);
    
    int checkEmptyColCluster(int view_index, int col_cluster);

    int chooseRandomElementFromCluster(int cluster, int *cluster_assignment, int assignment_size);

    std::string getTauRTauCString();

    void computeT(int view_index, int old_n_row_clusters, int new_n_row_clusters, int new_n_col_clusters);

    int checkRowClustConstraint(int row_cluster);

    int checkColClustConstraint(int view_index, int col_cluster);

    int findNonDominatedRow(std::vector<int> &row_indexes, double **lambdas,
            int source_cluster);

    double *computeEmulatedTauC(double **lambdas, int source_cluster, int destination_cluster);

    double *getTauC(double ***t, int temp_n_row_clusters);

    int findNonDominatedCol(std::vector<int> &col_indexes, double *lambdas,
            int source_cluster, int view_index);

    double computeEmulatedTauR(double *lambdas, int source_cluster, int destination_cluster, int view_index);

    double getTauR(double **t, int temp_n_col_clusters, int view_index);
    
public:

    /**
     * Constructor method
     *
     * @param dataset_matrices : pointer to the data matrices vector.
     * 		 The matrix contains one item per view, each item is a matrix
     * 		 with shape = (n_documents, n_features).
     *
     * @param min_num_row_cluster : int, minimum number of row clusters.
     */
    MultiViewCoClust(std::vector<Matrix> &dataset_matrices, int min_num_row_cluster) {
        num_views_ = dataset_matrices.size();
        dataset_ = dataset_matrices;
        Ic_ = new double[num_views_];
        Cc_ = new double[num_views_];

        min_n_row_cluster_ = min_num_row_cluster;
        random_generator_ = new CRandomMersenne(time(NULL));
        n_row_clusters_ = dataset_[0].n_row();
        n_col_clusters_ = new int[num_views_];
        for (int i = 0; i < num_views_; ++i) n_col_clusters_[i] = dataset_[i].n_col();

        row_assignment_ = new int[dataset_[0].n_row()];
        col_assignment_ = new int *[dataset_.size()];
        for (int i = 0; i < num_views_; ++i)
            col_assignment_[i] = new int[dataset_[i].n_col()];

        // nTot=nTot_pow=two_divideby_nTot=0;
        T_ = new double **[dataset_.size()];
        n_tot_ = new double[dataset_.size()];
        n_tot_square_ = new double[dataset_.size()];
        two_divided_by_n_tot_ = new double[dataset_.size()];
        for (int k = 0; k < num_views_; ++k) {
            T_[k] = dataset_[k].cloneMatrix();
            for (int i = 0; i < dataset_[k].n_row(); i++) {
                for (int j = 0; j < dataset_[k].n_col(); j++) n_tot_[k] += T_[k][i][j];
            }
            two_divided_by_n_tot_[k] = 2 / n_tot_[k];
            n_tot_square_[k] = pow(n_tot_[k], 2);
        }
        // initIandC();
    }

    /**
     * Constructor method
     *
     * @param dataset_matrices : pointer to the data matrices vector.
     * 			The matrix contains one item per view, each item is a matrix
     * 			with shape = (n_documents, n_features).
     *
     * @param min_num_row_cluster : int, minimum number of row clusters.
     *
     * @param output_files_prefix : string, the prefix to use for the output files' names.
     */
    MultiViewCoClust(std::vector<Matrix> &dataset_matrices, int min_num_row_cluster,
            std::string output_files_prefix) {
        min_n_row_cluster_ = min_num_row_cluster;
        num_views_ = dataset_matrices.size();
        dataset_ = dataset_matrices;

        Ic_ = new double[num_views_];
        Cc_ = new double[num_views_];

        min_n_row_cluster_ = min_num_row_cluster;

        // generate a 13 digits hexadecimal random number
        random_generator_ = new CRandomMersenne(time(NULL));

        // init nCoR to the number of documents
        n_row_clusters_ = dataset_[0].n_row();

        // init nCoC to the number of features for each view 
        n_col_clusters_ = new int[num_views_];
        for (int i = 0; i < num_views_; ++i) n_col_clusters_[i] = dataset_[i].n_col();

        // prepare lists for rows and cols assignments
        row_assignment_ = new int[dataset_[0].n_row()];
        col_assignment_ = new int *[dataset_.size()];
        for (int i = 0; i < num_views_; ++i)
            col_assignment_[i] = new int[dataset_[i].n_col()];

        // T is the clone of matrices 
        T_ = new double **[dataset_.size()];
        // n_tot is an array that contains for each matrix in T the sum of the elements 
        n_tot_ = new double[dataset_.size()];
        // n_tot_square_ is the element-wise power of nTot
        n_tot_square_ = new double[dataset_.size()];
        // two_divided_by_n_tot contains the value (2/e) for each e in n_tot
        two_divided_by_n_tot_ = new double[dataset_.size()];
        for (int view = 0; view < num_views_; ++view) {
            T_[view] = dataset_[view].cloneMatrix();
            for (int di = 0; di < dataset_[view].n_row(); di++) {
                for (int fi = 0; fi < dataset_[view].n_col(); fi++) {
                    n_tot_[view] += T_[view][di][fi];
                }
            }
            two_divided_by_n_tot_[view] = 2 / n_tot_[view];
            n_tot_square_[view] = pow(n_tot_[view], 2);
        }

        // save the output file prefix 
        output_filename_ = output_files_prefix;
    }

    void buildCoClu(int use_random_init, int n_iterations);

    void printDataset() {
        for (int k = 0; k < num_views_; ++k) {
            std::cout << "=============================================="
                    << std::endl;
            std::cout << "stampo matrice " << k << std::endl;
            dataset_[k].print();
            std::cout << "=============================================="
                    << std::endl;
        }
    }

    void printAssignment();

    void printContingencyT();

    void printView(int view_index);
    
    /**
     * Return a clone of the row assignment 
     * 
     * @return the row assignment
     */
    int *getRowAssignment() {
        int *ris = new int[dataset_[0].n_row()];
        for (int i = 0; i < dataset_[0].n_row(); i++) ris[i] = row_assignment_[i];
        return ris;
    }

    /**
     * Return a clone of the column assignment array for the specified view
     *
     * @param view_index, int, the considered view
     * @return the column assignment for the specified view
     */
    int *getColAssignment(int view_index) {
        int *ris = new int[dataset_[view_index].n_col()];
        for (int i = 0; i < dataset_[view_index].n_col(); i++)
            ris[i] = col_assignment_[view_index][i];
        return ris;
    }

    /**
     * Return the number of column clusters for the view 
     * 
     * @param view_index, the considered view
     * @return the number of column clusters for the view 
     */
    int get_n_col_clusters(int view_index) {
        return n_col_clusters_[view_index];
    }

    /**
     * Return the number of row clusters
     * 
     * @return the number of row clusters
     */
    int get_n_row_clusters() {
        return n_row_clusters_;
    }
    
    /**
     * Generate a random integer between 0 and max
     * 
     * @param max, the max int to generate
     * @return the random number
     */
    int randint(int max) {
        return random_generator_->IRandom(0, max - 1);
    }
    
    /**
     * Convert a numeric value to a string
     * 
     * @param num, the number
     * @return the string representation of num
     */
    std::string to_string(double num) {
        std::ostringstream oss;
        oss << num;
        return oss.str();
    }

    /**
     * Destructor of the class instance. 
     */
    ~MultiViewCoClust() {
        delete random_generator_;
        delete row_assignment_;
        for (int k = 0; k < num_views_; ++k) {
            delete col_assignment_[k];
        }
        delete col_assignment_;
        delete n_col_clusters_;
        delete n_tot_;
        delete n_tot_square_;
        delete two_divided_by_n_tot_;
        delete Ic_;
        delete Cc_;
        for (int k = 0; k < num_views_; ++k) {
            for (int i = 0; i < n_row_clusters_; i++) delete T_[k][i];
            delete T_[k];
        }
        delete T_;
    }
};

#endif
