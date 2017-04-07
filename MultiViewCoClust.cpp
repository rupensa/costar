/*
   The programm is the implementation of the Star-Structured Co-clustering
   approach.

   USAGE:

        ./multiViewCoClu -r num_row -n num_iteration -i configFile -o outputFile
   [-p] [-k min_number_of_cluster]

        configFile: A file that contains the information about the different
   views.
                                A line for each view
                                Each line is in the following format:
                                num_column matrixFileName

        outputFile: An output file without extension.
                                The algorithm produce the following files:
                                - 1 outputFile.row
                                - l outputFile_l.col (one for each view
   indicated in the configFile)


        optional argument

        -p : start with a random initialization with an initial number of
   row/column cluster lesser than the original dimension of the data

        -k : minimum number of cluster, this is a soft constraint. We bound the
   minimum number of row cluster, but not the maximum number of cluster



        Example:
                        ./multiViewCoClu -r 5 -n 100 -i config.txt -k 2 -o
   output

        config.txt:
                6 data.txt
                6 data1.txt

        results:
                output.row
                output_0.col
                output_1.col

        The results are presented in the following format:
                Each file contains a unique row.
                The row has as many elements as the number of elements in that
   dimension.
                Each number indicates the cluster assignment for the element.


        Example:
                output.row:
                                2 0 1 1 2
                5 elements
                element 0 and element 5 belong to cluster number 2
                element 1 belongs to cluster number 0
                element 2 and element 3 belong to cluster 1

        ---

   D. Ienco, C. Robardet, R.G. Pensa, R. Meo. Parameter-Less Co-Clustering for
   Star-Structured Heterogeneous Data.
   Data Min. Knowl. Discov. Vol. 26(2) 2013. pp 217-254. Springer.
 */

#include "MultiViewCoClust.h"

void MultiViewCoClust::updateIandC() {
    // Compute statistic to compute I and C
    // double totalN =0;
    double temp_ir = 0;
    double temp_cr = 0;

    // for each view 
    for (int view = 0; view < num_views_; ++view) {
        // init support lists to 0 
        double *tot_t_per_row_cluster = new double[n_row_clusters_];
        double *tot_t_square_per_row_cluster = new double[n_row_clusters_];
        for (int i = 0; i < n_row_clusters_; i++) tot_t_per_row_cluster[i]
                = tot_t_square_per_row_cluster[i] = 0;

        double *tot_t_per_col_cluster = new double[n_col_clusters_[view]];
        double *tot_t_square_per_col_cluster = new double[n_col_clusters_[view]];
        for (int i = 0; i < n_col_clusters_[view]; i++) tot_t_per_col_cluster[i]
                = tot_t_square_per_col_cluster[i] = 0;

        for (int rc = 0; rc < n_row_clusters_; rc++) {
            for (int cc = 0; cc < n_col_clusters_[view]; cc++) {
                // rowSum[i] is the sum of all data[document][.] 
                // where document is included in the i-th row cluster (considers all features)
                tot_t_per_row_cluster[rc] += T_[view][rc][cc];
                // colSum[j] is the sum of all data[.][feature] 
                // where feature is included in the j-th column cluster (considers all documents)
                tot_t_per_col_cluster[cc] += T_[view][rc][cc];

                tot_t_square_per_col_cluster[cc] += pow(T_[view][rc][cc], 2);
                tot_t_square_per_row_cluster[rc] += pow(T_[view][rc][cc], 2);
            }
        }

        // compute I values for columns
        double sum_i = 0;
        for (int cc = 0; cc < n_col_clusters_[view]; cc++) {
            sum_i += pow(tot_t_per_col_cluster[cc], 2);
        }
        Ic_[view] = 1 - (sum_i / pow(n_tot_[view], 2));

        // compute C values for columns
        double sum_c = 0;
        for (int rc = 0; rc < n_row_clusters_; rc++) {
            if (tot_t_per_row_cluster[rc] != 0) {
                sum_c += (tot_t_square_per_row_cluster[rc] / tot_t_per_row_cluster[rc]);
            }
        }
        Cc_[view] = 1 - (sum_c / n_tot_[view]);

        // compute I value for rows 
        double sum_ir = 0;
        for (int rc = 0; rc < n_row_clusters_; rc++) {
            sum_ir += pow(tot_t_per_row_cluster[rc], 2);
        }
        temp_ir += (sum_ir / pow(n_tot_[view], 2));

        // compute C value for rows
        double sum_cr = 0;
        for (int cc = 0; cc < n_col_clusters_[view]; cc++) {
            if (tot_t_per_col_cluster[cc] != 0) sum_cr += (tot_t_square_per_col_cluster[cc] / tot_t_per_col_cluster[cc]);
        }
        temp_cr += (sum_cr / n_tot_[view]);

        delete tot_t_per_row_cluster;
        delete tot_t_square_per_row_cluster;
        delete tot_t_per_col_cluster;
        delete tot_t_square_per_col_cluster;
    }

    Ir_ = num_views_ - temp_ir;
    Cr_ = num_views_ - temp_cr;
}

/**
 * Compute and print the tau_r and the tau_c values. 
 * Tau_r is tau(X|Y) where Y is the set of column partitions for each view 
 * Tau_c is tau(Yv|X) for each v in views
 * 
 * @return a string containing tau_r and tau_c values separated by spaces. 
 */
std::string MultiViewCoClust::getTauRTauCString() {
    std::string ris = "";
    double **tot_t_per_row_cluster = new double *[num_views_];
    double **tot_t_square_per_row_cluster = new double *[num_views_];

    double **tot_t_per_col_cluster = new double *[num_views_];
    double **tot_t_square_per_col_cluster = new double *[num_views_];

    for (int k = 0; k < num_views_; ++k) {
        tot_t_per_row_cluster[k] = new double[n_row_clusters_];
        tot_t_square_per_row_cluster[k] = new double[n_row_clusters_];
        for (int i = 0; i < n_row_clusters_; i++) tot_t_per_row_cluster[k][i] = tot_t_square_per_row_cluster[k][i] = 0;

        tot_t_per_col_cluster[k] = new double[n_col_clusters_[k]];
        tot_t_square_per_col_cluster[k] = new double[n_col_clusters_[k]];
        for (int i = 0; i < n_col_clusters_[k]; i++) tot_t_per_col_cluster[k][i] = tot_t_square_per_col_cluster[k][i] = 0;
    }

    for (int k = 0; k < num_views_; ++k) {
        for (int i = 0; i < n_row_clusters_; i++) {
            for (int j = 0; j < n_col_clusters_[k]; j++) {
                // To compute Co
                tot_t_square_per_col_cluster[k][j] += pow(T_[k][i][j], 2);
                // To compute C
                tot_t_per_row_cluster[k][i] += T_[k][i][j];
                tot_t_square_per_row_cluster[k][i] += pow(T_[k][i][j], 2);
                // To compute I
                tot_t_per_col_cluster[k][j] += T_[k][i][j];
            }
        }
    }

    double a = 0;
    double b = 0;
    for (int k = 0; k < num_views_; ++k) {
        for (int i = 0; i < n_row_clusters_; i++) {
            a += pow(tot_t_per_row_cluster[k][i], 2) / pow(n_tot_[k], 2);
        }

        for (int i = 0; i < n_row_clusters_; i++) {
            for (int j = 0; j < n_col_clusters_[k]; j++) {
                if (tot_t_per_col_cluster[k][j] != 0)
                    b += (pow(T_[k][i][j], 2) / (tot_t_per_col_cluster[k][j] * n_tot_[k]));
            }
        }
    }

    double tau_r = (b - a) / (num_views_ - a);
    
    ris = to_string(tau_r); 

    for (int view = 0; view < num_views_; ++view) {
        double a1 = 0;
        for (int cc = 0; cc < n_col_clusters_[view]; cc++) {
            a1 += pow(tot_t_per_col_cluster[view][cc], 2) / pow(n_tot_[view], 2);
        }

        double b1 = 0;
        for (int rc = 0; rc < n_row_clusters_; rc++) {
            for (int cc = 0; cc < n_col_clusters_[view]; cc++) {
                double den = (tot_t_per_row_cluster[view][rc] * n_tot_[view]);
                if (den != 0) {
                    b1 += (pow(T_[view][rc][cc], 2) / den);
                }
            }
        }

        double tau_c = (b1 - a1) / (1 - a1);
        ris = ris + " " + to_string(tau_c);
    }

    for (int k = 0; k < num_views_; ++k) {
        delete tot_t_per_row_cluster[k];
        delete tot_t_square_per_row_cluster[k];
        delete tot_t_per_col_cluster[k];
        delete tot_t_square_per_col_cluster[k];
    }
    delete[] tot_t_per_row_cluster;
    delete[] tot_t_square_per_row_cluster;
    delete[] tot_t_per_col_cluster;
    delete[] tot_t_square_per_col_cluster;
    return ris;
}

/**
 * Print the selected view.
 * 
 * @param view_index int, the index of the view to print
 */
void MultiViewCoClust::printView(int view_index) {
    for (int i = 0; i < dataset_[view_index].n_row(); i++) {
        for (int j = 0; j < dataset_[view_index].n_col(); j++) {
            std::cout << dataset_[view_index].getElem(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

/**
 * CoStar discrete initialization method. 
 * Each row (resp. column) belongs to a new row (resp. column) cluster.
 */
void MultiViewCoClust::discreteInitialization() {
    int index = 0;
    std::vector<std::string>::iterator it_str;

    // the number of documents
    int n_row = dataset_[0].n_row();
    // row_str is a vector of strings with n_rows elements. The string at index i
    // is the concatenation
    // of all features values in all views for the element i
    std::vector<std::string> row_str(n_row);
    // hashRow is an hashmap that maps a row string to a cluster id
    std::map<std::string, int> row_string_to_cluster;

    for (int view = 0; view < num_views_; ++view) {
        // clone the matrix related to view k
        double **temp = dataset_[view].cloneMatrix();

        std::map<std::string, int> column_string_to_cluster;

        // is a vector of n_feature elements for view k
        std::vector<std::string> col_str(dataset_[view].n_col());

        // foreach row
        for (int di = 0; di < n_row; di++) {
            std::string line = "";

            // foreach column in the current view
            for (int fi = 0; fi < dataset_[view].n_col(); fi++) {
                // oss always contains one single element dataset[k][i][j], is like a
                // cast to string
                std::string value_str = to_string(temp[di][fi]);

                // add the element to the string representing the j-th column of the
                // matrix
                col_str[fi] = col_str[fi] + value_str + " ";

                // add the element to the string of representing the i-th line of the
                // matrix
                line = line + value_str + " ";
            }

            // concatenate the row string for this view to the global document row
            // string
            std::string temp_str = row_str[di];
            row_str[di] = temp_str + line;

            delete temp[di];
        }
        delete temp;

        /** MANAGE INITIAL COMPRESSION OVER THE COLUMN **/
        //		for (index = 0; index < matrices[k].nCol();++index) {
        //      colAssignment[k][index] = index;
        //    }
        //		nCoC[k] = matrices[k].nCol();

        /**
         * INIT COLUMNS: each column is assigned to a cluster,
         * except for identical columns that are assigned to the same cluster.
         */

        // iterate on all elements in col_str
        for (index = 0, it_str = col_str.begin(); it_str != col_str.end();
                ++index, ++it_str) {
            // if the element col_string[i] is not in the col hashmap, add the pair
            // (col_string[i], i)
            if (column_string_to_cluster.find(*it_str) == column_string_to_cluster.end()) {
                column_string_to_cluster[*it_str] = index;
            }

            // the cluster assignment for the feature i is the column index or the
            // index of an identical column, if present
            col_assignment_[view][index] = column_string_to_cluster[*it_str];
        }

        // map_col is an hashmap with as key the cluster 'number', and values the
        // normalized index of the cluster
        // (in order to get sequential cluster indices)
        std::map<int, int> map_col;
        // for each column
        for (int i = 0; i < dataset_[view].n_col(); ++i) {
            if (map_col.find(col_assignment_[view][i]) == map_col.end()) {
                int actual_index = map_col.size();
                map_col[col_assignment_[view][i]] = actual_index;
            }
            col_assignment_[view][i] = map_col[col_assignment_[view][i]];
        }

        // save the number of clusters created
        n_col_clusters_[view] = map_col.size();
    }

    /**
     * INIT ROWS: each row is assigned to a row cluster,
     * except for identical rows that are assigned to the same row cluster.
     */
    for (index = 0, it_str = row_str.begin(); it_str != row_str.end();
            ++it_str, ++index) {
        if (row_string_to_cluster.find(*it_str) == row_string_to_cluster.end()) {
            row_string_to_cluster[*it_str] = index;
        }
        row_assignment_[index] = row_string_to_cluster[*it_str];
    }

    std::map<int, int> map_row;
    for (int i = 0; i < n_row; ++i) {
        if (map_row.find(row_assignment_[i]) == map_row.end()) {
            int actual_index = map_row.size();
            map_row[row_assignment_[i]] = actual_index;
        }
        row_assignment_[i] = map_row[row_assignment_[i]];
    }
    n_row_clusters_ = map_row.size();

    /**
     * Compute T contingency matrices 
     */
    for (int view = 0; view < num_views_; ++view) {
        for (int di = 0; di < dataset_[view].n_row(); di++) delete[] T_[view][di];
        delete[] T_[view];

        // clone the data matrix 
        double **temp = dataset_[view].cloneMatrix();

        // init a zero T[k] matrix 
        T_[view] = new double *[n_row_clusters_];
        for (int rc = 0; rc < n_row_clusters_; ++rc) {
            T_[view][rc] = new double[n_col_clusters_[view]];
            for (int cc = 0; cc < n_col_clusters_[view]; ++cc) {
                T_[view][rc][cc] = 0;
            }
        }

        // compute the matrix
        for (int di = 0; di < dataset_[view].n_row(); ++di) {
            int cl_row = row_assignment_[di];
            for (int fi = 0; fi < dataset_[view].n_col(); ++fi) {
                int cl_col_id = col_assignment_[view][fi];
                T_[view][cl_row][cl_col_id] += temp[di][fi];
            }
        }
    }
}

/**
 * Initialize the T matrix for the specified view (T[v]).
 * The T[v] matrix has shape = (n_row_clusters, n_col_clusters[v]).
 *
 * @param view : int, the index of the view to be considered
 * @param old_n_row_clusters : int, the old number of row clusters
 * @param new_n_row_clusters : int, the new number of row clusters
 *
 */
void MultiViewCoClust::computeT(int view, int old_n_row_clusters,
                                int new_n_row_clusters,
                                int new_n_col_clusters) {
    // delete existent rows
    for (int i = 0; i < old_n_row_clusters; ++i) delete T_[view][i];
    delete T_[view];

    // init new rows list
    T_[view] = new double *[new_n_row_clusters];
    // init each cell with 0
    for (int i = 0; i < new_n_row_clusters; ++i) {
        T_[view][i] = new double[new_n_col_clusters];
        for (int j = 0; j < new_n_col_clusters; ++j) T_[view][i][j] = 0;
    }

    for (int i = 0; i < dataset_[view].n_row(); ++i) {
        double *temp = dataset_[view].getRowGroupByColClust(i,
                                                            col_assignment_[view],
                                                            new_n_col_clusters);
        for (int j = 0; j < new_n_col_clusters; ++j)
            T_[view][row_assignment_[i]][j] += temp[j];
        delete temp;
    }
}

/**
 * Randomly initialize row clusters.
 * 
 * @return the number of created row clusters
 */
int MultiViewCoClust::randomInitRow() {
    // get the argument or 2 if is too low
    int min = (min_n_row_cluster_ > 2) ? min_n_row_cluster_ : 2;

    // 2 <= nCoR_new <= (number of documents in dataset / 2)
    int new_n_row_clusters = (dataset_[0].n_row() > 3)
            ? random_generator_->IRandomX(min, dataset_[0].n_row() / 2)
            : random_generator_->IRandomX(min, dataset_[0].n_row());

    // in this way I guarantee that for each cluster there is at least one element
    // for both row and column cluster
    for (int i = 0; i < new_n_row_clusters; ++i) row_assignment_[i] = i;
    for (int i = new_n_row_clusters; i < dataset_[0].n_row(); ++i)
        row_assignment_[i] = randint(new_n_row_clusters);

    return new_n_row_clusters;
}

/**
 * Randomly initialize the cluster assignment for columns of a specific view.
 *
 * @param viewIndex : int, the index of the view to consider
 * @return int, the number of column clusters created
 */
int MultiViewCoClust::randomInitCol(int view_index) {
    int new_n_col_clusters = (dataset_[view_index].n_col() > 3)
            ? random_generator_->IRandomX(2, dataset_[view_index].n_col() / 2)
            : random_generator_->IRandomX(2, dataset_[view_index].n_col());
    for (int i = 0; i < new_n_col_clusters; ++i) col_assignment_[view_index][i] = i;
    for (int i = new_n_col_clusters; i < dataset_[view_index].n_col(); ++i)
        col_assignment_[view_index][i] = randint(new_n_col_clusters);
    return new_n_col_clusters;
}

/**
 * Randomly initialize rows and columns clusters.
 */
void MultiViewCoClust::randomInitialization() {
    int new_n_row_clusters = randomInitRow();

    n_row_clusters_ = new_n_row_clusters;

    for (int k = 0; k < num_views_; ++k) {
        // update column clusters for view k
        int new_n_col_clusters = randomInitCol(k);
        n_col_clusters_[k] = new_n_col_clusters;
        // update the contingency matrix
        computeT(k, n_row_clusters_, new_n_row_clusters, new_n_col_clusters);

        // update matrices related to T (for view k)
        n_tot_[k] = 0;
        for (int i = 0; i < n_row_clusters_; ++i) {
            for (int j = 0; j < n_col_clusters_[k]; ++j) n_tot_[k] += T_[k][i][j];
        }
        n_tot_square_[k] = pow(n_tot_[k], 2);
    }
}

/**
 * Print the cluster assignments. 
 */
void MultiViewCoClust::printAssignment() {
    std::cout << "row assignment " << n_row_clusters_ << std::endl;
    for (int i = 0; i < dataset_[0].n_row(); ++i)
        std::cout << row_assignment_[i] << " ";
    std::cout << std::endl;

    for (int k = 0; k < num_views_; ++k) {
        std::cout << "col assignment " << n_col_clusters_[k] << std::endl;
        for (int i = 0; i < dataset_[k].n_col(); ++i)
            std::cout << col_assignment_[k][i] << " ";
        std::cout << std::endl;
    }
    getchar();
}

/**
 * Print the contingency matrix of clusters for each view. 
 */
void MultiViewCoClust::printContingencyT() {
    for (int view = 0; view < num_views_; ++view) {
        std::cout << "======================================" << std::endl;
        std::cout << "T[" << view << "] matrix" << std::endl;
        for (int rc = 0; rc < get_n_row_clusters(); rc++) {
            for (int cc = 0; cc < get_n_col_clusters(view); cc++) {
                std::cout << T_[view][rc][cc] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "======================================" << std::endl;
    }
    // getchar();
}

/**
 * Entry method to run the CoStar algorithm.
 *
 * randomPart : 0 or 1,
 *              1 to randomly initialize partitions,
 *              0 to use discrete partitions
 *
 * n_iter : int, the number of iterations to perform
 *
 */
void MultiViewCoClust::buildCoClu(int random_initialization, int n_iter) {
    clock_t start, end;
    // std::cout<<"runNam: "<<outputFile<<std::endl;

    start = clock();
    if (random_initialization) {
        randomInitialization();
        // std::cout<<"random partition"<<std::endl;
    }
    else {
        // std::cout<<"discrete partition"<<std::endl;
        discreteInitialization();
    }

    // std::cout<<"parto con "<<std::endl;
    // std::cout<<"ROW CLUSTER NUMBER: "<<getCoR()<<std::endl;
    // for (int k = 0; k < num_view; ++k)
    // std::cout<<"COLUMN CLUSTER "<<k<<" NUMBER: "<<getCoC(k)<<std::endl;

    srand(time(NULL));

    int counter = 0;

    updateIandC();

    std::ofstream iterations_file((output_filename_ + ".iter").c_str());
    iterations_file << counter << " " << getTauRTauCString() << std::endl;
    std::cout << counter << " " << getTauRTauCString() << std::endl;

    while (counter != n_iter) {
        rowPartitionCluster();
        for (int k = 0; k < num_views_; ++k) {
            columnPartitionCluster(k);
        }
        counter++;
        iterations_file << counter << " " << getTauRTauCString() << std::endl;
        // std::cout<< "Dino: " << counter<<std::endl;
        // std::cout<< counter << std::endl;

        // if (counter % 1000 == 0){
        //	std::cout<< counter << std::endl;
        // for (int k = 0; k < num_view; ++k) std::cout<<nCoC[k]<<" ";
        // end = clock();
        // std::cout<<" --- time "<<( end - start )  <<std::endl;
        // std::cout.flush();
        // start = clock();
        //}
    }
    end = clock();
    iterations_file.close();
    std::string time_name = output_filename_ + ".time";
    std::ofstream timeTau(time_name.c_str());
    double seconds = ((double) (end - start)) / CLOCKS_PER_SEC;
    timeTau << seconds << std::endl;
    timeTau.close();
    timeTau.close();

    std::ofstream row((output_filename_ + ".row").c_str());

    for (int i = 0; i < dataset_[0].n_row(); i++) row << row_assignment_[i] << " ";
    row << std::endl;
    row.close();

    for (int k = 0; k < num_views_; ++k) {
        std::stringstream out;
        out << output_filename_ << "_" << k << ".col";
        std::ofstream col(out.str().c_str());
        for (int i = 0; i < dataset_[k].n_col(); i++)
            col << col_assignment_[k][i] << " ";
        col << std::endl;
        col.close();
    }
}

/**
 * Given a cluster index and a cluster assignments, extract a random element
 * that belongs to the specified cluster
 * 
 * @param cluster : int, the cluster index
 * @param assignment : array, the clusters assignments
 * @param assignment_size : int, the size of the assignment array
 * @return int, the index of the extracted element
 */
int MultiViewCoClust::chooseRandomElementFromCluster(int cluster, int *assignment,
                                                     int assignment_size) {
    double *elements_of_cluster = new double[assignment_size];
    int n_elem = 0;
    for (int i = 0; i < assignment_size; i++) {
        if (assignment[i] == cluster) elements_of_cluster[n_elem++] = i;
    }
    int ris = elements_of_cluster[randint(n_elem)];
    delete elements_of_cluster;
    return ris;
}

/**
 * Partition on rows.
 */
void MultiViewCoClust::rowPartitionCluster() {
    // Work on the row partition
    // partition on the ROW
    double min_tau_r = 0;
    int cluster_min = -1;

    // choose random a row cluster
    int source_cluster = random_generator_->IRandomX(1, n_row_clusters_ * 100) % n_row_clusters_;
    // choose at random an element of row cluster random_cluster
    int random_element = chooseRandomElementFromCluster(source_cluster, row_assignment_, dataset_[0].n_row());

    double **lambda = new double *[num_views_];
    double **tot_t_per_col_cluster = new double *[num_views_];
    double *tot_lambda = new double[num_views_];

    // computes support vectors and matrices
    for (int view = 0; view < num_views_; ++view) {

        lambda[view] = dataset_[view].getRowGroupByColClust(random_element, col_assignment_[view], n_col_clusters_[view]);

        tot_t_per_col_cluster[view] = new double[n_col_clusters_[view]];
        for (int cc = 0; cc < n_col_clusters_[view]; ++cc) tot_t_per_col_cluster[view][cc] = 0;
        tot_lambda[view] = 0;
        for (int cc = 0; cc < n_col_clusters_[view]; ++cc) {
            tot_lambda[view] += lambda[view][cc];
            for (int rc = 0; rc < n_row_clusters_; ++rc) {
                tot_t_per_col_cluster[view][cc] += T_[view][rc][cc];
            }
        }
    }

    int ncluster_with_empty =
            (n_row_clusters_ + 1 > dataset_[0].n_row()) ? dataset_[0].n_row() : n_row_clusters_ + 1;

    std::vector<int> equal_solutions;

    for (int evaluated_cluster = 0; evaluated_cluster < ncluster_with_empty; evaluated_cluster++) {
        // if the cluster is not the same cluster from which I extract the element x
        // build a temporary assignment to evaluate this swap
        if (evaluated_cluster != source_cluster) {
            double delta_tau_r = deltaTauR(tot_t_per_col_cluster, tot_lambda, lambda, source_cluster, evaluated_cluster);
            if (delta_tau_r == min_tau_r) {
                equal_solutions.push_back(evaluated_cluster);
            }
            if (delta_tau_r < min_tau_r) {
                min_tau_r = delta_tau_r;
                cluster_min = evaluated_cluster;
                equal_solutions.clear();
                equal_solutions.push_back(cluster_min);
            }
        }
    }

    if (equal_solutions.size() > 1 && min_tau_r != 0) {
        // choose a row from the list of equal results 
        cluster_min = findNonDominatedRow(equal_solutions, lambda, source_cluster);
    }

    int go_on_normally = 1;
    if (n_row_clusters_ == min_n_row_cluster_) {
        // if the number of cluster is already equal to the min number of requested row clusters
        // check that the operation will not remove a row cluster otherwise skip this operation
        go_on_normally = checkRowClustConstraint(source_cluster);
    }

    if (min_tau_r != 0 && go_on_normally == 1) {
        row_assignment_[random_element] = cluster_min;
        modifyRowCluster(lambda, source_cluster, cluster_min);
        updateIandC();
    }
    for (int k = 0; k < num_views_; ++k) {
        delete lambda[k];
        delete tot_t_per_col_cluster[k];
    }
    delete tot_lambda;
    delete tot_t_per_col_cluster;
    delete lambda;
}

/**
 * Computes the tau_r value. 
 * 
 * @param t, the contingency matrix
 * @param temp_n_col_clusters int, the number of column clusters in the specified view
 * @param view_index int, the view index 
 * @return double, the computed value
 */
double MultiViewCoClust::getTauR(double **t, int temp_n_col_clusters, int view_index) {
    double *tot_t_per_row_cluster = new double[n_row_clusters_];
    // double *tot_t_square_per_row_cluster = new double[n_row_clusters_];

    double *tot_t_per_col_cluster = new double[temp_n_col_clusters];
    // double *tot_t_square_per_col_cluster = new double[temp_n_col_clusters];

    for (int i = 0; i < n_row_clusters_; i++) { 
        tot_t_per_row_cluster[i] = 0; // tot_t_square_per_row_cluster[i] = 0; 
    }

    for (int i = 0; i < temp_n_col_clusters; i++) {
        tot_t_per_col_cluster[i] = 0; // tot_t_square_per_col_cluster[i] = 0;
    }

    for (int i = 0; i < n_row_clusters_; i++) {
        for (int j = 0; j < temp_n_col_clusters; j++) {
            // To compute Co
            // tot_t_square_per_col_cluster[j] += pow(t[i][j], 2);
            // To compute C
            tot_t_per_row_cluster[i] += t[i][j];
            // tot_t_square_per_row_cluster[i] += pow(t[i][j], 2);
            // To compute I
            tot_t_per_col_cluster[j] += t[i][j];
        }
    }

    double a = 0;
    for (int i = 0; i < n_row_clusters_; i++) {
        a += pow(tot_t_per_row_cluster[i], 2) / pow(n_tot_[view_index], 2);
    }

    double b = 0;
    for (int i = 0; i < n_row_clusters_; i++) {
        for (int j = 0; j < temp_n_col_clusters; j++) {
            b += (pow(t[i][j], 2) / (tot_t_per_col_cluster[j] * n_tot_[view_index]));
        }
    }

    double tau_r = (b - a) / (1 - a);

    delete[] tot_t_per_row_cluster;
    // delete[] tot_t_square_per_row_cluster;
    delete[] tot_t_per_col_cluster;
    // delete[] tot_t_square_per_col_cluster;

    return tau_r;
}

/**
 * Computes the tau_c value. 
 * 
 * @param t, the contingency matrix
 * @param temp_n_row_clusters int, the number of row clusters
 * @return double, the computed value
 */
double *MultiViewCoClust::getTauC(double ***t, int temp_n_row_clusters) {
    double **tot_t_per_row_cluster = new double *[num_views_];
    double **tot_t_per_col_cluster = new double *[num_views_];

    for (int view = 0; view < num_views_; ++view) {
        tot_t_per_row_cluster[view] = new double[temp_n_row_clusters];
        for (int i = 0; i < temp_n_row_clusters; i++) {
            tot_t_per_row_cluster[view][i] = 0;
        }

        tot_t_per_col_cluster[view] = new double[n_col_clusters_[view]];
        for (int cc = 0; cc < n_col_clusters_[view]; cc++) {
            tot_t_per_col_cluster[view][cc] = 0;
        }
    }

    for (int view = 0; view < num_views_; ++view) {
        for (int rc = 0; rc < temp_n_row_clusters; rc++) {
            for (int cc = 0; cc < n_col_clusters_[view]; cc++) {
                // To compute C
                tot_t_per_row_cluster[view][rc] += t[view][rc][cc];
                // To compute I
                tot_t_per_col_cluster[view][cc] += t[view][rc][cc];
            }
        }
    }

    double *ris = new double[num_views_];

    for (int view = 0; view < num_views_; ++view) {
        double a1 = 0;
        for (int cc = 0; cc < n_col_clusters_[view]; cc++) {
            a1 += pow(tot_t_per_col_cluster[view][cc], 2) / pow(n_tot_[view], 2);
        }

        double b1 = 0;
        for (int rc = 0; rc < temp_n_row_clusters; rc++) {
            for (int cc = 0; cc < n_col_clusters_[view]; cc++) {
                double den = (tot_t_per_row_cluster[view][rc] * n_tot_[view]);
                if (den != 0) {
                    b1 += (pow(t[view][rc][cc], 2) / den);
                }
            }
        }
        ris[view] = (b1 - a1) / (1 - a1);
    }

    for (int k = 0; k < num_views_; ++k) {
        delete tot_t_per_row_cluster[k];
        delete tot_t_per_col_cluster[k];
    }
    delete[] tot_t_per_row_cluster;
    delete[] tot_t_per_col_cluster;

    return ris;
}

/**
 * Compute the emulated tau_r value, if a feature, described by the lambda 
 * vector, is moved from the source to the destination column cluster. 
 * 
 * @param lambda array of double, the list of changes for each row cluster 
 * if moving the specified element from the source to the destination column cluster
 * @param source_cluster int, the source column cluster
 * @param destination_cluster int, the destination column cluster
 * @param view_index int, the considered view
 * @return the computed tau_r value
 */
double MultiViewCoClust::computeEmulatedTauR(double *lambda, int source_cluster,
                                             int destination_cluster,
                                             int view_index) {
    int temp_n_col_clusters = n_col_clusters_[view_index];
    double **temp_t = new double *[n_row_clusters_];

    if (destination_cluster == n_col_clusters_[view_index]) {
        for (int rc = 0; rc < n_row_clusters_; rc++) {
            temp_t[rc] = new double[temp_n_col_clusters + 1];
            for (int cc = 0; cc < temp_n_col_clusters; cc++) {
                temp_t[rc][cc] = T_[view_index][rc][cc];
                // subtract from the original col cluster b of the element y its value
                if (cc == source_cluster) {
                    temp_t[rc][cc] -= lambda[rc];
                }
            }
            temp_t[rc][destination_cluster] = lambda[rc];
        }
        temp_n_col_clusters++;
    }
    else {
        // we move the object x from cluster b to cluster e
        for (int rc = 0; rc < n_row_clusters_; rc++) {
            temp_t[rc] = new double[temp_n_col_clusters];
            for (int cc = 0; cc < temp_n_col_clusters; cc++) {
                temp_t[rc][cc] = T_[view_index][rc][cc];
            }
        }

        for (int rc = 0; rc < n_row_clusters_; rc++) {
            temp_t[rc][destination_cluster] += lambda[rc];
            temp_t[rc][source_cluster] -= lambda[rc];
        }
        ///////checkEmptyColCluster////
        int *counter = new int[temp_n_col_clusters];
        for (int cc = 0; cc < temp_n_col_clusters; ++cc) counter[cc] = 0;
        for (int col_index = 0; col_index < dataset_[view_index].n_col(); ++col_index)
            counter[col_assignment_[view_index][col_index]]++;
        counter[source_cluster] = counter[source_cluster] - 1;
        counter[destination_cluster] = counter[destination_cluster] + 1;
        int is_cluster_empty = (counter[source_cluster] == 0) ? 1 : 0;
        delete counter;
        ///////////
        // if there is one empty cluster we compact the matrix, adjust the
        // assignment
        if (is_cluster_empty) {
            // if the row cluster is the empty one we don't assign the values to the
            // new tempT, for this reason
            // we check another index (k) and increment this index only when there is
            // an assignment to tempT
            for (int i = 0; i < n_row_clusters_; i++) {
                delete temp_t[i];
                temp_t[i] = new double[temp_n_col_clusters - 1];
                for (int j = 0, k = 0; j < temp_n_col_clusters; j++) {
                    if (j != source_cluster) {
                        temp_t[i][k] = T_[view_index][i][j];
                        k++;
                    }
                }
            }
            temp_n_col_clusters--;
        }
    }

    double ris = getTauR(temp_t, temp_n_col_clusters, view_index);

    for (int i = 0; i < n_row_clusters_; i++) delete temp_t[i];
    delete temp_t;
    return ris;
}

/**
 * Compute the emulated tau_c value, if an element, described by the lambda 
 * vector, is moved from the source to the destination row cluster. 
 * 
 * @param lambda array of double, the list of changes for each column cluster 
 * if moving the specified element from the source to the destination row cluster
 * @param source_cluster int, the source row cluster
 * @param destination_cluster int, the destination row cluster
 * @return the computed tau_c value
 */
double *MultiViewCoClust::computeEmulatedTauC(double **lambda,
                                              int source_cluster,
                                              int destination_cluster) {
    int temp_n_row_clusters = n_row_clusters_;
    double ***temp_t = new double **[num_views_];
    // we check if we must add a new row cluster
    if (destination_cluster == temp_n_row_clusters) {
        // std::cout<<"\t aggiungo un nuovo row cluster"<<std::endl;
        for (int view = 0; view < num_views_; ++view) {
            double **temp = new double *[temp_n_row_clusters + 1];

            for (int rc = 0; rc < temp_n_row_clusters; rc++) {
                temp[rc] = new double[n_col_clusters_[view]];
                for (int cc = 0; cc < n_col_clusters_[view]; cc++) {
                    temp[rc][cc] = T_[view][rc][cc];
                    // subtract from the original row cluster b of the element x its value
                    if (rc == source_cluster) {
                        temp[rc][cc] -= lambda[view][cc];
                    }
                }
            }
            temp[destination_cluster] = new double[n_col_clusters_[view]];
            for (int cc = 0; cc < n_col_clusters_[view]; cc++) {
                temp[destination_cluster][cc] = lambda[view][cc];
            }
            temp_t[view] = temp;
        }
        temp_n_row_clusters++;
    }
    else {
        // we move the object x from cluster b to cluster e
        for (int view = 0; view < num_views_; ++view) {
            double **temp = new double *[temp_n_row_clusters];
            for (int rc = 0; rc < temp_n_row_clusters; rc++) {
                temp[rc] = new double[n_col_clusters_[view]];
                for (int cc = 0; cc < n_col_clusters_[view]; cc++) {
                    temp[rc][cc] = T_[view][rc][cc];
                }
            }
            temp_t[view] = temp;
        }

        for (int view = 0; view < num_views_; ++view) {
            for (int cc = 0; cc < n_col_clusters_[view]; cc++) {
                temp_t[view][destination_cluster][cc] += lambda[view][cc];
                temp_t[view][source_cluster][cc] -= lambda[view][cc];
            }
        }

        ///////checkEmptyRowCluster////
        int *t_ = new int[temp_n_row_clusters];
        for (int i = 0; i < n_row_clusters_; ++i) t_[i] = 0;
        for (int i = 0; i < dataset_[0].n_row(); ++i) t_[row_assignment_[i]]++;
        t_[source_cluster] = t_[source_cluster] - 1;
        t_[destination_cluster] = t_[destination_cluster] + 1;
        int is_cluster_empty = (t_[source_cluster] == 0) ? 1 : 0;
        delete t_;
        ///////////

        // if there is one empty cluster we compact the matrix, adjust the
        // assignment
        if (is_cluster_empty) {
            for (int view = 0; view < num_views_; ++view) {
                double **temp = new double *[n_row_clusters_ - 1];
                // if the row cluster is the empty one we don't assign the values to the
                // new tempT, for this reason
                // we check another index (k) and increment this index only when there
                // is an assignment to tempT
                for (int rc = 0, new_rc = 0; rc < n_row_clusters_; rc++) {
                    if (rc != source_cluster) {
                        temp[new_rc] = new double[n_col_clusters_[view]];
                        for (int cc = 0; cc < n_col_clusters_[view]; cc++) {
                            temp[new_rc][cc] = T_[view][rc][cc];
                        }
                        new_rc++;
                    }
                }
                temp_t[view] = temp;
            }
            temp_n_row_clusters--;
        }
    }

    double *ris = getTauC(temp_t, temp_n_row_clusters);

    // clean allocated memory 
    for (int view = 0; view < num_views_; ++view) {
        for (int rc = 0; rc < temp_n_row_clusters; rc++) delete temp_t[view][rc];
        delete temp_t[view];
    }
    delete temp_t;

    // return the result 
    return ris;
}

/**
 * Select a non dominated column from the proposed equal solutions. 
 * 
 * @param equal_solutions vector<int>, the indices of the possible destination clusters
 * @param lambda, the lambda values for the selected feature
 * @param source_cluster, the current column cluster of the feature
 * @param view_index, the index of the considered view 
 * @return the selected destination column cluster
 */
int MultiViewCoClust::findNonDominatedCol(std::vector<int> &equal_solutions,
                                          double *lambda, int source_cluster,
                                          int view_index) {

    // find the destination cluster that maximizes the tau_r value
    int best_destination = equal_solutions[0];
    double tau_r_best = computeEmulatedTauR(lambda, source_cluster, best_destination, view_index);

    for (int i = 1; i < equal_solutions.size(); ++i) {
        int evaluated_destination = equal_solutions[i];
        double tau_r_evaluated =
                computeEmulatedTauR(lambda, source_cluster, evaluated_destination, view_index);

        if (tau_r_best < tau_r_evaluated) {
            best_destination = evaluated_destination;
            tau_r_best = tau_r_evaluated;
        }
    }
    return best_destination;
}

/**
 * Select a non dominated row from the proposed equal solutions. 
 * 
 * @param equal_solutions vector<int>, the indices of the possible destination clusters
 * @param lambda, the lambda values for the selected element
 * @param source_cluster, the current row cluster of the element
 * @return the selected destination row cluster
 */
int MultiViewCoClust::findNonDominatedRow(std::vector<int> &equal_solutions,
                                          double **lambda, int source_cluster) {

    std::vector<int>::iterator it;
    int best_destination = equal_solutions[0];
    double *tau_c_best = computeEmulatedTauC(lambda, source_cluster, best_destination);

    for (int i = 1; i < equal_solutions.size(); ++i) {
        int evaluated_destination = equal_solutions[i];
        double *tau_c_evaluated = computeEmulatedTauC(lambda, source_cluster, evaluated_destination);

        // count the number of views where the tau features for the current best is better
        // or where the evaluated solution is better, equals cases are not considered
        int best_count = 0;
        int eval_count = 0;
        for (int j = 0; j < num_views_; ++j) {
            if (tau_c_best[j] > tau_c_evaluated[j]) best_count++;
            if (tau_c_best[j] < tau_c_evaluated[j]) eval_count++;
        }

        if (best_count < eval_count) {
            best_destination = evaluated_destination;
            delete tau_c_best;
            tau_c_best = tau_c_evaluated;
        }
        else {
            delete tau_c_evaluated;
        }
    }
    delete tau_c_best;
    return best_destination;
}

/**
 * Check the violation of the non-empty row cluster constraint. 
 *       
 * @param start_cluster, the source cluster
 * @return 0 if the next moving create an empty cluster, 
 *         1 if the next moving doesn't create an empty cluster
 */
int MultiViewCoClust::checkRowClustConstraint(int start_cluster) {
    std::map<int, int> count;
    for (int i = 0; i < dataset_[0].n_row(); ++i) {
        if (count.find(row_assignment_[i]) == count.end()) {
            count[row_assignment_[i]] = 1;
        }
        else {
            count[row_assignment_[i]] += 1;
        }
    }
    return (count[start_cluster] == 1) ? 0 : 1;
}

/** 
 * Check the violation of the non-empty col cluster constraint on a specific view.  
 * 
 * @param view_index, int, the considered view
 * @param start_cluser, int, the source cluster 
 * @return 0 if the next moving create an empty cluster
 *         1 if the next moving doesn't create an empty cluster
 */
int MultiViewCoClust::checkColClustConstraint(int viewIndex, int start_cluster) {
    std::map<int, int> count;
    for (int i = 0; i < dataset_[viewIndex].n_col(); ++i) {
        if (count.find(col_assignment_[viewIndex][i]) == count.end()) {
            count[col_assignment_[viewIndex][i]] = 1;
        }
        else {
            count[col_assignment_[viewIndex][i]] += 1;
        }
    }
    return (count[start_cluster] == 1) ? 0 : 1;
}

/**
 * Partitions on columns. 
 * 
 * @param view_index, int, the considered view
 */
void MultiViewCoClust::columnPartitionCluster(int view_index) {
    double min_tau_c = 0;
    int destination_cluster = -1;
    int num = random_generator_->IRandomX(1, n_col_clusters_[view_index] * 100);
    // choose random a column cluster
    int source_cluster = num % n_col_clusters_[view_index];
    // choose at random an element of column cluster Cb
    int element = chooseRandomElementFromCluster(source_cluster, col_assignment_[view_index],
                                                 dataset_[view_index].n_col());
    double *lambda =
            dataset_[view_index].getColGroupByRowClust(element, row_assignment_, n_row_clusters_);

    int ncluster_with_empty = (n_col_clusters_[view_index] + 1 > dataset_[view_index].n_col())
            ? dataset_[view_index].n_col()
            : n_col_clusters_[view_index] + 1;

    double *tot_t_per_row_cluster = new double[n_row_clusters_];
    double tot_lambda = 0;
    for (int rc = 0; rc < n_row_clusters_; ++rc) {
        tot_lambda += lambda[rc];
        tot_t_per_row_cluster[rc] = 0;
        for (int cc = 0; cc < n_col_clusters_[view_index]; ++cc) {
            tot_t_per_row_cluster[rc] += T_[view_index][rc][cc];
        }
    }

    std::vector<int> equal_solutions;

    for (int evaluated_cluster = 0; evaluated_cluster < ncluster_with_empty; evaluated_cluster++) {
        // if the cluster is not the same cluster from which I extract the element y
        // build a temporary assignment to evaluate this swap
        if (evaluated_cluster != source_cluster) {
            double delta_tau_c = deltaTauC(view_index, tot_t_per_row_cluster, tot_lambda,
                                           lambda, source_cluster, evaluated_cluster);
            if (delta_tau_c == min_tau_c) equal_solutions.push_back(evaluated_cluster);

            if (delta_tau_c < min_tau_c) {
                min_tau_c = delta_tau_c;
                destination_cluster = evaluated_cluster;
            }
        }
    }
    delete tot_t_per_row_cluster;

    if (equal_solutions.size() > 1 && min_tau_c != 0) {
        destination_cluster = findNonDominatedCol(equal_solutions, lambda,
                                                  source_cluster, view_index);
    }

    int go_on_normally = 1;
    if (n_col_clusters_[view_index] == 2) {
        go_on_normally = checkColClustConstraint(view_index, source_cluster);
    }
    if (min_tau_c != 0 && go_on_normally == 1) {
        //		if (minTQ!=0){
        col_assignment_[view_index][element] = destination_cluster;
        modifyColumnCluster(view_index, lambda, source_cluster, destination_cluster);
        updateIandC();
    }
    delete lambda;
}

/**
 * Move a feature from the source to the destination column cluster
 *
 * @param view_index, int, the considered view
 * @param lambda, the lambda values for the feature to move 
 * @param source_cluster, int, the current column cluster
 * @param destination_cluster, the destination column cluster
 */
void MultiViewCoClust::modifyColumnCluster(int view_index, double *lambda,
                                           int source_cluster,
                                           int destination_cluster) {
    if (destination_cluster == n_col_clusters_[view_index]) {
        // the destination cluster is a new column cluster
        double **temp_t = new double *[n_row_clusters_];
        for (int rc = 0; rc < n_row_clusters_; rc++) {
            temp_t[rc] = new double[n_col_clusters_[view_index] + 1];
            for (int cc = 0; cc < n_col_clusters_[view_index]; cc++) {
                temp_t[rc][cc] = T_[view_index][rc][cc];
                // subtract from the original col cluster b of the element y its value
                if (cc == source_cluster) {
                    temp_t[rc][source_cluster] -= lambda[rc];
                }
            }
            temp_t[rc][destination_cluster] = lambda[rc];
            delete T_[view_index][rc];
        }
        delete T_[view_index];
        T_[view_index] = temp_t;
        n_col_clusters_[view_index]++;
    }
    else {
        // we move the object x from cluster b to cluster e
        for (int rc = 0; rc < n_row_clusters_; rc++) {
            T_[view_index][rc][destination_cluster] += lambda[rc];
            T_[view_index][rc][source_cluster] -= lambda[rc];
        }

        int is_cluster_empty = checkEmptyColCluster(view_index, source_cluster);
        // if there is one empty cluster we compact the matrix, adjust the
        // assignment
        if (is_cluster_empty) {
            double **tempT = new double *[n_row_clusters_];
            // if the row cluster is the empty one we don't assign the values to the
            // new tempT, for this reason
            // we check another index (k) and increment this index only when there is
            // an assignment to tempT
            for (int rc = 0; rc < n_row_clusters_; rc++) {
                tempT[rc] = new double[n_col_clusters_[view_index] - 1];
                for (int cc = 0, new_cc = 0; cc < n_col_clusters_[view_index]; cc++) {
                    if (cc != source_cluster) {
                        tempT[rc][new_cc] = T_[view_index][rc][cc];
                        new_cc++;
                    }
                }
                delete T_[view_index][rc];
            }
            delete T_[view_index];

            for (int f = 0; f < dataset_[view_index].n_col(); f++) {
                if (col_assignment_[view_index][f] > source_cluster) col_assignment_[view_index][f]--;
            }

            T_[view_index] = tempT;
            n_col_clusters_[view_index]--;
        }
    }
}

/**
 * Check if a row cluster is empty
 * 
 * @param cluster, int, the index of the row cluster to check
 * @return 1 if the cluster b is empty, 0 otherwise.
 */
int MultiViewCoClust::checkEmptyRowCluster(int cluster) {
    int *count = new int[n_row_clusters_];
    for (int rc = 0; rc < n_row_clusters_; ++rc) count[rc] = 0;
    for (int i = 0; i < dataset_[0].n_row(); ++i) count[row_assignment_[i]]++;
    int ris = (count[cluster] == 0) ? 1 : 0;
    delete count;
    return ris;
}

/**
 * Check if a column cluster is empty
 * 
 * @param view_index, int, the considered view
 * @param cluster, int, the index of the column cluster to check 
 * @return 1 if the cluster b in view viewIndex is empty, 0 otherwise.
 */
int MultiViewCoClust::checkEmptyColCluster(int view_index, int cluster) {
    int *count = new int[n_col_clusters_[view_index]];
    for (int cc = 0; cc < n_col_clusters_[view_index]; ++cc) count[cc] = 0;
    for (int i = 0; i < dataset_[view_index].n_col(); ++i)
        count[col_assignment_[view_index][i]]++;
    int ris = (count[cluster] == 0) ? 1 : 0;
    delete count;
    return ris;
}

/**
 * Moves an element from the source to the destination row cluster
 *
 * @param lambda, the lambda values for the feature to move 
 * @param source_row_cluster, int, the current row cluster
 * @param destination_cluster, the destination row cluster
 */
void MultiViewCoClust::modifyRowCluster(double **lambda, int source_row_cluster,
                                        int destination_cluster) {
    if (destination_cluster == n_row_clusters_) {
        // it is necessary to add a new row cluster
        for (int view = 0; view < num_views_; ++view) {
            double **tempT = new double *[n_row_clusters_ + 1];

            // clone T[k] in tempT
            for (int rc = 0; rc < n_row_clusters_; rc++) {
                tempT[rc] = new double[n_col_clusters_[view]];
                for (int cc = 0; cc < n_col_clusters_[view]; cc++) {
                    tempT[rc][cc] = T_[view][rc][cc];
                }
                delete T_[view][rc];
            }
            delete T_[view];

            tempT[destination_cluster] = new double[n_col_clusters_[view]];
            for (int cc = 0; cc < n_col_clusters_[view]; cc++) {
                // subtract from the original row cluster b the value of the moved element
                tempT[source_row_cluster][cc] -= lambda[view][cc];
                tempT[destination_cluster][cc] = lambda[view][cc];
            }

            // replace T[k] with the new value
            T_[view] = tempT;
        }

        // update the number of row clusters
        n_row_clusters_++;
    }
    else {
        // we move the object x from cluster b to cluster e
        for (int view = 0; view < num_views_; ++view) {
            for (int cc = 0; cc < n_col_clusters_[view]; cc++) {
                T_[view][destination_cluster][cc] += lambda[view][cc];
                T_[view][source_row_cluster][cc] -= lambda[view][cc];
            }
        }

        int is_source_empty = checkEmptyRowCluster(source_row_cluster);
        // if there is one empty cluster we compact the matrix, adjust the
        // assignment
        if (is_source_empty) {
            for (int view = 0; view < num_views_; ++view) {
                double **tempT = new double *[n_row_clusters_ - 1];
                // if the row cluster is the empty one we don't assign the values to the
                // new tempT, for this reason
                // we check another index (k) and increment this index only when there
                // is an assignment to tempT
                for (int rc = 0, new_rc = 0; rc < n_row_clusters_; rc++) {
                    if (rc != source_row_cluster) {
                        tempT[new_rc] = new double[n_col_clusters_[view]];
                        for (int j = 0; j < n_col_clusters_[view]; j++) tempT[new_rc][j] = T_[view][rc][j];
                        new_rc++;
                    }
                    delete T_[view][rc];
                }
                delete T_[view];
                T_[view] = tempT;
            }
            for (int row_index = 0; row_index < dataset_[0].n_row(); row_index++) {
                if (row_assignment_[row_index] > source_row_cluster) {
                    row_assignment_[row_index]--;
                }
            }
            n_row_clusters_--;
        }
    }
}

/**
 * Computes the deltaTauR value. 
 * 
 * @param tot_t_per_col_cluster, the sum of data group by column cluster
 * @param tot_lambda, the sum of lambda
 * @param lambda is the projection of the element x on the column cluster assignment
 * @param source_row_cluster is the row cluster that contains originally x
 * @param destination_row_cluster is the row cluster in which we move x
 * @return double, the computed value
 */
double MultiViewCoClust::deltaTauR(double **tot_t_per_col_cluster, double *tot_lambda,
                                   double **lambda, int source_row_cluster,
                                   int destination_row_cluster) {

    double *X_views = new double[num_views_];
    double *Y_views = new double[num_views_];
    for (int v = 0; v < num_views_; ++v) {
        X_views[v] = Y_views[v] = 0;
    }
    // we check if the new label of the element x is the empty set, this means
    // that we must add one row cluster to the original structure
    // if we must add a new row cluster n_{ej} and n_{e.} are equal to 0 elsewhere
    // we compute the original form
    if (destination_row_cluster == n_row_clusters_) {
        // the selected row cluster is the empty one
        for (int v = 0; v < num_views_; ++v) {
            for (int cc = 0; cc < n_col_clusters_[v]; cc++) {
                if (tot_t_per_col_cluster[v][cc] != 0) {
                    X_views[v] += (lambda[v][cc] / tot_t_per_col_cluster[v][cc]) *
                            (T_[v][source_row_cluster][cc] - lambda[v][cc]);
                }
                Y_views[v] += (-T_[v][source_row_cluster][cc] + lambda[v][cc]);
            }
        }
    }
    else {
        // the selected row cluster is an already existent one 
        // so we can also use the T-values for the new cluster
        for (int v = 0; v < num_views_; ++v) {
            for (int cc = 0; cc < n_col_clusters_[v]; cc++) {
                if (tot_t_per_col_cluster[v][cc] != 0) {
                    X_views[v] += (lambda[v][cc] / tot_t_per_col_cluster[v][cc]) *
                            (T_[v][source_row_cluster][cc] -
                             T_[v][destination_row_cluster][cc] - lambda[v][cc]);
                }
                Y_views[v] += (T_[v][destination_row_cluster][cc] -
                               T_[v][source_row_cluster][cc] + lambda[v][cc]);
            }
        }
    }

    double X = 0;
    double Y = 0;
    for (int k = 0; k < num_views_; ++k) {
        X += two_divided_by_n_tot_[k] * X_views[k];
        Y += ((2 * tot_lambda[k]) / n_tot_square_[k]) * Y_views[k];
    }

    double delta_tau_r = 0;

    if ((Ir_ * (Ir_ - Y)) != 0) {
        delta_tau_r = ((Ir_ * X) + (Cr_ * Y)) / (Ir_ * (Ir_ - Y));
    }

    delete X_views;
    delete Y_views;

    return delta_tau_r;
}

/**
 * Compute the deltaTauC value.
 * 
 * @param view_index, int, the considered view
 * @param tot_t_per_row_cluster, the sum of data group by row cluster
 * @param tot_lambda, the sum of lambda
 * @param lambda, array of n_row_clusters elements, the projection of the element x on the row cluster assignment
 * @param source_cluster, int, the col cluster that contains originally x
 * @param evaluated_destination_cluster, int, the col cluster in which we want to move x
 * 
 * @return double, the computed value
 */
double MultiViewCoClust::deltaTauC(int view_index, double *tot_t_per_row_cluster,
                                   double tot_lambda, double *lambda, int source_cluster,
                                   int evaluated_destination_cluster) {
    double X = 0;
    double Y = 0;

    // we check if the new label of the element x is the empty set, this means
    // that we must add one col cluster to the original structure
    // if we must add a new col cluster n_{ej} and n_{e.} are equal to 0 elsewhere
    // we compute the original form of the equation
    if (evaluated_destination_cluster == n_col_clusters_[view_index]) {
        for (int rc = 0; rc < n_row_clusters_; rc++) {
            if (tot_t_per_row_cluster[rc] != 0) {
                X += (lambda[rc] / tot_t_per_row_cluster[rc]) * (T_[view_index][rc][source_cluster] - lambda[rc]);
            }
            Y += (lambda[rc] - T_[view_index][rc][source_cluster]);
        }
    }
    else {
        for (int rc = 0; rc < n_row_clusters_; rc++) {
            if (tot_t_per_row_cluster[rc] != 0) {
                X += (lambda[rc] / tot_t_per_row_cluster[rc]) *
                        (T_[view_index][rc][source_cluster] - T_[view_index][rc][evaluated_destination_cluster] - lambda[rc]);
            }
            Y += (T_[view_index][rc][evaluated_destination_cluster] - T_[view_index][rc][source_cluster] + lambda[rc]);
        }
    }
    X = two_divided_by_n_tot_[view_index] * X;
    // X= (2/nTot) * X;
    Y = ((2 * tot_lambda) / n_tot_square_[view_index]) * Y;

    double delta_tau_c = 0;
    if ((Ic_[view_index] * (Ic_[view_index] - Y)) != 0)
        delta_tau_c = ((Ic_[view_index] * X) + (Cc_[view_index] * Y)) /
        (Ic_[view_index] * (Ic_[view_index] - Y));
    return delta_tau_c;
}