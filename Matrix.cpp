#include "Matrix.h"

/**
 * Computes the contingency matrix for clusters. Each cell <i,j> contains the 
 * sum of the cells in the original matrix that are assigned to the row cluster i
 * and the column cluster j. 
 * 
 * @param row_clusters_assignment : list, a cluster assignment for rows
 * @param col_clusters_assignment : list, a cluster assignment for columns
 * @param n_row_clusters : int, the number of row clusters
 * @param n_col_clusters : int, the number of column clusters
 * @return the computed contingency matrix
 */
double **Matrix::getDistributionMatrix(int *row_clusters_assignment, 
                                       int *col_clusters_assignment, 
                                       int n_row_clusters, int n_col_clusters) {
    double **T = new double *[n_row_clusters];
    for (int i = 0; i < n_row_clusters; i++) {
        T[i] = new double[n_col_clusters];
        for (int j = 0; j < n_col_clusters; j++) T[i][j] = 0;
    }
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            T[row_clusters_assignment[i]][col_clusters_assignment[j]] += matrix[i][j];
        }
    }
    return T;
}

/**
 * Fixed a document (row), computes for each column cluster the sum of the values (in the dataset) 
 *      of the features associated to that cluster.
 * 
 * @param row_index : int, the row to consider
 * @param col_assignment : list reference, the cluster assignment for columns
 * @param n_column_clusters : int, number of column cluster
 * @return a list containing for each column cluster the sum of features values for the document 'row'
 */
double *Matrix::getRowGroupByColClust(int row_index, int *col_assignment, int n_column_clusters) {
    double *ris = new double[n_column_clusters];
    for (int i = 0; i < n_column_clusters; i++) ris[i] = 0;
    for (int j = 0; j < ncol; j++) {
        ris[col_assignment[j]] += matrix[row_index][j];
    }
    return ris;
}

/**
 * Fixed a feature (column), computes for each row cluster the sum 
 *    of the values (in the dataset) of the rows associated to that cluster.
 * 
 * @param col_index : int, the column to consider
 * @param row_assignment : list reference, the document-cluster assignment for rows
 * @param n_row_clusters : int, number of row cluster
 * @return a list containing for each row cluster the sum of document values for the feature 'col'
 */
double *Matrix::getColGroupByRowClust(int col_index, int *row_assignment, int n_row_clusters) {
    double *ris = new double[n_row_clusters];
    for (int i = 0; i < n_row_clusters; i++) ris[i] = 0;
    for (int j = 0; j < nrow; j++) {
        ris[row_assignment[j]] += matrix[j][col_index];
    }
    return ris;
}

/**
 * Print the matrix.
 */
void Matrix::print() {
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}
