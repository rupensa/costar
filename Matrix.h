#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cstdlib>
#include <fstream>

class Matrix {
private:
    double **matrix;
    int nrow;
    int ncol;

    /**
     * Initialize a new empty double matrix 
     * @param row int, the number of rows
     * @param col int, the number of columns 
     */
    void initMatrix(int row, int col) {
        srand(time(NULL));
        this->nrow = row;
        this->ncol = col;
        this->matrix = new double *[row];
        for (int i = 0; i < nrow; i++) matrix[i] = new double[ncol];
    }
    
    /**
     * Fills the current matrix with random numbers from 1 to 100. 
     */
    void randMatrix() {
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) this->matrix[i][j] = rand() % 100 + 1;
        }
    }

public:

    Matrix(double **matrix_a, int nr, int nc)
    : matrix(matrix_a), nrow(nr), ncol(nc) {
    }
    
    /**
     * Create a matrix from file.
     * 
     * @param row, int number of rows 
     * @param col, int number of columns 
     * @param inFile, char* the path of the input file 
     */
    Matrix(int row, int col, const char *input_file) {
        initMatrix(row, col);
        this->nrow = row;
        this->ncol = col;
        std::ifstream in(input_file);
        while (in) {
            for (int i = 0; i < row; i++) {
                for (int j = 0; j < col; j++) in >> matrix[i][j];
            }
        }
        in.close();
    }
    
    /**
     * Create a matrix of random numbers. 
     * 
     * @param row int, the number of rows 
     * @param col int, the number of columns 
     */
    Matrix(int row, int col) {
        initMatrix(row, col);
        randMatrix();
    }
    
    /**
     * Create a matrix cloning the given matrix. 
     * 
     * @param matrixToClone the matrix to clone
     */
    Matrix(const Matrix &to_clone) {
        this->nrow = to_clone.nrow;
        this->ncol = to_clone.ncol;
        this->matrix = new double *[nrow];
        for (int i = 0; i < this->nrow; i++) {
            this->matrix[i] = new double[ncol];
            for (int j = 0; j < this->ncol; j++) {
                this->matrix[i][j] = to_clone.matrix[i][j];
            }
        }
    }

    double **getDistributionMatrix(int *row_assignment, int *colAssignment, int n_row_clusters, int nColClusters);

    double *getRowGroupByColClust(int row_index, int *colAssignment, int n_column_clusters);

    double *getColGroupByRowClust(int col_index, int *row_assignment, int n_row_clusters);
    
    /**
     * Get the number of rows 
     * @return the number of rows 
     */
    int n_row() {
        return nrow;
    }
    
    /**
     * Get the number of columns 
     * @return the number of columns 
     */
    int n_col() {
        return ncol;
    }

    /**
     * Clone the matrix and returns the copy
     * @return the clone matrix
     */
    double **cloneMatrix() {
        double **ris = new double *[nrow];
        for (int i = 0; i < nrow; i++) {
            ris[i] = new double[ncol];
            for (int j = 0; j < ncol; j++) ris[i][j] = matrix[i][j];
        }
        return ris;
    }
    
    /**
     * Edit an element of the matrix 
     * @param r, the row index
     * @param c, the col index 
     * @param value, the value to set
     */
    void setElem(int r, int c, double value) {
        matrix[r][c] = value;
    }

    /**
     * Get an element of the matrix
     * @param r, the row index 
     * @param c, the col index 
     * @return the element contained in the cell <r,c>
     */
    double getElem(int r, int c) {
        return matrix[r][c];
    }

    void print();
    
    /**
     * Delete the matrix
     */
    ~Matrix() {
        for (int i = 0; i < nrow; i++) {
            delete[] matrix[i];
        }
        delete matrix;
    }
};

#endif