#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <unistd.h>
#include <math.h>

using namespace std;

double rowSum(double **matrix, int num_col, int row) {
    double sum = 0.0;
    for (int i = 0; i < num_col; ++i) {
        sum += matrix[row][i];
    }
    return sum;
}

double colSum(double **matrix, int num_row, int col) {
    double sum = 0.0;
    for (int i = 0; i < num_row; ++i) {
        sum += matrix[i][col];
    }
    return sum;
}

double nmi(double **matrix, int n_row, int n_col, int num_points) {
    double ris = 0.0;
    for (int i = 0; i < n_row; ++i) {
        for (int j = 0; j < n_col; ++j) {
            double arg = matrix[i][j] * num_points /
                    (colSum(matrix, n_row, j) * rowSum(matrix, n_col, i));
            if (arg != 0) ris += matrix[i][j] * log(arg);
        }
    }
    double den1 = 0.0;
    for (int i = 0; i < n_row; ++i) {
        den1 +=
                rowSum(matrix, n_col, i) * log(rowSum(matrix, n_col, i) / num_points);
    }
    double den2 = 0.0;
    for (int i = 0; i < n_col; ++i) {
        den2 +=
                colSum(matrix, n_row, i) * log(colSum(matrix, n_row, i) / num_points);
    }
    return ris / sqrt(den1 * den2);
}

double ari(vector<int> &first, vector<int> &second) {
    // the number of pairs of elementsthat are in the same class I and in the same
    // cluster J
    double a = 0;
    // the number of pairs of elements that are in different class I and in
    // different cluster J
    double b = 0;
    // the number of pairs of elements that are in the same class I and in
    // different cluster J
    double c = 0;
    // the number of pairs of elements that are in different class I and in the
    // same cluster J
    double d = 0;
    // ARI = 2(ab - cd) / [(a+d)(d+b) + (a+c)(c+b)]
    vector<int>::iterator it1, it2;
    int i, j;
    for (int i = 0; i < first.size(); ++i) {
        for (int j = 0; j < second.size(); ++j) {
            if (i != j) {
                if (first[i] == first[j] && second[i] == second[j]) a++;
                if (first[i] != first[j] && second[i] != second[j]) b++;
                if (first[i] == first[j] && second[i] != second[j]) c++;
                if (first[i] != first[j] && second[i] == second[j]) d++;
            }
        }
    }
    double num = 2 * (a * b - c * d);
    double den = ((a + d) * (d + b) + (a + c) * (c + b));
    return num / den;
}

double **generateMatrix(vector<int> &first, vector<int> &second, int n_row,
                        int n_col) {
    vector<int>::iterator it1, it2;
    double **matrix = new double *[n_row];
    for (int i = 0; i < n_row; ++i) {
        matrix[i] = new double[n_col];
        for (int j = 0; j < n_col; ++j) matrix[i][j] = 0.0;
    }

    for (it1 = first.begin(), it2 = second.begin(); it1 != first.end();
            ++it1, ++it2) {
        matrix[*it1][*it2] += 1;
    }

    return matrix;
}

int main(int argc, char *argv[]) {
    string dato;
    string line;
    vector<int> first_assignment;
    vector<int> second_assignment;

    ifstream inF(argv[1]);

    /*
        int first = 0;
        while ( inF >> first){
                firstAssign.push_back(first);
        }
     */
    getline(inF, line);
    istringstream issF(line);
    while (getline(issF, dato, ' ')) {
        int first = atoi(dato.c_str());
        first_assignment.push_back(first);
    }

    inF.close();

    ifstream inS(argv[2]);
    int second = 0;
    while (inS >> second) {
        second_assignment.push_back(second);
    }

    inS.close();
    
    set<int> first_uniq;
    set<int> second_uniq;
    vector<int>::iterator it1, it2;
    for (it1 = first_assignment.begin(), it2 = second_assignment.begin();
            it1 != first_assignment.end(); ++it1, ++it2) {
        first_uniq.insert(*it1);
        second_uniq.insert(*it2);
    }
    int n_row = first_uniq.size();
    int n_col = second_uniq.size();
    double **matrix = generateMatrix(first_assignment, second_assignment, n_row, n_col);
    for (int i = 0; i < n_row; i++) {
        for (int j = 0; j < n_col; j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }

    cout << nmi(matrix, n_row, n_col, first_assignment.size()) << endl;
    cout << ari(first_assignment, second_assignment) << endl;
    return 0;
}