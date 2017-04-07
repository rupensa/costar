#include "Matrix.h"
#include <unistd.h>
#include <vector>
#include <string>
#include "MultiViewCoClust.h"

int main(int argc, char **argv) {
    int c;
    int nrow;
    int n_iter;
    int min_num_row_cluster = 2;
    int random_partition = 0;

    char *configFile = NULL;
    char *outputFile = NULL;
    opterr = 0;

    srand(time(NULL));

    if (argc == 1) {
        std::cout << "command -r nrow -n n_iter -i configFile -o outputFile [-p] "
                "[-k min_num_row_cluster]"
                << std::endl;
        return 1;
    }

    // parse arguments
    while ((c = getopt(argc, argv, "pr:i:o:n:k:")) != -1) switch (c) {
        case 'p':
            random_partition = 1;
            break;
        case 'r':
            nrow = atoi(optarg);
            break;
        case 'i':
            configFile = optarg;
            break;
        case 'o':
            outputFile = optarg;
            break;
        case 'n':
            n_iter = atoi(optarg);
            break;
        case 'k':
            min_num_row_cluster = atoi(optarg);
            break;
        case '?':
            std::cout << "command -r nrow -n n_iter -i configFile -o outputFile "
                    "[-p] [-k min_num_row_cluster]"
                    << std::endl;
            return 1;
        default:
            std::cout << "command -r nrow -n n_iter -i configFile -o outputFile "
                    "[-p] [-k min_num_row_cluster]"
                    << std::endl;
            return 1;
        }
    if (nrow == 0 || n_iter == 0) {
        std::cout << "command -r nrow -n n_iter -i configFile -o outputFile [-p] "
                "[-k min_num_row_cluster]"
                << std::endl;
        return 1;
    }

    // read the input files
    std::vector<std::string> feature_space;
    std::vector<int> col_feature_space;
    std::ifstream in(configFile);
    std::string file;
    int n_col = 0;
    while (in >> n_col) {
        in >> file;
        feature_space.push_back(file);
        col_feature_space.push_back(n_col);
    }

    // create view matrices
    std::vector<Matrix> matrix_list;
    for (int i = 0; i < feature_space.size(); ++i) {
        Matrix a(nrow, col_feature_space[i], feature_space[i].c_str());
        matrix_list.push_back(a);
    }

    std::string outputName = outputFile;
    // parameters of the constructor are:
    // * the matrices list,
    // * the minimum number of clusters to create on rows,
    // * the prefix for the output file names.
    MultiViewCoClust *b =
            new MultiViewCoClust(matrix_list, min_num_row_cluster, outputName);
    b->buildCoClu(random_partition, n_iter);

    delete b;
    return 0;
}
