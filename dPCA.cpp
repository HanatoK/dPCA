#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <cmath>
#include <fmt/format.h>

#define USE_EIGEN3
#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif

#define _USE_MATH_DEFINES

using std::vector;
using std::string;
using std::cout;
using std::endl;

void splitString(const std::string& data, const std::string& delim, std::vector<std::string>& dest) {
    if (delim.empty()) {
        dest.push_back(data);
        return;
    }
    size_t index = 0, new_index = 0;
    std::string tmpstr;
    while (index != data.length()) {
        new_index = data.find(delim, index);
        if (new_index != std::string::npos) tmpstr = data.substr(index, new_index - index);
        else tmpstr = data.substr(index, data.length());
        if (!tmpstr.empty()) {
            dest.push_back(tmpstr);
        }
        if (new_index == std::string::npos) break;
        index = new_index + 1;
    }
}

class DihedralData {
public:
    struct dPCAResult {
        vector<double> eigenvalues;
        vector<vector<double>> eigenvectors;
        vector<double> fractions;
    };
    DihedralData();
    DihedralData(const string& filename);
    void readFromFile(const string& filename);
    void writeRawToFile(const string& filename) const;
    vector<vector<double>> buildCovarianceMatrix() const;
    dPCAResult PCA() const;
    void writeProjection(const dPCAResult& result, const string& filename) const;
private:
    vector<vector<double>> m_sine;
    vector<vector<double>> m_cosine;
};

DihedralData::DihedralData() {}

DihedralData::DihedralData(const string& filename) {
    readFromFile(filename);
}

void DihedralData::readFromFile(const string& filename) {
    using std::ifstream;
    ifstream ifs(filename.c_str());
    string line;
    vector<string> tmp_fields;
    vector<double> tmp_sine;
    vector<double> tmp_cosine;
    while (std::getline(ifs, line)) {
        tmp_fields.clear();
        tmp_sine.clear();
        tmp_cosine.clear();
        splitString(line, " ", tmp_fields);
        // read the dihedral angles in degree and transform them into radians
        for (size_t i = 0; i < tmp_fields.size(); ++i) {
            tmp_cosine.push_back(std::cos(std::stod(tmp_fields[i]) / 180.0 * M_PI));
            tmp_sine.push_back(std::sin(std::stod(tmp_fields[i]) / 180.0 * M_PI));
        }
        m_sine.push_back(tmp_sine);
        m_cosine.push_back(tmp_cosine);
    }
}

void DihedralData::writeRawToFile(const string& filename) const {
    using std::ofstream;
    ofstream ofs(filename.c_str());
    if (m_sine.size() == 0) return;
    // get the number of variables
    const size_t N = m_sine[0].size();
    // the number of samples
    const size_t M = m_sine.size();
    string line;
    for (size_t i = 0; i < M; ++i) {
        line.clear();
        for (size_t j = 0; j < N; ++j) {
            line += fmt::format(" {:15.10f} {:15.10f}", m_cosine[i][j], m_sine[i][j]);
        }
        ofs << line << "\n";
    }
}

vector<vector<double>> DihedralData::buildCovarianceMatrix() const {
    if (m_sine.size() == 0) return vector<vector<double>>();
    // get the number of variables
    const size_t N = m_sine[0].size();
    // the number of samples
    const size_t M = m_sine.size();
    // compute the averages of variables
    vector<double> m_sine_mean(N, 0);
    vector<double> m_cosine_mean(N, 0);
    for (size_t j = 0; j < N; ++j) {
        for (size_t i = 0; i < M; ++i) {
            m_sine_mean[j] += m_sine[i][j];
            m_cosine_mean[j] += m_cosine[i][j];
        }
        m_sine_mean[j] /= M;
        m_cosine_mean[j] /= M;
    }
    // compute the covariance
    vector<vector<double>> covariance_matrix(N * 2, vector<double>(N * 2, 0));
    for (size_t j = 0; j < N; ++j) {
        for (size_t k = 0; k < N; ++k) {
            for (size_t i = 0; i < M; ++i) {
                const double left_even = m_cosine[i][j] - m_cosine_mean[j];
                const double right_even = m_cosine[i][k] - m_cosine_mean[k];
                const double left_odd = m_sine[i][j] - m_sine_mean[j];
                const double right_odd = m_sine[i][k] - m_sine_mean[k];
                covariance_matrix[2*j][2*k] += left_even * right_even;
                covariance_matrix[2*j+1][2*k] += left_odd * right_even;
                covariance_matrix[2*j][2*k+1] += left_even * right_odd;
                covariance_matrix[2*j+1][2*k+1] += left_odd * right_odd;
            }
            covariance_matrix[2*j][2*k] /= M;
            covariance_matrix[2*j+1][2*k] /= M;
            covariance_matrix[2*j][2*k+1] /= M;
            covariance_matrix[2*j+1][2*k+1] /= M;
        }
    }
    return covariance_matrix;
}

DihedralData::dPCAResult DihedralData::PCA() const {
    vector<vector<double>> covariance_matrix = buildCovarianceMatrix();
    dPCAResult result;
    if (covariance_matrix.empty()) return result;
#ifdef USE_EIGEN3
    using namespace Eigen;
    MatrixXf M(covariance_matrix.size(), covariance_matrix[0].size());
    for (size_t i = 0; i < covariance_matrix.size(); ++i) {
        for (size_t j = 0; j < covariance_matrix[i].size(); ++j) {
            M(i, j) = covariance_matrix[i][j];
        }
    }
    SelfAdjointEigenSolver<MatrixXf> eigensolver(M);
    if (eigensolver.info() != Success) abort();
    cout << "Eigenvalues:\n" << eigensolver.eigenvalues() << endl;
    cout << "Eigenvectors:\n" << eigensolver.eigenvectors() << endl;
    const auto tmp_eigenvalues = eigensolver.eigenvalues();
    const auto tmp_eigenvectors = eigensolver.eigenvectors();
    double factor = 0;
    for (size_t i = 0; i < covariance_matrix.size(); ++i) {
        factor += tmp_eigenvalues(i);
        result.eigenvalues.push_back(tmp_eigenvalues(i));
        vector<double> tmp_vector;
        for (size_t j = 0; j < covariance_matrix.size(); ++j) {
            tmp_vector.push_back(tmp_eigenvectors(i, j));
        }
        result.eigenvectors.push_back(tmp_vector);
    }
    cout << "Eigenvectors (converted to std::vector):\n";
    for (size_t i = 0; i < covariance_matrix.size(); ++i) {
        for (size_t j = 0; j < covariance_matrix.size(); ++j) {
            cout << result.eigenvectors[i][j] << " ";
        }
        cout << endl;
    }
    cout << "Fractions:\n";
    for (size_t i = 0; i < covariance_matrix.size(); ++i) {
        result.fractions.push_back(result.eigenvalues[i] / factor);
        cout << result.fractions.back() << endl;
    }
#endif
    return result;
}

void DihedralData::writeProjection(const dPCAResult& result, const string& filename) const {
    using std::ofstream;
    ofstream ofs(filename.c_str());
    const size_t num_components = result.eigenvalues.size();
    // the number of samples
    const size_t M = m_sine.size();
    // get the number of variables
    const size_t N = m_sine[0].size();
    // write the vectors
    ofs << "# Eigenvectors:\n";
    for (size_t i = 0; i < num_components; ++i) {
        ofs << fmt::format("# PC{} ", i + 1);
        for (size_t j = 0; j < N * 2; ++j) {
            ofs << fmt::format(" {:15.10f}", result.eigenvectors[j][i]);
        }
        ofs << "\n";
    }
    // write the eigenvalues
    ofs << "# Eigenvalues:\n";
    for (size_t i = 0; i < num_components; ++i) {
        ofs << fmt::format("# PC{} {:15.10f}\n", i + 1, result.eigenvalues[i]);
    }
    // write the fractions
    ofs << "# Component fractions:\n";
    for (size_t i = 0; i < num_components; ++i) {
        ofs << fmt::format("# PC{} {:15.10f}\n", i + 1, result.fractions[i]);
    }
    // write the projections
    ofs << "# Projections:\n";
    ofs << "# ";
    for (size_t i = 0; i < num_components; ++i) {
        ofs << fmt::format("PC{} ", i + 1);
    }
    ofs << "\n";
    for (size_t j = 0; j < M; ++j) {
        for (size_t l = 0; l < num_components; ++l) {
            // dot product (projection)
            double product = 0;
            for (size_t k = 0; k < N; ++k) {
                product += result.eigenvectors[2*k][l] * m_cosine[j][k];
                product += result.eigenvectors[2*k+1][l] * m_sine[j][k];
            }
            ofs << fmt::format(" {:15.10f}", product);
        }
        ofs << "\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) return 1;
    DihedralData x(argv[1]);
    x.writeRawToFile(string(argv[2]) + ".transformed");
    DihedralData::dPCAResult result = x.PCA();
    x.writeProjection(result, string(argv[2]) + ".projected");
    return 0;
}
