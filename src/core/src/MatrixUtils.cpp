#include "MatrixUtils.h"
#include <iostream>

void MatrixUtils::test_eigen_memory_order()
{
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rowMajorMatrix(3, 3);
    Eigen::MatrixXd colMajorMatrix(3, 3);
    const double _buffer[] = {
        1.1, 2.2, 3.345,
        4.1, 5.2, 6.345,
        7.1, 8.2, 9.345
    };
    memcpy(rowMajorMatrix.data(), _buffer, sizeof(_buffer));
    memcpy(colMajorMatrix.data(), _buffer, sizeof(_buffer));

    // Print the row-major matrix
    std::cout << "Row-Major Matrix:\n" << rowMajorMatrix << std::endl;
    std::cout << "Col-Major Matrix:\n" << colMajorMatrix << std::endl;

    std::cout << "row-Major memory::\n";
    for(int i = 0; i < 9; ++i) {
        std::cout << rowMajorMatrix.data()[i] << " ";
    }
    std::cout << "\n";

    std::cout << "Col-Major memory::\n";
    for(int i = 0; i < 9; ++i) {
        std::cout << colMajorMatrix.data()[i] << " ";
    }
    std::cout << "\n";

}

void MatrixUtils::test_save_and_load()
{
#if 0
    Eigen::MatrixXcd _matrix(2, 3);
    _matrix << std::complex<double>(1., 1.), std::complex<double>(1., 2.), std::complex<double>(1., 3.), std::complex<double>(2., 1.), std::complex<double>(2., 2.), std::complex<double>(2., 3.);
    const auto _filename = "../tests/matrix_2by3_cd.bin";
    save_to_file(_filename, _matrix);

    std::cout << "Mat ref = \n" << _matrix << std::endl;
    
    auto _mat_loaded = std::get<Eigen::MatrixXcd>(load_from_file(_filename));
    std::cout << "Mat loaded = \n" << _mat_loaded << std::endl;
    
    assertm(0 == memcmp(_mat_loaded.data(), _matrix.data(), _matrix.size() * sizeof(_matrix(0, 0))), "save and load mismatch");
#endif
}

Eigen::VectorXd MatrixUtils::linspace(double start, double end, int samples)
{
    Eigen::VectorXd vec(samples);
    vec.setLinSpaced(start, end);
    return std::move(vec);
}


MatrixVariant MatrixUtils::load_from_file(const std::string& _name_file) {
    assertm(!_name_file.empty(), "filename is empty");
    if(':' == _name_file[0])
    {
        auto _res = get_resource(std::string(_name_file.begin() + 2, _name_file.end()));
        return load_from_resource(_res);
    }
    else
    {
        return load_from_file_on_disk(_name_file);
    }
}

MatrixVariant MatrixUtils::load_from_resource(const Resource& _resource) {
    // Read the header: rows, columns, and type code
    int rows, cols, type_code;
    rows = *reinterpret_cast<const int*>(_resource.start);
    cols = *reinterpret_cast<const int*>(_resource.start + sizeof(int));
    type_code = *reinterpret_cast<const int*>(_resource.start + 2 * sizeof(int));
    const auto _ptr_matrix = reinterpret_cast<const uint8_t*>(_resource.start + 3 * sizeof(int));

    // Select the matrix type based on type_code and load the data
    switch(type_code) {
        case 0:
            return load_matrix_from_buffer<uint8_t>(_ptr_matrix, rows, cols);
        case 1:
            return load_matrix_from_buffer<int8_t>(_ptr_matrix, rows, cols);
        case 2:
            return load_matrix_from_buffer<uint16_t>(_ptr_matrix, rows, cols);
        case 3:
            return load_matrix_from_buffer<int16_t>(_ptr_matrix, rows, cols);
        case 4:
            return load_matrix_from_buffer<int32_t>(_ptr_matrix, rows, cols);
        case 5:
            return load_matrix_from_buffer<float>(_ptr_matrix, rows, cols);
        case 6:
            return load_matrix_from_buffer<double>(_ptr_matrix, rows, cols);
        case 14:
            return load_matrix_from_buffer<std::complex<double>>(_ptr_matrix, rows, cols);
        default:
            throw std::runtime_error("Unsupported type code: " + std::to_string(type_code));
    }
}

MatrixVariant MatrixUtils::load_from_file_on_disk(const std::string& _name_file) {
    assertm(!_name_file.empty(), "filename is empty");
    // Open the file in binary mode
    std::ifstream file(_name_file, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + _name_file);
    }

    // Read the header: rows, columns, and type code
    int rows, cols, type_code;
    file.read(reinterpret_cast<char*>(&rows), sizeof(int));
    file.read(reinterpret_cast<char*>(&cols), sizeof(int));
    file.read(reinterpret_cast<char*>(&type_code), sizeof(int));
    if(!file) {
        throw std::runtime_error("Failed to read header from file: " + _name_file);
    }

    // Select the matrix type based on type_code and load the data
    switch(type_code) {
        case 0:
            return load_matrix<uint8_t>(file, rows, cols);
        case 1:
            return load_matrix<int8_t>(file, rows, cols);
        case 2:
            return load_matrix<uint16_t>(file, rows, cols);
        case 3:
            return load_matrix<int16_t>(file, rows, cols);
        case 4:
            return load_matrix<int32_t>(file, rows, cols);
        case 5:
            return load_matrix<float>(file, rows, cols);
        case 6:
            return load_matrix<double>(file, rows, cols);
        case 14:
            return load_matrix<std::complex<double>>(file, rows, cols);
        default:
            throw std::runtime_error("Unsupported type code: " + std::to_string(type_code));
    }
}

MatrixCompareResult::MatrixCompareResult(const bool &_match, const double &_sum_diff, const double &_max_diff, const double& _min_diff):
    m_match(_match), m_sum_diff(_sum_diff), m_max_diff(_max_diff), m_min_diff(_min_diff)
{
}
