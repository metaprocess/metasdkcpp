#include "MatrixUtils.h"
#include <iostream>
#include <algorithm>

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


MatrixVariant MatrixUtils::load_from_resource_mat(const std::string& _resource_path, const std::string& _name_variable) {
    // Get the resource
    auto _res = get_resource(_resource_path);
    assertm(_res.start != nullptr, ("Could not find resource: " + _resource_path).c_str());
    assertm(_res.size > 0, ("Resource is empty: " + _resource_path).c_str());

    // Create a temporary file to write the resource data
    // Note: matio doesn't support reading directly from memory, so we need to use a temp file
    std::string temp_filename = "/tmp/matio_resource_" + _resource_path;
    std::replace(temp_filename.begin(), temp_filename.end(), '/', '_');
    
    // Write resource to temporary file
    std::ofstream temp_file(temp_filename, std::ios::binary);
    if (!temp_file.is_open()) {
        throw std::runtime_error("Could not create temporary file for MAT resource: " + temp_filename);
    }
    temp_file.write(reinterpret_cast<const char*>(_res.start), _res.size);
    temp_file.close();
    
    if (temp_file.fail()) {
        throw std::runtime_error("Failed to write resource to temporary file: " + temp_filename);
    }

    // Open the MAT file from the temporary file
    mat_t *matfp = Mat_Open(temp_filename.c_str(), MAT_ACC_RDONLY);
    if (!matfp) {
        // Clean up temp file on error
        std::remove(temp_filename.c_str());
        throw std::runtime_error("Could not open MAT file from resource: " + _resource_path);
    }

    auto _ret = load_from_mat_file(temp_filename, _name_variable);
    std::remove(temp_filename.c_str());
    return _ret;

    #if 0
    // Read the variable
    matvar_t *matvar = Mat_VarRead(matfp, _name_variable.c_str());
    assertm(matvar, ("Could not read variable '" + _name_variable + "' from MAT resource: " + _resource_path).c_str());

    // Check if it's a 2D matrix
    assertm(matvar->rank == 2, ("Variable '" + _name_variable + "' is not a 2D matrix (rank = " + std::to_string(matvar->rank) + ")").c_str());

    size_t rows = matvar->dims[0];
    size_t cols = matvar->dims[1];
    MatrixVariant result;

    try {
        // Handle different data types
        switch (matvar->class_type) {
            case MAT_C_DOUBLE: {
                if (matvar->data_type == MAT_T_DOUBLE) {
                    if (matvar->isComplex) {
                        // Complex double matrix
                        mat_complex_split_t* complex_data = static_cast<mat_complex_split_t*>(matvar->data);
                        double* real_part = static_cast<double*>(complex_data->Re);
                        double* imag_part = static_cast<double*>(complex_data->Im);
                        
                        Eigen::MatrixXcd matrix(rows, cols);
                        for (size_t i = 0; i < rows; ++i) {
                            for (size_t j = 0; j < cols; ++j) {
                                size_t idx = i + j * rows;
                                matrix(i, j) = std::complex<double>(real_part[idx], imag_part[idx]);
                            }
                        }
                        result = matrix;
                    } else {
                        // Double matrix
                        Eigen::MatrixXd matrix(rows, cols);
                        memcpy(matrix.data(), matvar->data, rows * cols * sizeof(double));
                        result = matrix;
                    }
                }
                break;
            }
            case MAT_C_SINGLE: {
                if (matvar->data_type == MAT_T_SINGLE) {
                    // Float matrix
                    Eigen::MatrixXf matrix(rows, cols);
                    memcpy(matrix.data(), matvar->data, rows * cols * sizeof(float));
                    result = matrix;
                }
                break;
            }
            case MAT_C_INT32: {
                if (matvar->data_type == MAT_T_INT32) {
                    // Int32 matrix
                    Eigen::MatrixXi matrix(rows, cols);
                    memcpy(matrix.data(), matvar->data, rows * cols * sizeof(int32_t));
                    result = matrix;
                }
                break;
            }
            case MAT_C_INT16: {
                if (matvar->data_type == MAT_T_INT16) {
                    // Int16 matrix
                    Eigen::Matrix<int16_t, Eigen::Dynamic, Eigen::Dynamic> matrix(rows, cols);
                    memcpy(matrix.data(), matvar->data, rows * cols * sizeof(int16_t));
                    result = matrix;
                }
                break;
            }
            case MAT_C_UINT16: {
                if (matvar->data_type == MAT_T_UINT16) {
                    // Uint16 matrix
                    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> matrix(rows, cols);
                    memcpy(matrix.data(), matvar->data, rows * cols * sizeof(uint16_t));
                    result = matrix;
                }
                break;
            }
            case MAT_C_INT8: {
                if (matvar->data_type == MAT_T_INT8) {
                    // Int8 matrix
                    Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic> matrix(rows, cols);
                    memcpy(matrix.data(), matvar->data, rows * cols * sizeof(int8_t));
                    result = matrix;
                }
                break;
            }
            case MAT_C_UINT8: {
                if (matvar->data_type == MAT_T_UINT8) {
                    // Uint8 matrix
                    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> matrix(rows, cols);
                    memcpy(matrix.data(), matvar->data, rows * cols * sizeof(uint8_t));
                    result = matrix;
                }
                break;
            }
            case MAT_C_UINT32: {
                if (matvar->data_type == MAT_T_UINT32) {
                    // Uint32 matrix
                    Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic> matrix(rows, cols);
                    memcpy(matrix.data(), matvar->data, rows * cols * sizeof(uint32_t));
                    result = matrix;
                }
                break;
            }
            default: {
                Mat_VarFree(matvar);
                Mat_Close(matfp);
                assertm(false, ("Unsupported MAT variable class type: " + std::to_string(matvar->class_type)).c_str());
            }
        }
    } catch (const std::exception& e) {
        Mat_VarFree(matvar);
        Mat_Close(matfp);
        throw;
    }

    // Clean up
    Mat_VarFree(matvar);
    Mat_Close(matfp);
    
    // Remove temporary file
    std::remove(temp_filename.c_str());

    return result;
    #endif
}

MatrixVariant MatrixUtils::load_mat(const std::string &_name_file, const std::string &_name_variable)
{
    assertm(!_name_file.empty(), "filename is empty");
    if(':' == _name_file[0])
    {
        auto _str = std::string(_name_file.begin() + 2, _name_file.end());
        return load_from_resource_mat(_str, _name_variable);
    }
    else
    {
        return load_from_mat_file(_name_file, _name_variable);
    }
}

MatrixVariant MatrixUtils::load_from_mat_file(const std::string &_name_file, const std::string &_name_variable)
{
    // Open the MAT file
    mat_t *matfp = Mat_Open(_name_file.c_str(), MAT_ACC_RDONLY);
    assertm(matfp, ("Could not open MAT file: " + _name_file).c_str());

    // Read the variable
    matvar_t *matvar = Mat_VarRead(matfp, _name_variable.c_str());
    assertm(matvar, ("Could not read variable '" + _name_variable + "' from MAT file: " + _name_file).c_str());

    // Check if it's a 2D matrix
    assertm(matvar->rank == 2, ("Variable '" + _name_variable + "' is not a 2D matrix (rank = " + std::to_string(matvar->rank) + ")").c_str());

    size_t rows = matvar->dims[0];
    size_t cols = matvar->dims[1];
    MatrixVariant result;

    try {
        // Handle different data types
        switch (matvar->class_type) {
            case MAT_C_DOUBLE: {
                if (matvar->data_type == MAT_T_DOUBLE) {
                    if (matvar->isComplex) {
                        // Complex double matrix
                        mat_complex_split_t* complex_data = static_cast<mat_complex_split_t*>(matvar->data);
                        double* real_part = static_cast<double*>(complex_data->Re);
                        double* imag_part = static_cast<double*>(complex_data->Im);
                        
                        Eigen::MatrixXcd matrix(rows, cols);
                        for (size_t i = 0; i < rows; ++i) {
                            for (size_t j = 0; j < cols; ++j) {
                                size_t idx = i + j * rows;
                                matrix(i, j) = std::complex<double>(real_part[idx], imag_part[idx]);
                            }
                        }
                        result = matrix;
                    } else {
                        // Double matrix
                        Eigen::MatrixXd matrix(rows, cols);
                        memcpy(matrix.data(), matvar->data, rows * cols * sizeof(double));
                        result = matrix;
                    }
                }
                break;
            }
            case MAT_C_SINGLE: {
                if (matvar->data_type == MAT_T_SINGLE) {
                    // Float matrix
                    Eigen::MatrixXf matrix(rows, cols);
                    memcpy(matrix.data(), matvar->data, rows * cols * sizeof(float));
                    result = matrix;
                }
                break;
            }
            case MAT_C_INT32: {
                if (matvar->data_type == MAT_T_INT32) {
                    // Int32 matrix
                    Eigen::MatrixXi matrix(rows, cols);
                    memcpy(matrix.data(), matvar->data, rows * cols * sizeof(int32_t));
                    result = matrix;
                }
                break;
            }
            case MAT_C_INT16: {
                if (matvar->data_type == MAT_T_INT16) {
                    // Int16 matrix
                    Eigen::Matrix<int16_t, Eigen::Dynamic, Eigen::Dynamic> matrix(rows, cols);
                    memcpy(matrix.data(), matvar->data, rows * cols * sizeof(int16_t));
                    result = matrix;
                }
                break;
            }
            case MAT_C_UINT16: {
                if (matvar->data_type == MAT_T_UINT16) {
                    // Uint16 matrix
                    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> matrix(rows, cols);
                    memcpy(matrix.data(), matvar->data, rows * cols * sizeof(uint16_t));
                    result = matrix;
                }
                break;
            }
            case MAT_C_INT8: {
                if (matvar->data_type == MAT_T_INT8) {
                    // Int8 matrix
                    Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic> matrix(rows, cols);
                    memcpy(matrix.data(), matvar->data, rows * cols * sizeof(int8_t));
                    result = matrix;
                }
                break;
            }
            case MAT_C_UINT8: {
                if (matvar->data_type == MAT_T_UINT8) {
                    // Uint8 matrix
                    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> matrix(rows, cols);
                    memcpy(matrix.data(), matvar->data, rows * cols * sizeof(uint8_t));
                    result = matrix;
                }
                break;
            }
            case MAT_C_UINT32: {
                if (matvar->data_type == MAT_T_UINT32) {
                    // Uint32 matrix
                    Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic> matrix(rows, cols);
                    memcpy(matrix.data(), matvar->data, rows * cols * sizeof(uint32_t));
                    result = matrix;
                }
                break;
            }
            default: {
                Mat_VarFree(matvar);
                Mat_Close(matfp);
                assertm(false, ("Unsupported MAT variable class type: " + std::to_string(matvar->class_type)).c_str());
            }
        }
    } catch (const std::exception& e) {
        Mat_VarFree(matvar);
        Mat_Close(matfp);
        throw;
    }

    // Clean up
    Mat_VarFree(matvar);
    Mat_Close(matfp);

    return result;
}
