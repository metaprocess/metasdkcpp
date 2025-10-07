#ifndef __MATRIX_UTILS_H__
#define __MATRIX_UTILS_H__

#include "Eigen/Dense"
#include <unsupported/Eigen/CXX11/Tensor>
#include <fstream>
#include <iostream>
#include <CommonConstants.h>
#include <variant>
#include <limits>
#include <stdexcept>      // For std::runtime_error
#include <matio.h>
#include RESOURCES_HEADER


// Row-major matrix
// using TypeRowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
// // Column-major matrix (default for MatrixXd)
// using TypeColMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

// using TypeRowMatrixXComplex = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
class MatrixCompareResult
{
public:
    MatrixCompareResult(const bool& _match, const double& _sum_diff, const double& _max_diff, const double& _min_diff = 0);
    // friend std::ostream& operator<<(std::ostream& output, const MatrixCompareResult& _result);
    friend std::ostream& operator<<(std::ostream& output, const MatrixCompareResult& _result)
    {
        output << "Compare Result is " << (_result.match() ? "(MATCH)." : "(!UNMATCH!), ");
        if(!_result.match())
        {
            output << "max diff = (" << _result.max_diff() << "), sum diff = (" << _result.sum_diff() << "), min diff = (" << _result.min_diff() << ")";
        }
        return output;
    }
    bool match() const {return m_match;}
    double max_diff() const {return m_max_diff;}
    double sum_diff() const {return m_sum_diff;}
    double min_diff() const {return m_min_diff;}
private:
    bool m_match{0};
    double m_sum_diff{0};
    double m_max_diff{0};
    double m_min_diff{0};
};


#define DEG2RAD(A) (A*M_PI/180.)
#define RAD2DEG(A) (A*180./M_PI)

using MatrixVariant = std::variant<
    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>,
    Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic>,
    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>,
    Eigen::Matrix<int16_t, Eigen::Dynamic, Eigen::Dynamic>,
    Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic>,
    Eigen::MatrixXi,
    Eigen::MatrixXf,
    Eigen::MatrixXd,
    Eigen::MatrixXcd
>;

class MatrixUtils
{
public:
    MatrixUtils(){}
    ~MatrixUtils(){}

    static void test_eigen_memory_order();
    static void test_save_and_load();
    
    template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    static void save_to_file(const std::string& _name_file, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> _matrix);

    template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    static void release(Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& _matrix);
    
    // template<class TYPE>
    // static Eigen::Matrix<TYPE, Eigen::Dynamic, Eigen::Dynamic> load_from_file(const std::string& _name_file, const int& _rows, const int& _cols);

    static MatrixVariant load_from_file(const std::string& _name_file);
    static MatrixVariant load_from_resource(const Resource& _resource);
    static MatrixVariant load_from_file_on_disk(const std::string& _name_file);
    
    template<typename TYPE>
    static TYPE load_value_from_file(const std::string& _name_file);
    
    template<typename TYPE>
    static TYPE load_value_from_file_on_disk(const std::string& _name_file);
    
    template<typename TYPE>
    static TYPE load_value_from_resource(const Resource& _resource);
    
    static Eigen::VectorXd linspace(double start, double end, int samples);

    template<typename _Scalar>
    static MatrixCompareResult compare_matrices(const Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>& _mat1, const Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>& _mat2);

    template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    static void save_to_mat_file(const std::string& filename, const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& matrix, const std::string& varname = "matrix");

    template<typename Scalar, int NumDims>
    static void save_to_mat_file(const std::string& filename, const Eigen::Tensor<Scalar, NumDims>& tensor, const std::string& varname = "tensor");

    static MatrixVariant load_mat(const std::string& _name_file, const std::string& _name_variable);
    
private:
    static MatrixVariant load_from_mat_file(const std::string& _name_file, const std::string& _name_variable);
    static MatrixVariant load_from_resource_mat(const std::string& _resource_path, const std::string& _name_variable);
    template <typename Scalar>
    static Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> load_matrix(std::istream& file, int rows, int cols);
    
    template <typename Scalar>
    static Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> load_matrix_from_buffer(const void* _stream, int rows, int cols);
};

template<typename _Scalar>
MatrixCompareResult MatrixUtils::compare_matrices(const Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>& _mat1, const Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>& _mat2)
{
    // Step 1: Check if dimensions match
    if (_mat1.rows() != _mat2.rows() || _mat1.cols() != _mat2.cols()) {
        std::cout << "Matrices have different dimensions: (" << _mat1.rows() << "x" << _mat1.cols() 
                << ") vs (" << _mat2.rows() << "x" << _mat2.cols() << ")" << std::endl;
        return MatrixCompareResult(false, -1, -1);
    }

    // Step 2: Handle empty matrices
    if (_mat1.size() == 0) {
        // Both matrices are empty and have the same dimensions
        std::cout << "Matrices are exactly EQUAL\n";
        return MatrixCompareResult(true, 0, 0);
    }

    // Step 3: Check if matrices are equal and return the result
    bool are_equal = (_mat1 == _mat2);
    if(are_equal)
    {
        return MatrixCompareResult(true, 0, 0);
    }
    // Step 4: Compute the absolute difference matrix
    const _Scalar* _ptr1 = _mat1.data();
    const _Scalar* _ptr2 = _mat2.data();
    auto _size = _mat1.size();
    auto _idx = _mat1.size();
    auto _idx_max = _mat1.size();
    auto _idx_min = _mat1.size();
    _idx_max = 0;
    _idx_min = 0;
    double _max_diff{0};
    double _min_diff{std::numeric_limits<double>::max()};
    double _sum_diff{0};
    double _diff{0};
    for(_idx = 0; _idx < _size; _idx++)
    {
        _diff = std::abs(_ptr1[_idx] - _ptr2[_idx]);
        _sum_diff += _diff;
        if(_diff > _max_diff)
        {
            _max_diff = _diff;
        }
        else if(_diff < _min_diff)
        {
            _min_diff = _diff;
        }
    }
    // auto diff = (_mat1 - _mat2).cwiseAbs();

    // for(int _r = 0; _r < 10; _r++)
    // {
    //     for(int _c = 0; _c < 5; _c++)
    //     {
    //         std::cout << diff(_r, _c) << " ";
    //     }
    //     std::cout << "\n";
    // }

    // Step 5: Find maximum and minimum differences and their positions
    // Eigen::Index max_row, max_col, min_row, min_col;
    // auto max_diff = diff.maxCoeff(&max_row, &max_col);
    // auto min_diff = diff.minCoeff(&min_row, &min_col);

    // Step 6: Print the results to stdout
    // std::cout << "Maximum difference: " << max_diff << " at (" << max_row << ", " << max_col << ")" << std::endl;
    // std::cout << "Minimum difference: " << min_diff << " at (" << min_row << ", " << min_col << ")" << std::endl;

    return MatrixCompareResult(false, _sum_diff, _max_diff, _min_diff);
}

template<typename TYPE>
TYPE MatrixUtils::load_value_from_file(const std::string &_name_file)
{
    assertm(!_name_file.empty(), "filename is empty");
    if(':' == _name_file[0])
    {
        auto _res = get_resource(std::string(_name_file.begin() + 2, _name_file.end()));
        return load_value_from_resource<TYPE>(_res);
    }
    else
    {
        return load_value_from_file_on_disk<TYPE>(_name_file);
    }
}

template<typename TYPE>
TYPE MatrixUtils::load_value_from_file_on_disk(const std::string &_name_file)
{
    assertm(!_name_file.empty(), "filename is empty");
    std::ifstream file(_name_file, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + _name_file);
    }
    TYPE value;
    file.read(reinterpret_cast<char*>(&value), sizeof(TYPE));
    return value;
}

template<typename TYPE>
TYPE MatrixUtils::load_value_from_resource(const Resource& _resource)
{
    assertm(sizeof(TYPE) == _resource.size, "resource size with size of datatype must match");
    return *(reinterpret_cast<const TYPE*>(_resource.start));
}

template <typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixUtils::load_matrix(std::istream& file, int rows, int cols) {
    // Create a matrix with the specified dimensions
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> mat(rows, cols);
    
    // Calculate the expected data size in bytes
    size_t data_size = static_cast<size_t>(rows) * cols * sizeof(Scalar);
    
    // Read the matrix data directly into the matrix's memory
    if(!file.read(reinterpret_cast<char*>(mat.data()), data_size)) {
        throw std::runtime_error("Failed to read matrix data");
    }
    
    return std::move(mat);
}

template <typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixUtils::load_matrix_from_buffer(const void* _stream, int rows, int cols) {
    // Create a matrix with the specified dimensions
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> mat(rows, cols);
    
    // Calculate the expected data size in bytes
    size_t data_size = static_cast<size_t>(rows) * cols * sizeof(Scalar);
    
    // Read the matrix data directly into the matrix's memory
    memcpy(mat.data(), _stream, data_size);
    
    return std::move(mat);
}


// template<class TYPE>
// Eigen::Matrix<TYPE, Eigen::Dynamic, Eigen::Dynamic> MatrixUtils::load_from_file(const std::string& _name_file, const int& _rows, const int& _cols)
// {
//     Eigen::Matrix<TYPE, Eigen::Dynamic, Eigen::Dynamic> _ret;
//     std::ifstream _file(_name_file, std::ios_base::binary);
//     if(!_file.is_open())
//     {
//         throw std::runtime_error("opening file '" + _name_file + "' failed.");
//         return _ret;
//     }
//     _ret.conservativeResize(_rows, _cols);

//     _file.read(reinterpret_cast<char*>(_ret.data()), _ret.size() * sizeof(TYPE));
    
//     return std::move(_ret);
// }

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline void MatrixUtils::save_to_file(const std::string &_name_file, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> _matrix)
{
    assertm(!_name_file.empty(), "filename is empty");
    std::ofstream _file(_name_file, std::ios::binary);
    assertm(_file.is_open(), "could not open filename '" + _name_file + "' for writing.");
    int _rows = _matrix.rows();
    int _cols = _matrix.cols();
    int _type = 0;
    if (std::is_same<_Scalar, uint8_t>::value) {
        _type = 0;
    } else if (std::is_same<_Scalar, int8_t>::value) {
        _type = 1;
    } else if (std::is_same<_Scalar, uint16_t>::value) {
        _type = 2;
    } else if (std::is_same<_Scalar, int16_t>::value) {
        _type = 3;
    } else if (std::is_same<_Scalar, int32_t>::value) {
        _type = 4;
    } else if (std::is_same<_Scalar, float>::value) {
        _type = 5;
    } else if (std::is_same<_Scalar, double>::value) {
        _type = 6;
    } else if (std::is_same<_Scalar, std::complex<double>>::value) {
        _type = 14;
    } else {
        throw std::invalid_argument("Unsupported matrix type");
    }
    _file.write(reinterpret_cast<const char*>(&_rows), sizeof(_rows));
    _file.write(reinterpret_cast<const char*>(&_cols), sizeof(_cols));
    _file.write(reinterpret_cast<const char*>(&_type), sizeof(_type));
    _file.write(reinterpret_cast<const char*>(_matrix.data()), _matrix.size() * sizeof(_Scalar));
    _file.close();
    if(_file.fail()) {
        throw std::runtime_error("Error occurred while writing to file: " + _name_file);
    }
    // std::cerr << "matrix '" << _name_file << "':\n" << _matrix << "\n";
}

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline void MatrixUtils::release(Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& _matrix)
{
    _matrix.resize(0);
    _matrix.shrinkToFit();
}

template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void MatrixUtils::save_to_mat_file(const std::string& filename, const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& matrix, const std::string& varname)
{
    mat_t *matfp = Mat_CreateVer(filename.c_str(), nullptr, MAT_FT_MAT5);
    if (!matfp) {
        throw std::runtime_error("Could not create MAT file: " + filename);
    }

    size_t dims[2] = {static_cast<size_t>(matrix.rows()), static_cast<size_t>(matrix.cols())};
    enum matio_types data_type;
    enum matio_classes class_type;
    int isComplex = 0;

    if (std::is_same<_Scalar, double>::value) {
        data_type = MAT_T_DOUBLE;
        class_type = MAT_C_DOUBLE;
    } else if (std::is_same<_Scalar, float>::value) {
        data_type = MAT_T_SINGLE;
        class_type = MAT_C_SINGLE;
    } else if (std::is_same<_Scalar, int32_t>::value) {
        data_type = MAT_T_INT32;
        class_type = MAT_C_INT32;
    } else if (std::is_same<_Scalar, std::complex<double>>::value) {
        data_type = MAT_T_DOUBLE;
        class_type = MAT_C_DOUBLE;
        isComplex = 1;
    } else if (std::is_same<_Scalar, int16_t>::value) {
        data_type = MAT_T_INT16;
        class_type = MAT_C_INT16;
    } else if (std::is_same<_Scalar, uint16_t>::value) {
        data_type = MAT_T_UINT16;
        class_type = MAT_C_UINT16;
    } else if (std::is_same<_Scalar, int8_t>::value) {
        data_type = MAT_T_INT8;
        class_type = MAT_C_INT8;
    } else if (std::is_same<_Scalar, uint8_t>::value) {
        data_type = MAT_T_UINT8;
        class_type = MAT_C_UINT8;
    } else {
        Mat_Close(matfp);
        throw std::invalid_argument("Unsupported scalar type for MAT file");
    }

    matvar_t *matvar = nullptr;
    if (isComplex) {
        Eigen::Matrix<double, _Rows, _Cols> real_part = matrix.real();
        Eigen::Matrix<double, _Rows, _Cols> imag_part = matrix.imag();
        mat_complex_split_t complex_data;
        complex_data.Re = real_part.data();
        complex_data.Im = imag_part.data();
        matvar = Mat_VarCreate(varname.c_str(), class_type, data_type, 2, dims, &complex_data, MAT_F_COMPLEX);
    } else {
        matvar = Mat_VarCreate(varname.c_str(), class_type, data_type, 2, dims, const_cast<_Scalar*>(matrix.data()), 0);
    }

    if (!matvar) {
        Mat_Close(matfp);
        throw std::runtime_error("Could not create MAT variable");
    }

    if (Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB) != 0) {
        Mat_VarFree(matvar);
        Mat_Close(matfp);
        throw std::runtime_error("Could not write variable to MAT file");
    }

    Mat_VarFree(matvar);
    Mat_Close(matfp);
}

template<typename Scalar, int NumDims>
void MatrixUtils::save_to_mat_file(const std::string& filename, const Eigen::Tensor<Scalar, NumDims>& tensor, const std::string& varname)
{
    mat_t *matfp = Mat_CreateVer(filename.c_str(), nullptr, MAT_FT_MAT5);
    if (!matfp) {
        throw std::runtime_error("Could not create MAT file: " + filename);
    }

    size_t dims[NumDims];
    size_t total_elements = 1;
    for (int i = 0; i < NumDims; ++i) {
        dims[i] = static_cast<size_t>(tensor.dimension(i));
        total_elements *= dims[i];
    }

    enum matio_types data_type;
    enum matio_classes class_type;
    int isComplex = 0;

    if (std::is_same<Scalar, double>::value) {
        data_type = MAT_T_DOUBLE;
        class_type = MAT_C_DOUBLE;
    } else if (std::is_same<Scalar, float>::value) {
        data_type = MAT_T_SINGLE;
        class_type = MAT_C_SINGLE;
    } else if (std::is_same<Scalar, int32_t>::value) {
        data_type = MAT_T_INT32;
        class_type = MAT_C_INT32;
    } else if (std::is_same<Scalar, std::complex<double>>::value) {
        data_type = MAT_T_DOUBLE;
        class_type = MAT_C_DOUBLE;
        isComplex = 1;
        // Note: For tensors, complex handling is complex, may need to reshape or something
    } else if (std::is_same<Scalar, int16_t>::value) {
        data_type = MAT_T_INT16;
        class_type = MAT_C_INT16;
    } else if (std::is_same<Scalar, uint16_t>::value) {
        data_type = MAT_T_UINT16;
        class_type = MAT_C_UINT16;
    } else if (std::is_same<Scalar, int8_t>::value) {
        data_type = MAT_T_INT8;
        class_type = MAT_C_INT8;
    } else if (std::is_same<Scalar, uint8_t>::value) {
        data_type = MAT_T_UINT8;
        class_type = MAT_C_UINT8;
    } else {
        Mat_Close(matfp);
        throw std::invalid_argument("Unsupported scalar type for tensor MAT file");
    }

    matvar_t *matvar = nullptr;
    // For simplicity, assume non-complex for tensors; complex would need flattening
    if (isComplex) {
        Mat_Close(matfp);
        throw std::invalid_argument("Complex tensor support not implemented yet");
    } else {
        matvar = Mat_VarCreate(varname.c_str(), class_type, data_type, NumDims, dims, const_cast<Scalar*>(tensor.data()), 0);
        if (!matvar) {
            Mat_Close(matfp);
            throw std::runtime_error("Could not create MAT variable for tensor");
        }
    }

    if (Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE) != 0) {
        Mat_VarFree(matvar);
        Mat_Close(matfp);
        throw std::runtime_error("Could not write tensor variable to MAT file");
    }

    Mat_VarFree(matvar);
    Mat_Close(matfp);
}

#endif //__MATRIX_UTILS_H__
