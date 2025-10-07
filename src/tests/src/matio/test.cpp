// tests/test_example.cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "matio.h"
#include "MatrixUtils.h"
#include "CommonConstants.h"
#include <ElapsedTimer.h>

MatrixCompareResult test_matio()
{
    auto _input_image = std::get<Eigen::MatrixXd>(MatrixUtils::load_mat(":/src/matio/resources/test_linspace.mat", "a"));

    int _rows = 10000;
    int _cols = 1;
    auto _ref = Eigen::VectorXd::LinSpaced(_rows, 1., 1000);
    auto _ref_mat = Eigen::MatrixXd(_rows, _cols);
    _ref_mat.col(0) = _ref;



    // int _window_size = MatrixUtils::load_value_from_file<double>(":/src/lee_filter/resources/lee_filter_window_size.bin");
    // int _ENL = MatrixUtils::load_value_from_file<double>(":/src/lee_filter/resources/lee_filter_ENL.bin");
    // ElapsedTimer _timer;
    // auto _ret_filter = LeeFilter(_input_image, _window_size, _ENL);
    // _timer.print_and_restart("lee_filter");
    // // MatrixUtils::save_to_file("c_lee_filter_output_image.bin"}), _ret_filter)
    // auto _matrix_ref_output = std::get<Eigen::MatrixXcd>(MatrixUtils::load_from_file(":/src/lee_filter/resources/lee_filter_output_image.bin"));

    return MatrixUtils::compare_matrices(_input_image, _ref_mat);
    // std::cout << _ret << std::endl;
    
    // bool _test_result = (_ret.match() || (_ret.sum_diff() <= CONST_ACCEPTABLE_SUM_DIFF_ERROR && _ret.max_diff() <= CONST_ACCEPTABLE_MAX_DIFF_ERROR));
    // return _test_result;
}

// Traditional test style
TEST_CASE("TEST_MATIO", "[matrix]") {
    // REQUIRE(true == test_lee_filter());
    // auto _result = MatrixCompareResult(false, 1, 2);
    
    SECTION("DIFF WITH REF") {
        auto _result = test_matio();
        REQUIRE_THAT(_result.max_diff(), Catch::Matchers::WithinAbs(0, CONST_ACCEPTABLE_MAX_DIFF_ERROR));
        REQUIRE_THAT(_result.sum_diff(), Catch::Matchers::WithinAbs(0, CONST_ACCEPTABLE_SUM_DIFF_ERROR));
    }

    // SECTION("SUM DIFF") {
    //     REQUIRE(_result.sum_diff() <= CONST_ACCEPTABLE_SUM_DIFF_ERROR);
    // }
}

// // BDD style
// SCENARIO("String manipulation", "[strings]") {
//     GIVEN("An empty string") {
//         std::string s;
        
//         WHEN("Appending text") {
//             s += "hello";
            
//             THEN("Size increases") {
//                 REQUIRE(s.size() == 5);
//         }
//     }
// }
