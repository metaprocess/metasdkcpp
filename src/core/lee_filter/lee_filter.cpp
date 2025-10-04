#include "lee_filter.hpp"
#include "Utils.h"
#include "MatrixUtils.h"
#include <ElapsedTimer.h>
#include RESOURCES_HEADER

Eigen::MatrixXcd LeeFilter(const Eigen::MatrixXcd& input_image, int window_size, int ENL) {
    // Step 1: Handle zero values in the input image
    // Replace any complex zero (0 + 0i) with a small value (1e-20 + 0i)
    Eigen::MatrixXcd modified_input = input_image.unaryExpr([](const std::complex<double>& c) {
        return c == std::complex<double>(0.0, 0.0) ? std::complex<double>(1e-20, 0.0) : c;
    });

    // Step 2: Extract magnitude and phase
    // If the input is complex, the filter operates on the magnitude and preserves the phase
    Eigen::MatrixXd magnitude = modified_input.cwiseAbs();
    Eigen::MatrixXd phase = modified_input.unaryExpr([](const std::complex<double>& c) {
        return std::arg(c);
    });

    // Step 3: Compute noise variation coefficient
    double noise_variation_coef = 1.0 / static_cast<double>(ENL);

    // Step 4: Pad the magnitude matrix with zeros
    int pad_size = window_size / 2;
    int rows = magnitude.rows();
    int cols = magnitude.cols();
    Eigen::MatrixXd padded_magnitude = Eigen::MatrixXd::Zero(rows + 2 * pad_size, cols + 2 * pad_size);
    padded_magnitude.block(pad_size, pad_size, rows, cols) = magnitude;

    // Step 5: Initialize the filtered magnitude matrix
    Eigen::MatrixXd filtered_magnitude = Eigen::MatrixXd::Zero(rows, cols);

    // Step 6: Apply the Lee filter to each pixel
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // Extract the local window from the padded magnitude
            auto window = padded_magnitude.block(i, j, window_size, window_size);

            // Compute local mean
            double local_mean = window.mean();

            // Compute local variance: var = mean(x^2) - mean(x)^2
            // double local_variance = (window.array() - local_mean).square().sum() /  - static_cast<double>(window.size());
            double N = window.size();
            double local_variance = ((window.array() - local_mean).square()).sum() / (N - 1);

            // Compute local variation coefficient with a small epsilon to avoid division by zero
            double local_variation_coef = local_variance / (local_mean * local_mean/* + 1e-20*/);

            // Compute weighting function, clipped to [0, infinity) to prevent negative weights
            double weighting_function = 1.0 - noise_variation_coef / local_variation_coef;

            // Compute filtered pixel value
            double filtered_pixel = magnitude(i, j) * weighting_function + local_mean * (1 - weighting_function);

            // Store the filtered pixel
            filtered_magnitude(i, j) = filtered_pixel;
        }
    }

    // Step 7: Reconstruct the output image with original phase
    // Combine filtered magnitude with the original phase using element-wise operations
    Eigen::MatrixXcd output_image = filtered_magnitude.binaryExpr(phase, [](double mag, double p) {
        return mag * std::exp(std::complex<double>(0.0, p));
    });

    return std::move(output_image);
}

