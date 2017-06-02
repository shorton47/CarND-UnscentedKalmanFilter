#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
    VectorXd rmse(4);
    VectorXd residual;
    rmse << 0,0,0,0;
    
    // Check the validity of inputs: estimation vector size should not be zero
    if (estimations.size() == 0 || (estimations.size() != ground_truth.size())) {
        cout << "Tools::CalculateRMSE: Error - estimation & ground truth vectors empty or not same size.";
        return rmse;
    }
    
    // Sum squared residuals
    for (int i=0; i < estimations.size(); ++i){
        residual = estimations[i] - ground_truth[i];
        residual = residual.array()*residual.array();  // Coefficient-wise multiplication
        rmse += residual;
    }
    
    // Calculate square root of means
    rmse = rmse/estimations.size();
    rmse = rmse.array().sqrt();
    
    return rmse;
    
    
    
    
    
    
    
    
    
    
    return rmse;
}
