#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


//
// <TODO> Probably need to move all this to .h file!!!!!!
//




//---
// Constructor - Initializes Unscented Kalman filter
//---
UKF::UKF() {

    //---
    // BOOLEAN controls
    //---
    // Initially set to false, set to true in first call of ProcessMeasurement
    is_initialized_ = false;

    // If this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // If this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    //---
    // State Space
    //--
    
    n_x_ = 5;    // State dimension
    n_aug_ = 7;  // Augmented state dimension
    num_sigma_pts_ = 2 * n_aug_ + 1;   // Number of sigma points (following convention)
    lambda_ = 3 - n_aug_;   // Sigma point spreading parameter
    
    // Class variables not currently in .h file... <TODO> CHECK ON THIS
    stepnum_ = 0;             // Keep track of each step in KF
    previous_timestamp_ = 0;  // Keep track of previous time to calc delta time between measurement points
    
    // Number of sigma points (following convention)
    //int num_sigma_pts;
    
    // Sigma point spreading parameter
    //double lambda;
    
    // Epsilon for near zero tests
    //const double EPS_ = 1.0e-15;
    
   
    
    // More time
    //double dt,dt2;
    
    
    // Augmented State vector
    VectorXd x_aug = VectorXd(7);
    MatrixXd P_aug = MatrixXd(7,7);
    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(7, 15);
    
    
    
    
    
    
    // Initial State vector
    x_ = VectorXd(n_x_);
    x_ << 0.0, 0.0, 0.0, 0.0, 0.0;
    x_.fill(0.0);
    
    x_aug_ << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    x_aug_.fill(0.0);
    
    // Initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);
    P_ << 1.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 1.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 1.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 1.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 1.0;
    
    // Augmented State covariance matrix
    P_aug_ << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
             0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
    
    ///<TODO> FIGURE THIS OUT* predicted sigma points matrix
    //MatrixXd Xsig_pred_;
    Xsig_pred_.fill(0.0);

    // Create sigma point matrix
    // Note: currently (7,15)
    Xsig_aug_ = MatrixXd(n_aug_, num_sigma_pts_);

    // Weights for sigma points
    VectorXd weights = VectorXd(num_sigma_pts_);
    weights(0) = lambda_ /(lambda_ + n_aug_);
    for (int i=1; i<num_sigma_pts_; i++) {
        weights(i) = 0.5 / (n_aug_ + lambda_);
    }
    
    ///* time when the state is true, in us <TODO> (??? - from .h file)
    time_us_ = 0L;

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 30;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 30;

    // Laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15;

    // Laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;

    // Radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;

    // Radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.03;

    // Radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.3;
    
    

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}


//---
// Destructor
//---
UKF::~UKF() {}


//---
// ProcessMeasurement Method - called for every point from the simulator
//---
/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  
    /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
    //---
    // Step 0. - Initialization of Measurements (For 1st incoming measurement point only)
    //---
    if (!is_initialized_) {
        //x_ = VectorXd(4);
        //x_ << 1, 1, 1, 1;  // * Initialize the state x_ with the first measurement.
        
        // Get initial RADAR measurements from measurement_pack into state variable ekf_.x_
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            int vec_size = meas_package.raw_measurements_.size();
            if (vec_size != 3) {
                cout << "ProcessMeasurement: Error - 1st RADAR data point is incorrect size=" << vec_size;
                return;
            }
            
            // Set the state with the initial location and zero velocity
            double rho    = meas_package.raw_measurements_(0);
            double phi    = meas_package.raw_measurements_(1);
            double rhodot = meas_package.raw_measurements_(2);
            
            double px = rho*cos(phi);  // Convert from polar to cartesian
            double py = rho*sin(phi);
            
            // Initialize
            x_ << px , py, 0.0, 0.0, 0.0;                // Can only derive position (not velocity) from rho & phi
            //<TODO> FIGURE THIS OUT
            //Hj_ = tools.CalculateJacobian(ekf_.x_);  // Initialize Jacobian Hj_ w/ 1st point
            
            // Debug
            cout << "----------\n";
            cout << "UKF: Step #" << stepnum_ << "\n";
            cout << "Init measurement is RADAR: rho,phi,rhodot=: " << rho << " " << phi << " " << rhodot << endl;
            cout << "Init measurement is RADAR: px,py: " << px << " "<< py << endl;
        }
        // Get LASER measurements from measurement_pack into state variable ekf_.x_
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            int vec_size = meas_package.raw_measurements_.size();
            if (vec_size != 2) {
                cout << "ProcessMeasurement: Error - 1st LASER data point incorrect size=" << vec_size;
                return;
            }
            
            // Set the state with the initial location and zero velocity
            double px = meas_package.raw_measurements_(0);
            double py = meas_package.raw_measurements_(1);
            
            // Initialize
            x_ << px , py, 0.0, 0.0, 0.0;  // Only have position (not velocity) at this point
            
            // Debug
            cout << "----------\n";
            cout << "UKF: Step #" << stepnum_ << "\n";
            cout << "Init measurement is LASER (LIDAR): px,py: " << px << " " << py << " " << endl;
        }
        
        // Save initial time stamp for dt calculations
        previous_timestamp_ = meas_package.timestamp_;
        cout << "1st timestamp=" << previous_timestamp_ << endl;
        
        // Done initializing, no need to predict or update
        is_initialized_ = true;
        
        return;
    } // if !is_initialized
    
    
    //---
    // Major Step 1. - Kalman Filter Prediction Step
    //---
    stepnum_ += 1;
    
    // Update time step
    double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;  // dt converted to sec from microseconds
    previous_timestamp_ = meas_package.timestamp_;                            // Save for next dt calculation
    
    // Debug ONLY
    cout << "----------\n";
    cout << "UKF: Step #" << stepnum_ << " Delta time=" << dt << "\n";
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        cout << "New incoming point is LASER" << endl;
        double px = meas_package.raw_measurements_(0);
        double py = meas_package.raw_measurements_(1);
        cout << "New incoming point: px,py= " << px << " " << py <<  endl;
    } else {
        cout << "New incoming point is RADAR" << endl;
        double rho    = meas_package.raw_measurements_(0);
        double phi    = meas_package.raw_measurements_(1);
        double rhodot = meas_package.raw_measurements_(2);
        cout << "New incoming point: rho,phi,rhodot= " << rho << " " << phi << " " << rhodot << endl;
        
    }
    cout << "Start x_ = \n" << x_ << endl;
    cout << "Start P_ = \n" << P_ << endl;
    
    
    Prediction(dt);
    
    
    //---
    // Major Step 2. - Kalman Filter Update (Measurement) Step
    //---
    
    
    
    
    
    return;
} // ProcessMeasurement


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    
    //---
    // Step 1. - Generate Sigma Points (Augmented with process noise)
    //---
    
    //create augmented State mean
    x_aug_.head(n_x_) = x_; // where n is the number of elements from first element, and y is an input vector of that size.
    x_aug_(n_aug_-2) = 0;  // Assume Process noise is zero <TODO> CHECK THIS
    x_aug_(n_aug_-1) = 0;
    std::cout << "x_aug=\n" << x_aug_ << "\n";
    
    //create augmented covariance matrix
    P_aug_.fill(0.0);                               // Initialize full matrix to zero
    P_aug_.topLeftCorner(n_x_, n_x_) = P_;             // Standard P covariance matx in top-left
    P_aug_(n_aug_-2,n_aug_-2) = std_a_*std_a_;          // Augment bot-right Q matrix w/ variances
    P_aug_(n_aug_-1,n_aug_-1) = std_yawdd_*std_yawdd_;  // Augment bot-right Q matrix w/ variances
    std::cout << "P_aug=\n" << P_aug_ << "\n";
    
    //create square root of augmented covariance (P) matrix
    MatrixXd L = P_aug_.llt().matrixL();
    std::cout << "L=\n" << L << "\n";
    
    //---
    // Create augmented sigma points
    //---
    // 1st column of sigma point matrix
    Xsig_aug_.col(0)  = x_aug_;
    
    // Set remaining sigma points
    for (int i = 0; i < n_aug_; i++)
    {
        Xsig_aug_.col(i+1)        = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);
        Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
    }
    cout << "Xsig_aug=\n" << Xsig_aug_ << "\n";
    // <TODO> - PROBABLY SHOULD MOVE INTO DIFFERENT SEGMENT AS THIS IS INVARIENT WITH TIME

    //---
    // Step 2. - Predict Sigma Points
    //---
    //double delta_t = 0.1; //time diff in sec
    //create matrix with predicted sigma points as columns
    
    
    //Put in constructor
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    

    
    // Xsig_aug is x at K, Xsig_pred is x at k+1
    double delta_t2 = delta_t*delta_t;
    
    // For each column vector sigma point j of Xsig_aug
    for (int j=0; j<num_sigma_pts_; j++) {
        
        double px     = Xsig_aug_(0,j);
        double py     = Xsig_aug_(1,j);
        double v      = Xsig_aug_(2,j);
        double phi    = Xsig_aug_(3,j);
        double phidot = Xsig_aug_(4,j);
        double nua    = Xsig_aug_(5,j);
        double nuphia = Xsig_aug_(6,j);
        
        // Check is phidot is near zero to avoid div by zero
        if (fabs(phidot) > __DBL_EPSILON__) {
            
            // Process equations w/ sigma points
            Xsig_pred_(0,j) =  px    + (v/phidot) * (sin(phi+(phidot*delta_t)) - sin(phi)) \
                                    + 0.5*delta_t2*cos(phi)*nua;
            Xsig_pred_(1,j) =  py    + (v/phidot) * (-cos(phi+(phidot*delta_t)) + cos(phi)) \
                                    + 0.5*delta_t2*sin(phi)*nua;
            Xsig_pred_(2,j) = v      + 0 + delta_t*nua;
            Xsig_pred_(3,j) = phi    + phidot*delta_t + 0.5*delta_t2*nuphia;
            Xsig_pred_(4,j) = phidot + 0 + delta_t*nuphia;
            
        } else {
            
            // Process equations w/ sigma points when phidot zero
            Xsig_pred_(0,j) = px     + v*cos(phi)*delta_t + 0.5*delta_t2*cos(phi)*nua;
            Xsig_pred_(1,j) = py     + v*sin(phi)*delta_t + 0.5*delta_t2*sin(phi)*nua;
            Xsig_pred_(2,j) = v      + 0 + delta_t*nua;
            Xsig_pred_(3,j) = phi    + 0 + 0.5*delta_t2*nuphia;
            Xsig_pred_(4,j) = phidot + 0 + delta_t*nuphia;
            
        }
        
    }
    cout << "Xsig_pred=\n" << Xsig_pred_ << "\n";
    
    
    
    
    
    
    //---
    // Step 3. - Predict State Mean & State Covariance matrix of Kalman Filter Equations
    //---
    
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
