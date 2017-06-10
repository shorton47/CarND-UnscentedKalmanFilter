//--------------------------------------------------------------------------------------------------------------
//
//
// CTRV Model !! 5 component state vector turning rate and velocity are constant. Augmented state is 7 compoentns

// Conventions:
// i = row, j = column
// _ means a Class level variable

// Commented out "angle normalization as I am not convinced it is needed. Double precision calc should be ok.
//  Lab is woring fine without it. I tried it both ways and not any better with all the checks.
//
//--------------------------------------------------------------------------------------------------------------


#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <iomanip>

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

//----------
// Constructor - Initializes Unscented Kalman filter
//----------
UKF::UKF() {

    //---
    // BOOLEAN program controls
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
    
    n_x_   = 5;  // State dimension
    n_aug_ = 7;  // Augmented state dimension
    
    num_sigma_pts_ = (2 * n_aug_) + 1;  // Number of sigma points (following convention)
    lambda_ = 3 - n_aug_;               // Sigma point spreading parameter
    
    // Class variables not currently in .h file... <TODO> CHECK ON THIS
    stepnum_ = 0;             // Keep track of each step in KF
    previous_timestamp_ = 0;  // Keep track of previous data points time step to calc delta time between measurement points
    

    // Initial State vector
    x_ = VectorXd(n_x_);
    //x_ << 0.0, 0.0, 0.0, 0.0, 0.0;
    x_.fill(0.0);
    
    x_aug_ = VectorXd(7);
    //x_aug_ << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    x_aug_.fill(0.0);
    
    // Initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);
    P_.Identity(n_x_, n_x_);
    //P_ << 1.0, 0.0, 0.0, 0.0, 0.0,
    //      0.0, 1.0, 0.0, 0.0, 0.0,
    //      0.0, 0.0, 1.0, 0.0, 0.0,
    //      0.0, 0.0, 0.0, 1.0, 0.0,
    //      0.0, 0.0, 0.0, 0.0, 1.0;
    
    // Augmented State covariance matrix
    P_aug_ = MatrixXd(7,7);
    P_aug_.Identity(n_aug_, n_aug_);
    //P_aug_ << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    //          0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    //          0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    //          0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    //          0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    //          0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    //          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
    
    ///<TODO> FIGURE THIS OUT* predicted sigma points matrix
    //MatrixXd Xsig_pred_;
    //Xsig_pred_.fill(0.0); //?????
    
    // Create sigma point matrix
    // Note: currently (7,15)
    //create sigma point matrix
    
    
    // Predicted Sigma Points
    Xsig_pred_ = MatrixXd(n_x_, num_sigma_pts_);
    Xsig_pred_.fill(0.0);

    // Generate Sigma Points ( Augmented)
    //MatrixXd Xsig_aug_ = MatrixXd(7, 15);
    Xsig_aug_ = MatrixXd(n_aug_, num_sigma_pts_);
    Xsig_aug_.fill(0.0);
    
    
    // Pre-compute fixed weights for sigma points
    // <TODO> Gonna need to move to .h!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    weights_ = VectorXd(num_sigma_pts_);
    weights_(0) = lambda_ /(lambda_ + n_aug_);
    for (int i=1; i<num_sigma_pts_; i++) {
        weights_(i) = 0.5 / (n_aug_ + lambda_);
    }
    

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 2.0; // W=0.5 1.0 0.9 2.0 was 30

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.875; // W=.875 W=.25 0.50 0.4 3.0 was 30

    //---
    // Measurement device noise(s) provided by manufacturer
    //---
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
    
    
    //set measurement dimension, radar can measure rho, phi, and rho_dot
    n_z_ = 3;
    
    // Measurement noise covariance matrix for RADAR
    R_ = MatrixXd(n_z_,n_z_);
    R_ << std_radr_*std_radr_,                       0,                     0,
                            0, std_radphi_*std_radphi_,                     0,
                            0,                       0, std_radrd_*std_radrd_;
    
    // Process Noise - 2 elements, stda,
    n_Q_ = 2;
    Q_ = MatrixXd(n_Q_,n_Q_);
    Q_ << std_a_*std_a_,                     0,
                      0, std_yawdd_*std_yawdd_;

    // Measurement noise covariance matrix
    R_laser_ = MatrixXd(2,2);
    R_laser_ << std_laspx_*std_laspx_,                     0,
                                    0, std_laspy_*std_laspy_;
    

}


//----------
// Destructor
//----------
UKF::~UKF() {}


//---
// ProcessMeasurement Method - called for every point from the simulator
//---
/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  
    //---
    // Step 0. - Initialization of Measurements (For 1st incoming measurement point only)
    //---
    if (!is_initialized_) {
        
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
            
            // Debug
            //if (DEBUG_) cout << "----------\n";
            cout << "----- UKF: Step #" << stepnum_ << "-----\n";
            cout << "Init measurement is RADAR: rho,phi,rhodot=: " << rho << " " << phi << " " << rhodot << endl;
            cout << "Init measurement is RADAR: px,py: " << px << " "<< py << endl;

        } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            
            // Get LASER measurements from measurement_pack into state variable ekf_.x
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
    // Step 1. - Kalman Filter Prediction & Update Step based on Boolean Control
    //---
    stepnum_ += 1;
    
    // RADAR
    if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        
        // Update time step
        double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;  // dt converted to sec from microseconds
    
        cout << "---------- UKF Predict w/ RADAR: Step #" << stepnum_ << " Delta time=" << dt << "\n";
        double rho    = meas_package.raw_measurements_(0);
        double phi    = meas_package.raw_measurements_(1);
        double rhodot = meas_package.raw_measurements_(2);
        cout << "New incoming point: rho,phi,rhodot= " << rho << " " << phi << " " << rhodot << endl;

        // Kalman Predicition & Update
        Prediction(dt);
        UpdateRadar(meas_package);
        
        previous_timestamp_ = meas_package.timestamp_;  // Save for next dt calculation
        
    // LASER
    } else if (use_laser_ && (meas_package.sensor_type_ == MeasurementPackage::LASER)) {
        
        // Update time step
        double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;  // dt converted to sec from microseconds
        
        cout << "---------- UKF Predict w/ LASER: Step #" << stepnum_ << " Delta time=" << dt << "\n";
        double px = meas_package.raw_measurements_(0);
        double py = meas_package.raw_measurements_(1);
        cout << "New incoming point: px,py= " << px << " " << py <<  endl;
        
        // Kalman Predicition & Update
        Prediction(dt);
        UpdateLaser(meas_package);
        
        previous_timestamp_ = meas_package.timestamp_;  // Save for next dt calculation

        
    } else {
        if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            cout << "---------- UKF: Step #" << stepnum_ << " Skipping LASER point" << "-----\n";
        } else {
            cout << "---------- UKF: Step #" << stepnum_ << " Skipping RADAR point" << "-----\n";
        }
    }
    
    return;
} // ProcessMeasurement


//---------------
// Predicition Step of Unscented Kalman Filter (UKF). Predicts the State Mean and State Covariance Matrix using Sigma Points.
// * @param {double} delta_t the change in time (in seconds) between the last measurement and this one.
//---------------
void UKF::Prediction(double delta_t) {
  
    VectorXd x_diff(num_sigma_pts_);
    
    
    //----------
    // #1. Generate Sigma Points (augmented with process noise)
    //----------
   
    // Create Augmented State Mean
    x_aug_.head(n_x_) = x_;  // x_ is the previous state (before the new measurement update that has just come in)
    x_aug_(n_aug_-2)  = 0;   // Assume acceleration process noise is zero
    x_aug_(n_aug_-1)  = 0;   // Assume angle acceleration process noise is zero
    
    //Create Augmented Covariance Matrix
    P_aug_.fill(0.0);
    P_aug_.topLeftCorner(n_x_, n_x_) = P_;  // State covariance matrix P is top-left
    P_aug_.bottomRightCorner(2,2) = Q_;     // Process noise covariance matrix Q is bottom-right
    
    // Square root of augmented covariance matrix
    MatrixXd L = P_aug_.llt().matrixL();
    
    // Create Sigma Point Matrix (Augmented)
    Xsig_aug_.col(0)  = x_aug_;  // 1st col of sigma points in matrix is State vector
    for (int i = 0; i < n_aug_; i++) {
        Xsig_aug_.col(i+1)        = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);  // Col vector add
        Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
    }
    
    PrintEigenColVector("Pred: x_aug=",&x_aug_);
    PrintEigenMatrix("P_aug=",&P_aug_);
    PrintEigenMatrix("L=",&L);
    PrintEigenMatrix("Xsig_aug=",&Xsig_aug_);
    
    
    //----------
    // #2. Predict Sigma Points
    //----------
    double delta_t2 = delta_t*delta_t;
    
    // For each column vector of Xsig_pred from Xsig_aug
    for (int j=0; j<num_sigma_pts_; j++) {
        
        double px     = Xsig_aug_(0,j);
        double py     = Xsig_aug_(1,j);
        double v      = Xsig_aug_(2,j);
        double phi    = Xsig_aug_(3,j);
        double phidot = Xsig_aug_(4,j);
        double nua    = Xsig_aug_(5,j);
        double nuphia = Xsig_aug_(6,j);
        
        // Check if phidot is near zero to avoid div by zero
        if (fabs(phidot) > __DBL_EPSILON__) {
            
            // Process equations w/ sigma points
            Xsig_pred_(0,j) = px     + (v/phidot) * (sin(phi+(phidot*delta_t)) - sin(phi)) \
                                     + 0.5*delta_t2*cos(phi)*nua;
            Xsig_pred_(1,j) = py     + (v/phidot) * (-cos(phi+(phidot*delta_t)) + cos(phi)) \
                                     + 0.5*delta_t2*sin(phi)*nua;
            Xsig_pred_(2,j) = v      + 0.0 + delta_t*nua;
            Xsig_pred_(3,j) = phi    + phidot*delta_t + 0.5*delta_t2*nuphia;
            Xsig_pred_(4,j) = phidot + 0.0 + delta_t*nuphia;
        } else {
            Xsig_pred_(0,j) = px     + v*cos(phi)*delta_t + 0.5*delta_t2*cos(phi)*nua;
            Xsig_pred_(1,j) = py     + v*sin(phi)*delta_t + 0.5*delta_t2*sin(phi)*nua;
            Xsig_pred_(2,j) = v      + 0.0 + delta_t*nua;
            Xsig_pred_(3,j) = phi    + 0.0 + 0.5*delta_t2*nuphia;
            Xsig_pred_(4,j) = phidot + 0.0 + delta_t*nuphia;
        }
    }
    PrintEigenMatrix("Xsig_pred_=",&Xsig_pred_);
    
    
    //----------
    // #3. Predict State Mean Vector (x_) & State Covariance Matrix (P_) w/ Unscented Kalman Filter Equations
    //----------
    
    // Calc Predicted State Mean x_ (weighted)
    x_.fill(0.0);  // Note: current x_ has already been captured in sigma points, so zeroing out for new sum
    for (int j=0; j<num_sigma_pts_; j++) {
        x_ +=  weights_(j) * Xsig_pred_.col(j);  // Col vector mult & add
    }
   
    // Calc Predicted State Covariance Matrix P_ (weighted)
    P_.fill(0.0);  // Note: current P_ has already been captured in sigma points, so zeroing out for new sum
    for (int j = 0; j <num_sigma_pts_; j++) {
        
        x_diff = Xsig_pred_.col(j) - x_ ;  // State Vector difference (Pred - Current)

        // Angle normalization for phidot (? not needed)
        //while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        //while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        P_ += weights_(j) * x_diff * x_diff.transpose();
    }

    // Debug (prediction of x_ & P_ before measurement update)
    PrintEigenColVector("x_=",&x_);
    PrintEigenMatrix("P_=",&P_);
    
    return;
} // Prediction


//---------------
// Update Step - Radar. Updates the State Mean x_ and the State Covariance Matrix P_.
//@param {MeasurementPackage} meas_package
//---------------
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  
    MatrixXd Zsig = MatrixXd(n_z_, num_sigma_pts_);  // Create matrix for sigma points in measurement space
    MatrixXd S    = MatrixXd(n_z_, n_z_);            // Measurement covariance matrix S
    MatrixXd Tc   = MatrixXd(n_x_, n_z_);            // Measurement difference cross correlation matrix
    MatrixXd K    = MatrixXd(n_z_, n_z_);            // Kalman Gain matrix
    
    VectorXd z       = VectorXd(n_z_);  // New mean measurements
    VectorXd z_pred  = VectorXd(n_z_);  // Mean predicted measurements
    VectorXd delta_z = VectorXd(n_z_);  // Delta means


    //---
    // #1. Predict Measurement using predicted measurement sigma points
    //---
    
    for (int j=0; j<num_sigma_pts_; j++) {
        
        // Transform predicted sigma points (from prediction step) into measurement space
        double px     = Xsig_pred_(0,j);
        double py     = Xsig_pred_(1,j);
        double v      = Xsig_pred_(2,j);
        double yaw    = Xsig_pred_(3,j);
        //double yawdot = Xsig_pred_(4,j);
        
        double v1   = v*cos(yaw);
        double v2   = v*sin(yaw);
        double term = sqrt(px*px + py*py);
        
        // Measurement model
        Zsig(0,j) = term;                        // rho
        Zsig(1,j) = atan2(py,px);                // phi NOTE: ATAN2 handles zero case & normalizes (-pi,pi)
        if (term > __DBL_EPSILON__) {
            Zsig(2,j) = (px*v1 + py*v2) / term;  //rhodot
        } else {
            Zsig(2,j) = 0.0;
        }
    }
   
    
    //---
    // #2. Calculate mean predicted measurement
    //---
    z_pred.fill(0.0);
    for (int j=0; j<num_sigma_pts_; j++) {
        z_pred += weights_(j) * Zsig.col(j);
    }
    
    
    //---
    // #3. Calculate measurement covariance matrix S
    //---
    S.fill(0.0);
    for (int j=0; j<num_sigma_pts_; j++) {
        
        //residual
        VectorXd z_diff = VectorXd(n_z_);
        z_diff = Zsig.col(j) - z_pred;
        
        // Angle normalization for phi (? not needed)
        //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
            
        S += (weights_(j) * z_diff * z_diff.transpose()) ;
    }
    
    // Add "measurement noise" to covariance matrix once at end (noise is invarient in this case)
    S = S + R_;
    PrintEigenMatrix("S=",&S);
    
    
    //---
    // #4. Update State Mean and Covariance (RADAR) using Unscented Kalman Filter (UKF) Equations
    //---
    
    // Calculate cross correlation matrix Tc
    Tc.fill(0.0);
    for (int j = 0; j < num_sigma_pts_; j++) {
        
        // Residual
        VectorXd z_diff = Zsig.col(j) - z_pred;
        
        //angle normalization check
        //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        // State vector difference
        VectorXd x_diff = Xsig_pred_.col(j) - x_;
        
        //angle normalization check
        //while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        //while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        Tc += weights_(j) * x_diff * z_diff.transpose();
    }
    
    // Calculate Kalman gain K
    K = Tc * S.inverse();
    
    // New incoming measurement & delta
    
    z = meas_package.raw_measurements_;  // In sensor space
    delta_z = z - z_pred;                // residual (New Measurement - Predicited)
    
    //angle normalization
    //while (delta_z(1)> M_PI) delta_z(1)-=2.*M_PI;
    //while (delta_z(1)<-M_PI) delta_z(1)+=2.*M_PI;
    
    // Update State Mean & Covariance matrix
    x_ = x_ + K * delta_z;
    P_ = P_ - K*S*K.transpose();
    
    // Debug (Update of x_ & P_)
    PrintEigenMatrix("K=",&K);
    PrintEigenColVector("Update: x_=",&x_);
    PrintEigenMatrix("P_=",&P_);
    
    
    //---
    // Step Prologue. - Calculate Radar NIS
    //---
    double radar_NIS = delta_z.transpose() * S.inverse() * delta_z;
    cout << "Radar NIS=" << radar_NIS << endl;
    
    return;
}  // UpdateRadar


//---------------
// Update Step - LASER. Updates the State Mean x_ and the State Covariance Matrix P_ w/ incoming measurement.
// @param {MeasurementPackage} meas_package
//---------------
void UKF::UpdateLaser(MeasurementPackage meas_package) {
    

    int n_z_laser = 2;
    MatrixXd Zsig = MatrixXd(n_z_laser, num_sigma_pts_);  // Create matrix for sigma points in measurement space
    MatrixXd S    = MatrixXd(n_z_laser, n_z_laser);       // Measurement covariance matrix S
    MatrixXd Tc   = MatrixXd(n_x_, n_z_laser);            // Measurement difference cross correlation matrix
    MatrixXd K    = MatrixXd(n_z_, n_z_laser);            // Kalman Gain matrix
    
    VectorXd z       = VectorXd(n_z_laser);  // New mean measurements
    VectorXd z_pred  = VectorXd(n_z_laser);  // Mean predicted measurements
    VectorXd delta_z = VectorXd(n_z_laser);  // Delta means
    
    
    //---
    // #1. Predict Measurement using predicted measurement sigma points
    //---
    
    for (int j=0; j<num_sigma_pts_; j++) {
        
        // Transform predicted sigma points (from prediction step) into measurement space
        double px = Xsig_pred_(0,j);
        double py = Xsig_pred_(1,j);
        
        // Measurement model (straight through)
        Zsig(0,j) = px;  // p sub x
        Zsig(1,j) = py;  // p sub y
    }
    //PrintEigenMatrix("UpdateLaser: Zsig=",&Zsig);

    
    //---
    // #2. Calculate mean predicted measurement
    //---
    z_pred.fill(0.0);
    for (int j=0; j<num_sigma_pts_; j++) {
        z_pred += weights_(j) * Zsig.col(j);
    }
    //PrintEigenColVector("zpred=",&z_pred);
    
    //---
    // #3. Calculate measurement covariance matrix S
    //---
    S.fill(0.0);
    for (int j=0; j<num_sigma_pts_; j++) {
        
        //residual
        VectorXd z_diff = VectorXd(n_z_laser);
        z_diff = Zsig.col(j) - z_pred;
        
        // Angle normalization for phi (? not needed)
        //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        S += weights_(j) * z_diff * z_diff.transpose();
    }
    
    // Add "measurement noise" to covariance matrix once at end (noise is invarient in this case)
    S = S + R_laser_;
    PrintEigenMatrix("S=",&S);
    
    
    //---
    // #4. Update State Mean and Covariance (RADAR) using Unscented Kalman Filter (UKF) Equations
    //---
    
    // Calculate cross correlation matrix Tc
    Tc.fill(0.0);
    for (int j = 0; j < num_sigma_pts_; j++) {
        
        // Residual
        VectorXd z_diff = Zsig.col(j) - z_pred;
        
        //angle normalization check
        //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        // State vector difference
        VectorXd x_diff = Xsig_pred_.col(j) - x_;
        
        //angle normalization check
        //while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        //while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        Tc += weights_(j) * x_diff * z_diff.transpose();
    }
    
    // Calculate Kalman gain K
    K = Tc * S.inverse();
    
    // New incoming measurement & delta
    
    z = meas_package.raw_measurements_;  // In sensor space
    delta_z = z - z_pred;                // residual (New Measurement - Predicited)
    
    //angle normalization
    //while (delta_z(1)> M_PI) delta_z(1)-=2.*M_PI;
    //while (delta_z(1)<-M_PI) delta_z(1)+=2.*M_PI;
    
    // Update State Mean & Covariance matrix
    x_ = x_ + K * delta_z;
    P_ = P_ - K*S*K.transpose();
    
    // Debug (Update of x_ & P_)
    PrintEigenMatrix("K=",&K);
    PrintEigenColVector("Update: x_=",&x_);
    PrintEigenMatrix("P_=",&P_);
    
    
    //---
    // Step Prologue. - Calculate Radar NIS
    //---
    double laser_NIS = delta_z.transpose() * S.inverse() * delta_z;
    cout << "Laser (Lidar) NIS=" << laser_NIS << endl;
    
    return;
}  // UpdateLaser


//----------
// Utility Methods for UKF
//----------
void UKF::PrintEigenMatrix(string label, MatrixXd *mat) {

    cout << label << endl;
    for (int i=0; i< mat->rows(); i++) {
        for (int j=0; j< mat->cols(); j++) {
            cout << setw(12) << setprecision(4) << (*mat)(i,j);
        }
        cout << endl;
    }
    return;
}

//
void UKF::PrintEigenColVector(string label, VectorXd *vector) {
    
    cout << label << endl;
    for (int i=0; i< vector->rows(); i++) {
            cout << setw(14) << setprecision(8) << (*vector)(i);
        cout << endl;
    }
    return;
}



