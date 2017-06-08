#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <iomanip>

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

//---
// Constructor - Initializes Unscented Kalman filter
//---
UKF::UKF() {

    
    //std::setw(9);
    //std::setprecision(4);
    
    //---
    // BOOLEAN program controls
    //---
    // Initially set to false, set to true in first call of ProcessMeasurement
    is_initialized_ = false;

    // If this is false, laser measurements will be ignored (except during init)
    use_laser_ = false;

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
    


    // More time
    //double dt,dt2;
    
    
    // Augmented State vector
   
  
    
    
    
    
    
    
    // Initial State vector
    x_ = VectorXd(n_x_);
    //x_ << 0.0, 0.0, 0.0, 0.0, 0.0;
    x_.fill(0.0);
    
    x_aug_ = VectorXd(7);
    //x_aug_ << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    x_aug_.fill(0.0);
    
    // Initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);
    P_ << 1.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 1.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 1.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 1.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 1.0;
    
    // Augmented State covariance matrix
    P_aug_ = MatrixXd(7,7);
    P_aug_ << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
    
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
    std_a_ = 0.50; // 1.0 0.9 2.0 was 30

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = .25; // 0.50 0.4 3.0 was 30

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
    
    // Measurement noise covariance matrix
    R_ = MatrixXd(n_z_,n_z_);
    R_ << std_radr_*std_radr_,                       0,                     0,
                            0, std_radphi_*std_radphi_,                     0,
                            0,                       0, std_radrd_*std_radrd_;
    
    

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
            //x_aug_ << px, py, 0.0, 0.0, 0.0, 0.0, 0.0;
            //<TODO> FIGURE THIS OUT
            //Hj_ = tools.CalculateJacobian(ekf_.x_);  // Initialize Jacobian Hj_ w/ 1st point
            
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
            //x_aug_ << px, py, 0.0, 0.0, 0.0, 0.0, 0.0;
            
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
    //double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;  // dt converted to sec from microseconds
    //previous_timestamp_ = meas_package.timestamp_;                            // Save for next dt calculation
    
    // Debug ONLY
    /*
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
    cout << "Start xaug = \n" << x_aug_ << endl;
    cout << "Start Paug = \n" << P_aug_ << endl;
    */
    
    //Prediction(dt);
    
    
    //---
    // Major Step 2. - Kalman Filter Update (Measurement) Step
    //---

    if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        
        // Update time step
        double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;  // dt converted to sec from microseconds
        //previous_timestamp_ = meas_package.timestamp_;

        cout << "----------\n";
        cout << "UKF Predict w/ Radar: Step #" << stepnum_ << " Delta time=" << dt << "\n";
        double rho    = meas_package.raw_measurements_(0);
        double phi    = meas_package.raw_measurements_(1);
        double rhodot = meas_package.raw_measurements_(2);
        cout << "New incoming point: rho,phi,rhodot= " << rho << " " << phi << " " << rhodot << endl;

        Prediction(dt);
        UpdateRadar(meas_package);
        
        previous_timestamp_ = meas_package.timestamp_;   // Save for next dt calculation
        
        
    } else if (use_laser_ && (meas_package.sensor_type_ == MeasurementPackage::LASER)) {
        
        // Update time step
        double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;  // dt converted to sec from microseconds
        
        cout << "----------\n";
        cout << "UKF Predict & LASER: Step #" << stepnum_ << " Delta time=" << dt << "\n";
        double px = meas_package.raw_measurements_(0);
        double py = meas_package.raw_measurements_(1);
        cout << "New incoming point: px,py= " << px << " " << py <<  endl;
        
        Prediction(dt);
        UpdateLidar(meas_package);
        
        previous_timestamp_ = meas_package.timestamp_;                            // Save for next dt calculation

        
    } else {
        cout << "----------\n";
        if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            cout << "----- UKF: Step #" << stepnum_ << " Skipping LASER point" << "-----\n";
        } else {
            cout << "----- UKF: Step #" << stepnum_ << " Skipping RADAR point" << "-----\n";

        }
        
        //cout << "No UKF update: x_= \n" << x_ << endl;
        //cout << "No UKF update: P_= \n" << P_ << endl;
    }
    
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
  
   Make sure to update F & Q before running Kalman equations
   
   
   */
    
    //----------
    // # 1. - Generate Sigma Points (augmented with process noise)
    //----------
    // CTRV Model !! 5 component state vector turning rate and velocity are constant. Augmented state is 7 compoentns
    // Create Augmented State mean
    x_aug_.head(n_x_) = x_;                    // x_ is the previous state (before the new measurement update that has just come in)
    //x_aug_(n_aug_-2) = std_a_*std_a_;          // Assume acceleration process noise is fixed
    //x_aug_(n_aug_-1) = std_yawdd_*std_yawdd_ ; // Assume angle acceleration process noise is fixed
    x_aug_(n_aug_-2) = 0;          // Assume acceleration process noise is fixed
    x_aug_(n_aug_-1) = 0; // Assume angle acceleration process noise is fixed
    
    std::cout << "x_aug=\n" << x_aug_ << "\n";
    
    //create augmented covariance matrix
    P_aug_.fill(0.0);                               // Initialize full matrix to zero
    P_aug_.topLeftCorner(n_x_, n_x_) = P_;             // Standard P covariance matx is top-left
    P_aug_(n_aug_-2, n_aug_-2) = std_a_*std_a_;          // Augment bot-right Q matrix w/ variances
    P_aug_(n_aug_-1, n_aug_-1) = std_yawdd_*std_yawdd_;  // Augment bot-right Q matrix w/ variances
    std::cout << "P_aug=\n" << P_aug_ << "\n";
    
    //create square root of augmented covariance (P) matrix
    MatrixXd L = P_aug_.llt().matrixL();
    std::cout << "L=\n" << L << "\n";
    
    //---
    // Step 1. - Create Sigma Points (augmented)
    //---
    // 1st column of sigma point matrix
    Xsig_aug_.col(0)  = x_aug_;
    
    // Set remaining sigma points
    for (int i = 0; i < n_aug_; i++)
    {
        Xsig_aug_.col(i+1)        = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);
        Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
    }
    
    //cout << "Xsig_aug=\n" << Xsig_aug_ << "\n";
    PrintEigenMatrix("Xsig_aug=",&Xsig_aug_);
    
    /*
    cout << "Xsig_aug=\n";
    for (int i=0; i< Xsig_aug_.rows(); i++) {
        for (int j=0; j<Xsig_aug_.cols(); j++) {
            cout << setw(6) << setprecision(3) << Xsig_aug_(i,j);
        }
        cout << endl;
    }
    */
    
    
    
    //----------
    // #2. - Predict Sigma Points
    //----------
    
    //double delta_t = 0.1; //time diff in sec
    //create matrix with predicted sigma points as columns
    
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
            
            cout<< "phidot not near zero & eps=" << phidot << " " << __DBL_EPSILON__ << "\n";
            // Process equations w/ sigma points
            Xsig_pred_(0,j) = px     + (v/phidot) * (sin(phi+(phidot*delta_t)) - sin(phi)) \
                                     + 0.5*delta_t2*cos(phi)*nua;
            Xsig_pred_(1,j) = py     + (v/phidot) * (-cos(phi+(phidot*delta_t)) + cos(phi)) \
                                     + 0.5*delta_t2*sin(phi)*nua;
            Xsig_pred_(2,j) = v      + 0 + delta_t*nua;
            Xsig_pred_(3,j) = phi    + phidot*delta_t + 0.5*delta_t2*nuphia;
            Xsig_pred_(4,j) = phidot + 0 + delta_t*nuphia;
            
        } else {
            
            // Process equations w/ sigma points when phidot zero
            cout<< "phidot near zero (!) & eps=" << phidot << " " << __DBL_EPSILON__ << "\n";
            Xsig_pred_(0,j) = px     + v*cos(phi)*delta_t + 0.5*delta_t2*cos(phi)*nua;
            Xsig_pred_(1,j) = py     + v*sin(phi)*delta_t + 0.5*delta_t2*sin(phi)*nua;
            Xsig_pred_(2,j) = v      + 0 + delta_t*nua;
            Xsig_pred_(3,j) = phi    + 0 + 0.5*delta_t2*nuphia;
            Xsig_pred_(4,j) = phidot + 0 + delta_t*nuphia;
            
        }
        
    }
    cout << "Xsig_pred=\n" << Xsig_pred_ << "\n";
    
    

    //----------
    // # 3 - Predict State Mean Vector (x_) & State Covariance Matrix (P_) w/ Unscented Kalman Filter Equations
    //----------
    
    // Calc Predicted State Mean (weighted)
    x_.fill(0.0); // Note: current x_ has already been captured in the sigman points, so just zeroing out for new sum
    for (int j=0; j<num_sigma_pts_; j++) {  //iterate over sigma points (cols of matrix)
        x_ +=  weights_(j) * Xsig_pred_.col(j);
    }
    
    //---
    // Step 3B. - Calc Predicted State Covariance Matrix
    //---
    P_.fill(0.0);  // Note: current P_ has already been captured in the sigman points, so just zeroing out for new sum

    VectorXd x_diff(num_sigma_pts_);
    x_diff.fill(0.0);
    for (int j = 0; j <num_sigma_pts_; j++) {  //iterate over sigma points
        
        // State Vector difference (Pred - current)
        x_diff = Xsig_pred_.col(j) - x_ ;
        
        // CHECK angle normalization?? (for phi)
        //while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        //while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        P_ += weights_(j) * x_diff * x_diff.transpose() ;
    }

    // Debug - This is Prediction of x_ & P_ given state before measurement update
    cout << "Predict x_=\n" << x_ << "\n";
    cout << "Predict P_=\n" << P_ << "\n";

    return;
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
    
    
    //---
    // Step 1. - Predict Measurement using predicted measurement & existing sigma points
    //---
    
    
    //---
    // Step 2. - Update State w/ Kalman Filter Equations
    //---
    
    /*
    // Update the state by using Kalman Filter equations
    z_pred = H_ * x_;
    y = z - z_pred;
    Ht = H_.transpose();
    PHt = P_ * Ht;
    S = H_ * PHt + R_;  //S = H_ * P_ * Ht + R_;
    Si = S.inverse();
    K = PHt * Si;       // PHt = P_ * Ht;
    
    // New estimates (x_,P_)
    x_ = x_ + (K * y);
   
    long x_size = x_.size();
    I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
    
    //Debug
    std::cout << "After Update x_=\n" << x_ << "\n";
    std::cout << "After Update P_=\n" << P_ << "\n";
    */
    
    
    
    
    
    
    //---
    // Step Prologue. - Calculate Lidar NIS
    //---
    
    
    
    return;
}  // UpdateLidar



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
    
    
    //---
    // Step 1. - Predict Measurement using predicted measurement & existing sigma points
    //---

    //create example matrix with predicted sigma points
    //MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
    
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
    
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z_);
    
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z_,n_z_);
    
   
    // i=row, j=col
    
    //int num_sigma_pts = 2 * n_aug + 1;
    for (int j=0; j<num_sigma_pts_; j++) {
        
        // Step 1. Transform Predicted Sigma Points into measurement space (using Sigma Points from Prediction Step!!!)
        // 1-A. State vector
        double px     = Xsig_pred_(0,j);
        double py     = Xsig_pred_(1,j);
        double v      = Xsig_pred_(2,j);
        double yaw    = Xsig_pred_(3,j);
        //double yawdot = Xsig_pred_(4,j);
        
        double v1 = v*cos(yaw);
        double v2 = v*sin(yaw);
        
        // 1-B. Measurement model
        // NOTE: Need to check PX zero and px=py=zero !!!!!
        Zsig(0,j) = sqrt(px*px + py*py);                                    //rho
        Zsig(1,j) = atan2(py,px);                                           //phi NOTE: Need to add zero check
        Zsig(2,j) = (px*v1 + py*v2) / sqrt(px*px + py*py);  //rhodot
    }

    //---
    // Step 2. Calculate mean predicted measurement
    //---
    z_pred.fill(0.0);
    for (int j=0; j<num_sigma_pts_; j++) {
        z_pred += weights_(j) * Zsig.col(j);
    }
        
    //---
    // Step 3. Calculate measurement covariance matrix S
    //---
    S.fill(0.0);
    for (int j=0; j<num_sigma_pts_; j++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(j) - z_pred;
            
        //Note: Check angle normalization
        //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
            
        S += (weights_(j) * z_diff * z_diff.transpose()) ;
    }
    
    // Add "measurement noise" to covariance matrix once at end (its invarient)
    S = S + R_; // NOW ITS RIGHT. JUST ONCE ??????? IS this right - at every sigma point this gets updated????????
    
    //---
    // Step 4. - Update State w/ RADAR using  Kalman Filter Equations
    //---
    
    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z_);
    
    //---
    // Step 4A. - Calculate cross correlation matrix
    //---
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int j = 0; j < num_sigma_pts_; j++) {  //2n+1 simga points
        
        //residual
        VectorXd z_diff = Zsig.col(j) - z_pred;
        
        //angle normalization check
        //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(j) - x_;
        
        //angle normalization check
        //while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        //while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        Tc += weights_(j) * x_diff * z_diff.transpose();
    }
    
    //---
    // Step 4B. - Calculate Kalman gain K
    //---
    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    //---
    // Step 4C. - Update State Mean and Covariance matrix of the Unscented Kalman Filter equations
    //---
    
    //create example vector for incoming radar measurement
    VectorXd z = VectorXd(n_z_);
    // BIG EXPERIMENT
    //z << Zsig(0,0), Zsig(1,0) , Zsig(2,0); // 1st Sigma Point is incoming Radar point converted to rho, phi, phidot
    
    z = meas_package.raw_measurements_;  //Try using latest raw input measurement for z
    
    //residual
    VectorXd z_diff = z - z_pred;  // <TODO> WHAT IS Z ????????????
    cout << "z,z_pred,diff=" << z << " " << z_pred << " " << z_diff << "\n";
    
    
    //angle normalization
    //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    // Update State Mean & Covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
    
    
    //---
    // Step Prologue. - Calculate Radar NIS
    //---
    double radar_NIS = z_diff.transpose() * S.inverse() * z_diff;
    cout << "Radar NIS=" << radar_NIS << endl;
    
    return;
    
}  // UpdateRadar

//
// Utility Method to print a MatrixXd to cout
//
void UKF::PrintEigenMatrix(string label, MatrixXd *mat) {

    cout << label << endl;
    for (int i=0; i< mat->rows(); i++) {
        for (int j=0; j< mat->cols(); j++) {
            cout << setw(6) << setprecision(3) << (*mat)(i,j);
        }
        cout << endl;
    }
    return;
}
