#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF
{
public:
  bool is_initialized_;  // set to true in first call of ProcessMeasurement
  bool use_laser_;
  bool use_radar_;

  int n_x_;             // State dimension
  int n_xp1_;           // n_x_ + 1
  VectorXd x_;          // [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd dx_;         // temporary vector
  MatrixXd P_;          // state covariance matrix

  // sigma points
  int n_aug_;           // Augmented state dimension
  int n_s_;             // augmented sigma points
  double lambda_;       // Sigma point spreading parameter
  double term_c_;       // constant term
  VectorXd x_aug_;      // augmented state
  MatrixXd P_aug_;      // augmented covariance matrix
  MatrixXd Xsig_aug_;   // augmented sigma point matrix
  MatrixXd Xsig_pred_;  // predicted sigma points matrix
  VectorXd weights_;    // Weights of sigma points


  long long time_us_;   // time when the state is true, in us

  double std_a_;        // Process noise std longitudinal acceleration in m/s^2
  double std_yawdd_;	// Process noise std yaw acceleration in rad/s^2

  double std_laspx_;    // Laser measurement noise standard deviation position1 in m
  double std_laspy_;    // Laser measurement noise standard deviation position2 in m
  double std_radr_;     // Radar measurement noise standard deviation radius in m
  double std_radphi_;   // Radar measurement noise standard deviation angle in rad
  double std_radrd_ ;   // Radar measurement noise standard deviation radius change in m/s

  MatrixXd H_laser_;    // LIDAR observations
  MatrixXd H_laser_t_;  // transpose
  MatrixXd R_laser_;    // LIADR measurement covariance
  Eigen::MatrixXd I_;   // identity

  MatrixXd Zsig_;       // radar sigma points matrix
  VectorXd z_;          // radar observation
  MatrixXd S_;          // radar sigma points matrix
  MatrixXd Tc_;         // radar cross correlation matrix
  MatrixXd R_radar_;    // radar observation noise covariance

  UKF();
  virtual ~UKF();

  // param meas_package is the latest measurement data of either radar or laser
  void ProcessMeasurement(MeasurementPackage meas_package);

private:
  // Predicts sigma points, x, P. delta_t Time between k and k+1 in s
  void Prediction(double delta_t);

  // Updates x and P using a laser measurement, param meas_package at k+1
  void UpdateLidar(MeasurementPackage meas_package);

  // Updates x and P using a radar measurement, param meas_package at k+1
  void UpdateRadar(MeasurementPackage meas_package);

  // first measurement
  void initialize(const MeasurementPackage &measurement_pack);

  // normalise angle
  void norm_angle(double &theta);
  void norm_angle(VectorXd &vec, int index);

  long count;
};

#endif
