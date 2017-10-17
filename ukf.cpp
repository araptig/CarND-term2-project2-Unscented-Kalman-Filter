#include "ukf.h"
#include "Eigen/Dense"
//#include <math.h>
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

UKF::UKF()
{//constructor

  {// initialize boolean flags
    use_laser_   	= true;            	// false ==> ignore laser measurements
    use_radar_ 		= true;            	// false ==> ignore lidar measurements
    is_initialized_ = false;      		// not initialized yet
    count        	= 0;				// for debugging
  }

  {// initialize x & P
    n_x_   = 5;							// state is of size 5
    n_xp1_ = n_x_ + 1;					// for convenience
    x_     = VectorXd(n_x_);            // state vector
    dx_    = VectorXd(n_x_);            // temporary vector
    P_     = MatrixXd(n_x_, n_x_);      // covariance
  }


  {// process noise: ****OPTIMIZE****
 	  std_a_     = .75;           		// std longitudinal acceleration in m/s^2
 	  std_yawdd_ = M_PI/6;      		// std yaw acceleration in rad/s^2
 	  Q_          = MatrixXd(2, 2);
 	  Q_.fill(0.0);
 	  Q_(0,0)     = std_a_*std_a_;
 	  Q_(1,1)     = std_yawdd_* std_yawdd_;
   }

  {// initialize augmented x & P
	  n_aug_ = n_x_ + 2; ;               // augmented state dimension (with process noise)
	  x_aug_ = VectorXd(n_aug_);         // augmented state vector
	  x_aug_.fill(0.0);

	  P_aug_ = MatrixXd(n_aug_,n_aug_);  // augmented covariance
	  P_aug_.fill(0.0);
	  P_aug_.bottomRightCorner(2,2)      = Q_;
  }

  {//initialize sigma point matrices
    n_s_    = 2*n_aug_ +1;             // sigma points
    lambda_ = 3 - n_aug_;              // hyper parameter
    term_c_ = sqrt(lambda_+n_aug_);    // precomputed entity for sigma point generation

    Xsig_aug_  = MatrixXd(n_aug_, n_s_);
    Xsig_pred_ = MatrixXd(n_x_, n_s_);
  }

  {// set weights
	  weights_       = VectorXd(n_s_);
	  double temp    = 1.0/ double(lambda_ + n_aug_);
      weights_(0)    = lambda_ * temp;
      temp           = 0.5*temp;
      for (int i=1; i<n_s_; i++)
      {
          weights_(i) = temp;
      }
  }


  {// set laser
	  int n_l = 2;                      // observations
	  H_laser_ = MatrixXd(n_l, n_x_);	// observation for laser
	  H_laser_ << 1, 0, 0, 0, 0,
	 		      0, 1, 0, 0, 0;
	  H_laser_t_ = H_laser_.transpose();// for convenience

	  R_laser_ = MatrixXd(n_l, n_l);    // observation noise covariance
	  std_laspx_  = 0.15;  				// noise std position1 in m
	  std_laspy_  = 0.15;  				// noise std position2 in m
	  R_laser_ << std_laspx_*std_laspx_, 0,
	 	 		      0,                 std_laspy_*std_laspy_;
	  I_ = MatrixXd::Identity(n_x_,n_x_); // for convenience
  }

  {// set radar
	  int n_r  = 3;						// observations
	  Zsig_    = MatrixXd(n_r, n_s_);   // sigma points
	  z_       = VectorXd(n_r);         // observation vector
	  S_       = MatrixXd(n_r, n_r);    // temporary matrix
	  Tc_      = MatrixXd(n_x_, n_r);   // temporary matrix

	  std_radr_   = 0.3;   // Radar measurement noise std radius in m
	  std_radphi_ = 0.03;  // Radar measurement noise std angle in rad
	  std_radrd_ = 0.3;    // Radar measurement noise std radius change in m/s
	  R_radar_ = MatrixXd(n_r, n_r);
	  R_radar_.fill(0.0);
	  R_radar_(0,0) = (std_radr_*std_radr_);
	  R_radar_(1,1) = (std_radphi_*std_radphi_);
	  R_radar_(2,2) = (std_radrd_*std_radrd_);
  }
}

UKF::~UKF() {}

// new measurement available
void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{// ProcessMeasurement
  // Complete this function! Make sure you switch between lidar and radar measurements.

  // initialization
  if (!is_initialized_)
  {//initialize
	initialize(meas_package);
    return;
  }//initialize

  {// prediction
	  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
	  Prediction(dt);
	  time_us_ = meas_package.timestamp_;                  // set new reference
  }// prediction

  {//update
	  if ( (meas_package.sensor_type_ == MeasurementPackage::RADAR) && (use_radar_==true))
	  {//radar
		  cout << "radar" << endl;
		  UpdateRadar(meas_package);
	  }//radar

	  if ( (meas_package.sensor_type_ == MeasurementPackage::LASER) && (use_laser_==true))
	  {// laser
		  cout << "laser" << endl;
		  UpdateLidar(meas_package);
	  }// laser
  }//update

  {// print
	  count ++;
	  cout << "x_ = " <<  double(x_(0)) << " "
			          <<  double(x_(1)) << " "
					  <<  double(x_(2)) << " "
					  <<  double(x_(3)) << " "
					  <<  double(x_(4)) << endl;
	  cout << "P_ = " << double(P_(0,0)) << " "
			          << double(P_(1,1)) << " "
					  << double(P_(2,2)) << " "
					  << double(P_(3,3)) << " "
					  << double(P_(4,4)) << endl;
  }
}// ProcessMeasurement


void UKF::Prediction(double delta_t)
{// prediction

  // (I) generate sigma points

  {//set augmented x and P
      x_aug_.head(5)                     = x_;
      P_aug_.topLeftCorner(n_x_,n_x_)    = P_;
  }

  {// sigma-points
	  //create square root matrix
	  MatrixXd A = P_aug_.llt().matrixL();

	  //create augmented sigma points
	  Xsig_aug_.col(0) = x_aug_;
	  for (int i = 0; i < n_aug_; i++)
	  {
		  VectorXd temp = term_c_ * A.col(i);
		  Xsig_aug_.col(i+1)        = x_aug_ + temp;
		  Xsig_aug_.col(i+1+n_aug_) = x_aug_ - temp;
	  }
  }

  // (II) predict sigma points @ (k+1)
  double dt2 = 0.5 * delta_t * delta_t;
  for (int i=0; i<n_s_; i++)
  {// for each sigma point
        VectorXd x_aug = Xsig_aug_.col(i);			// pick a sigma point
        double x_pos = Xsig_aug_(0,i);
        double y_pos = Xsig_aug_(1,i);
        double vel   = Xsig_aug_(2,i);
        double psi   = Xsig_aug_(3,i);
        double psi_d = Xsig_aug_(4,i);
        double n_a   = Xsig_aug_(5,i);
        double n_p   = Xsig_aug_(6,i);

        double cos_p =  cos(psi);
        double sin_p =  sin(psi);

        dx_(0) = dt2 * cos_p* n_a;  // noise
        dx_(1) = dt2 * sin_p* n_a;  // noise
        dx_(2) = delta_t * n_a;
        dx_(3) = delta_t * psi_d + dt2*n_p;
        dx_(4) = delta_t * n_p;

        if(fabs(psi_d) < .002)
        {
            dx_(0) += delta_t * vel * cos_p;
            dx_(1) += delta_t * vel * sin_p;
        }
        else
        {
           double temp_1 = vel/psi_d;
           double temp_2 = psi + delta_t*psi_d;
           norm_angle(temp_2);
           dx_(0) +=  temp_1 * (sin(temp_2) - sin_p);
           dx_(1) +=  temp_1 * (-cos(temp_2) + cos_p);
        }

        Xsig_pred_.col(i) = x_aug.head(n_x_) + dx_;
    }// for each sigma point

  // (III) predict x & P

  {//predicted state @ (k+1)
    x_ = weights_(0)*Xsig_pred_.col(0);
    for (int i=1; i<n_s_; i++)
    {
        x_ +=  weights_(i) * Xsig_pred_.col(i);
    }
  }

  {//predicted covariance @ (k+1)
    P_.fill(0);
    for (int i=0; i<n_s_; i++)
    {
        VectorXd temp = Xsig_pred_.col(i) - x_;
        norm_angle(temp,3);
        P_ +=  weights_(i) * temp * temp.transpose();
    }
  }
}// prediction


void UKF::UpdateLidar(MeasurementPackage meas_package)
{// UpdateLidar
	VectorXd y = meas_package.raw_measurements_ - H_laser_ * x_;
	MatrixXd S = H_laser_ * P_ * H_laser_t_ + R_laser_;
	MatrixXd K = P_ * H_laser_t_ * S.inverse();
	x_ = x_ + K * y;
	P_ = (I_ - K*H_laser_)*P_;
  //You'll also need to calculate the lidar NIS.
}// UpdateLidar


void UKF::UpdateRadar(MeasurementPackage meas_package)
{// UpdateRadar
	// Xsig_pred_ --> Zsig_
	for (int i=0; i<n_s_; i++)
	{
	  double rho  = sqrt(double(Xsig_pred_(0,i)*Xsig_pred_(0,i) + Xsig_pred_(1,i)*Xsig_pred_(1,i)));

	  double phi  = 0.0;
	  double rho_d= 0.0;
	  if(rho > 0.01)
	  {
	    phi   = atan2(double(Xsig_pred_(1,i)), double(Xsig_pred_(0,i)));
	    rho_d = (Xsig_pred_(0,i)*cos(double(Xsig_pred_(3,i))) + Xsig_pred_(1,i)*sin(double(Xsig_pred_(3,i))));
	    rho_d *= (Xsig_pred_(2,i)/rho);
	  }
	  Zsig_(0,i) = rho;
	  Zsig_(1,i) = phi;
	  Zsig_(2,i) = rho_d;
	}

	//calculate mean predicted measurement
	z_.fill(0);
	for (int i=0; i<n_s_; i++)
	{
	  z_ +=  weights_(i) * Zsig_.col(i);
	}

	//calculate measurement covariance matrix S
	S_.fill(0);
	for (int i=0; i<n_s_; i++)
	{
	      VectorXd temp = Zsig_.col(i) - z_;
	      norm_angle(temp,1);
	      S_ +=  weights_(i) * temp * temp.transpose();
	}
	// add observation covariance
	S_ += R_radar_;

	// set Tc
	Tc_.fill(0.0);
	for(int i=0; i<n_s_; i++)
	{
	      VectorXd temp_x = Xsig_pred_.col(i) - x_;
	      norm_angle(temp_x,3);

	      VectorXd temp_z = Zsig_.col(i) - z_;
	      norm_angle(temp_z,1);

	      Tc_ +=  weights_(i) * temp_x * temp_z.transpose();
	  }

	  //calculate Kalman gain K;
	  MatrixXd K = Tc_ * S_.inverse();

	  //update state mean and covariance matrix
	  VectorXd z_diff = meas_package.raw_measurements_ - z_;
	  norm_angle(z_diff,1);
	  x_ += (K*z_diff);
	  P_ -= (K*S_*K.transpose());
}// UpdateRadar


void UKF::norm_angle(double &theta)
{
	while (theta > M_PI)
		theta -= 2.*M_PI;

	while (theta < -M_PI)
		theta += 2.*M_PI;
}

void UKF::norm_angle(VectorXd &vec, int index)
{
	while (vec(index) > M_PI)
			vec(index) -= 2.*M_PI;

	while (vec(index) < -M_PI)
			vec(index) += 2.*M_PI;
}

void UKF::initialize(const MeasurementPackage &measurement_pack)
{//initialize
	x_.fill(0.0);
    P_.fill(0.0);

	if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
	{//laser
		x_(0) = measurement_pack.raw_measurements_(0);
		x_(1) = measurement_pack.raw_measurements_(1);
		P_(0,0) = R_laser_(0,0);
		P_(1,1) = R_laser_(1,1);
		cout << "laser" << endl;
	}//laser
	else
	{//radar
		float rho      = measurement_pack.raw_measurements_(0);
		double phi     = measurement_pack.raw_measurements_(1);
		norm_angle(phi);
		float rho_dot  = measurement_pack.raw_measurements_(2);
		x_(0)      = rho * cos(phi);
		x_(1)      = rho * sin(phi);
		P_(0,0)    = 0.3;
		P_(1,1)    = 0.3;
		cout << "radar" << endl;
	}//radar

	{// covariance rest
		P_(2,2) = 20;
		P_(3,3) = 2;
		P_(4,4) = 1;
	}

	is_initialized_ = true;
	time_us_ = measurement_pack.timestamp_;         // remember time
	cout << x_ << endl;
	cout << P_ << endl;
	cout << "UKF initialized: " << endl;
} // initialize


