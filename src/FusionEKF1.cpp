#include "FusionEKF.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor
 */
FusionEKF::FusionEKF() {
	is_initialized_ = false;
	previous_timestamp_ = 0;

	R_laser_ = MatrixXd(2, 2);
	R_radar_ = MatrixXd(3, 3);
	H_laser_ = MatrixXd(2, 4);
	Hj_ = MatrixXd(3, 4);

	// measurement covariance matrix - laser
	R_laser_ << 2, 0,
							0, 2;

	// measurement covariance matrix - radar
	R_radar_ << 0.1, 0, 0,
							0, 0.001, 0,
							0, 0, 0.1;

	H_laser_ << 1, 0, 0, 0,
							0, 1, 0, 0;

	ekf_.F_ = MatrixXd(4, 4);

	ekf_.F_ << 1, 0, 1, 0,
						 0, 1, 0, 1,
						 0, 0, 1, 0,
						 0, 0, 0, 1;
}

FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

	/*****************************************************************************
	 *  Prediction
	 ****************************************************************************/
	if (!is_initialized_) {
		// first measurement
		cout << "EKF: " << endl;
		ekf_.x_ = VectorXd(4);
		ekf_.x_ << 1, 1, 1, 1;

		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
		{
			ekf_.x_(0) = measurement_pack.raw_measurements_(0) * cos(measurement_pack.raw_measurements_(1));
			ekf_.x_(1) = measurement_pack.raw_measurements_(0) * sin(measurement_pack.raw_measurements_(1));
		}
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
		{
			ekf_.x_(0) = measurement_pack.raw_measurements_(0);
			ekf_.x_(1) = measurement_pack.raw_measurements_(1);
		}

		ekf_.P_ = MatrixXd(4, 4);
		ekf_.P_ << 1, 0, 0, 0,
							 0, 1, 0, 0,
							 0, 0, 1, 0,
							 0, 0, 0, 1;

		previous_timestamp_ = measurement_pack.timestamp_;
		is_initialized_ = true;

		// done initializing, no need to predict or update
		return;
	}

	/*****************************************************************************
	 *  Prediction
	 ****************************************************************************/
	// compute the time elapsed between the current and previous measurements
	// delta_t - expressed in seconds
	float delta_t = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
	float delta_t2 = delta_t * delta_t;
	float delta_t3 = delta_t2 * delta_t;
	float delta_t4 = delta_t3 * delta_t;
	previous_timestamp_ = measurement_pack.timestamp_;

	
	ekf_.F_(0, 2) = delta_t;
	ekf_.F_(1, 3) = delta_t;


	//Process covariance matrix
	//(17,7)
	//1021,511 good for text-1
	float sigma_ax = 1021;//13;
	float sigma_ay = 511;//21;
	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ << (delta_t4 / 4 )* sigma_ax, 0, delta_t3 / 2 * sigma_ax, 0,
						 0, (delta_t4 / 4 )* sigma_ay, 0, delta_t3 / 2 * sigma_ay,
						 (delta_t3 / 2) * sigma_ax, 0, delta_t2 * sigma_ax, 0,
						 0, (delta_t3 / 2) * sigma_ay, 0, delta_t2 * sigma_ay;

	ekf_.Predict();

	/*****************************************************************************
	 *  Update
	 ****************************************************************************/
	/* TODO:
	 * 1. Define H, R for each sensor
	 * 2. Change H, R, based on the sensor type
	 * 3. Apply KF
	 * 4. Define timestamps
	 */

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
	{
		// state parameters
		float x = ekf_.x_(0);
		float y = ekf_.x_(1);

		if (x == 0 && y == 0)
		{
			// avoid NaN values, just predict the current state
			return;
		}

		MatrixXd Hj = tools.CalculateJacobian(ekf_.x_);

		ekf_.H_ = Hj;
		ekf_.R_ = R_radar_;

		ekf_.UpdateEKF(measurement_pack.raw_measurements_,ekf_);
	}
	else
	{
		ekf_.H_ = H_laser_;
		ekf_.R_ = R_laser_;
		ekf_.Update(measurement_pack.raw_measurements_);
	}

	// print the output
	cout << "x_= " << ekf_.x_ << endl;
	cout << "P_= " << ekf_.P_ << endl;
}
