#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}


VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check estimations not empty
	if (estimations.size() == 0)
	{
	    cout << "estimation vector empty" << endl;
	    return(rmse);
	}

	//  check estimation  size equals ground truth size
	if (estimations.size() != ground_truth.size())
	{
	    cout << "estimation vector size does not match ground_truth" << endl;
	    return(rmse);
	}

    VectorXd error(4);
	for(int i=0; i < estimations.size(); ++i)
	{
	    error = estimations[i]-ground_truth[i];
	    error = error.array()*error.array();
        rmse = rmse + error;
	}
	rmse = rmse / float(estimations.size());
	rmse = rmse.array().sqrt();
	//cout << "mse = " << rmse << endl;
	return(rmse);
}
