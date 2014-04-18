/*
* Copyright (c) 2009, 2010, 2011, 2012 Brendon J. Brewer.
*
* This file is part of DNest3.
*
* DNest3 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DNest3 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with DNest3. If not, see <http://www.gnu.org/licenses/>.
*/

#include "MyModel.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include "Data.h"
#include <cmath>
#include <gsl/gsl_sf_gamma.h>

using namespace std;
using namespace DNest3;

const Data& MyModel::data = Data::get_instance();

MyModel::MyModel()
:bursts(4, 100, false, ClassicMassInf1D(data.get_t_min(), data.get_t_max(),
				1E-3*data.get_y_mean(), 1E3*data.get_y_mean()))
,mu(data.get_t().size())
{

}

void MyModel::calculate_mu()
{
	const vector< vector<double> >& components = bursts.get_components();
	const vector<double>& t = data.get_t();

	double scale;

	// Put in the non-delayed image
	for(size_t i=0; i<mu.size(); i++)
	{
		mu[i] = background;
		for(size_t j=0; j<components.size(); j++)
		{
			scale = components[j][2];
			if(t[i] > components[j][0])
				scale *= components[j][3];

			mu[i] += components[j][1]
					*exp(-abs(t[i] - components[j][0])/
						scale);
		}
	}

	// Repeat, with the delayed image
	for(size_t i=0; i<mu.size(); i++)
	{
		for(size_t j=0; j<components.size(); j++)
		{
			scale = components[j][2];
			if(t[i] > components[j][0])
				scale *= components[j][3];

			mu[i] += mag_ratio*components[j][1]
				*exp(-abs(t[i] -
					(components[j][0] - time_delay)
					)/scale);
		}
	}
}

void MyModel::fromPrior()
{
	background = exp(log(1E-3) + log(1E3)*randomU())*data.get_y_mean();
	bursts.fromPrior();

	double t_range = data.get_t_max() - data.get_t_min();
	time_delay = t_range*(randomU() - 0.5);
	mag_ratio = exp(randn());

	calculate_mu();
}

double MyModel::perturb()
{
	double logH = 0.;

	if(randomU() <= 0.5)
	{
		int which = randInt(3);

		if(which == 0)
		{
			background = log(background/data.get_y_mean());
			background += log(1E3)*pow(10., 1.5 - 6.*randomU())
						*randn();
			background = mod(background - log(1E-3), log(1E3))
						+ log(1E-3);
			background = exp(background)*data.get_y_mean();
		}
		else if(which == 1)
		{
			double t_range = data.get_t_max() - data.get_t_min();
			time_delay /= t_range;
			time_delay += randh();
			time_delay = mod(time_delay + 0.5, 1.) - 0.5;
			time_delay *= t_range;
		}
		else
		{
			mag_ratio = log(mag_ratio);
			logH -= -0.5*pow(mag_ratio, 2);
			mag_ratio += randh();
			logH += -0.5*pow(mag_ratio, 2);
			mag_ratio = exp(mag_ratio);
		}
	}
	else
	{
		logH += bursts.perturb();
	}

	calculate_mu();

	return logH;
}

double MyModel::logLikelihood() const
{
	const vector<int>& y = data.get_y();

	double logl = 0.;
	for(size_t i=0; i<y.size(); i++)
		logl += -mu[i] + y[i]*log(mu[i]) - gsl_sf_lngamma(y[i] + 1.);

	return logl;
}

void MyModel::print(std::ostream& out) const
{
	out<<background<<' '<<time_delay<<' '<<mag_ratio<<' ';
	bursts.print(out);
	for(size_t i=0; i<mu.size(); i++)
		out<<mu[i]<<' ';
}

string MyModel::description() const
{
	return string("background, time delay, mag ratio");
}

