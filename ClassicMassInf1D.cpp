#include "ClassicMassInf1D.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cmath>

using namespace DNest3;

ClassicMassInf1D::ClassicMassInf1D(double x_min, double x_max,
					double amp_min, double amp_max)
:x_min(x_min)
,x_max(x_max)
,amp_min(amp_min)
,amp_max(amp_max)
{

}

void ClassicMassInf1D::fromPrior()
{
	min_amp = exp(log(amp_min) + log(amp_max/amp_min)*randomU());
	min_width = exp(log(1E-3*(x_max - x_min)) + log(1E3)*randomU());

	alpha_amp = 0.1 + 4.9*randomU();
	alpha_width = 0.1 + 4.9*randomU();

	a = -10. + 10.*randomU();
	b = 2.*randomU();
}

double ClassicMassInf1D::perturb_parameters()
{
	double logH = 0.;

	int which = randInt(6);

	if(which == 0)
	{
		min_amp = log(min_amp);
		min_amp += log(amp_max/amp_min)*pow(10., 1.5 - 6.*randomU())*randn();
		min_amp = mod(min_amp - log(amp_min), log(amp_max/amp_min)) + log(amp_min);
		min_amp = exp(min_amp);
	}
	if(which == 1)
	{
		min_width = log(min_width/(x_max - x_min));
		min_width += log(1E3)*pow(10., 1.5 - 6.*randomU())*randn();
		min_width = mod(min_width - log(1E-3), log(1E3)) + log(1E-3);
		min_width = (x_max - x_min)*exp(min_width);
	}
	if(which == 2)
	{
		a += 20.*randh();
		a = mod(a + 10., 20.) - 10.;
	}
	if(which == 3)
	{
		b += 2.*randh();
		b = mod(b, 2.);
	}

	return logH;
}

double ClassicMassInf1D::log_pdf(const std::vector<double>& vec) const
{
	if(vec[0] < x_min || vec[0] > x_max || vec[1] < 0. || vec[2] < 0.
		|| log(vec[3]) < (a-b) || log(vec[3]) > (a + b))
		return -1E300;

	return -log(mu) - vec[1]/mu - log(mu_widths) - vec[2]/mu_widths
			- log(2.*b*vec[3]);
}

void ClassicMassInf1D::from_uniform(std::vector<double>& vec) const
{
	vec[0] = x_min + (x_max - x_min)*vec[0];
	vec[1] = -mu*log(1. - vec[1]);
	vec[2] = -mu_widths*log(1. - vec[2]);
	vec[3] = exp(a - b + 2.*b*vec[3]);
}

void ClassicMassInf1D::to_uniform(std::vector<double>& vec) const
{
	vec[0] = (vec[0] - x_min)/(x_max - x_min);
	vec[1] = 1. - exp(-vec[1]/mu);
	vec[2] = 1. - exp(-vec[2]/mu_widths);
	vec[3] = (log(vec[3]) + b - a)/(2.*b);
}

void ClassicMassInf1D::print(std::ostream& out) const
{
	out<<mu<<' '<<mu_widths<<' '<<a<<' '<<b<<' ';
}

