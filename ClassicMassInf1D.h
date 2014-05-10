#ifndef _ClassicMassInf1D_
#define _ClassicMassInf1D_

#include "Distributions/Distribution.h"

class ClassicMassInf1D:public Distribution
{
	private:
		// Limits
		double x_min, x_max;
		double amp_min, amp_max;

		// Lower limits for amplitudes and widths
		double min_amp, min_width;

		// Pareto slopes for amplitudes and widths
		double alpha_amp, alpha_width;

		// Uniform for log-skews
		double a, b; // Midpoint and half-width

		double perturb_parameters();

	public:
		ClassicMassInf1D(double x_min, double x_max,
					double amp_min, double amp_max);

		void fromPrior();

		double log_pdf(const std::vector<double>& vec) const;
		void from_uniform(std::vector<double>& vec) const;
		void to_uniform(std::vector<double>& vec) const;

		void print(std::ostream& out) const;
};

#endif

