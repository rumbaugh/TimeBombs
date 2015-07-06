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

#ifndef _MyModel_
#define _MyModel_

#include "dnest3/Model.h"
#include "Data.h"
#include "RJObject.h"
#include "ClassicMassInf1D.h"
#include <vector>

class MyModel:public DNest3::Model
{
	private:
		// Reference to the data
		static const Data& data;

		// A flat background level
		double background;

		// The bursts
		RJObject<ClassicMassInf1D> bursts;

		// Time delay
		double time_delay;

		// Magnification ratio
		double mag_ratio;

		// Poisson mean
		std::vector<double> mu;

		// Calculates mu from scratch
		void calculate_mu();

		// Calculates mu from scratch, ignoring second image
		void calculate_mu_single();

	public:
		MyModel();

		// Generate the point from the prior
		void fromPrior();

		// Metropolis-Hastings proposals
		double perturb();

		// Likelihood function
		double logLikelihood() const;

		// Print to stream
		void print(std::ostream& out) const;

		// Return string with column information
		std::string description() const;
};

#endif

