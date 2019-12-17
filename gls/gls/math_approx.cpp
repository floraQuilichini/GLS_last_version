#include "stdafx.h"
#include "math_approx.h"

 tanh_lookup_table::tanh_lookup_table(Scalar min_abscisse, Scalar max_abscisse, int nb_samples)
{
	 min_abscisse_ = min_abscisse;
	 max_abscisse_ = max_abscisse;
	 nb_samples_ = nb_samples;

	step_ = (max_abscisse - min_abscisse) / (nb_samples - 1);

	lookup_table_.reserve(nb_samples);
	for (int i = 0; i < nb_samples; i++)
	{
		Scalar current_abscisse = min_abscisse + i*step_;
		Scalar tanh_value = tanh(current_abscisse);
		lookup_table_[i] = tanh_value;
	}
}


 void tanh_lookup_table::read_lookup_table(std::vector<Scalar>& abscisses, std::vector<Scalar>& tanh_values)
 {
	 size_t size_vec = abscisses.size();
	 tanh_values.reserve(size_vec);

	 // convert abscisse values to int
	 for (int i = 0; i < size_vec; i++)
	 {
		 Scalar value = abscisses[i];
		 if (value < min_abscisse_)
			 tanh_values[i] = -1;
		 else
		 {
			 if (value > max_abscisse_)
				 tanh_values[i] = 1;
			 else
			 {
				 int corresponding_int = round((value - min_abscisse_) / step_);
				 tanh_values[i] = corresponding_int;
			 }
		 }
	}
 }