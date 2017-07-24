#include "BMP.h"

#include <string>
#include <cstring>
#include <vector>
#include <cstdlib>

#include <stdint.h>
#include <complex>
#include <cmath>

using namespace std;

const int bmp_width = 1000;
const int bmp_height = 1000;
const double bmp_width_2 = bmp_width / 2.0;
const double bmp_height_2 = bmp_height / 2.0;
const double aspect = bmp_width / bmp_height;

const complex<double> center(-0.5, 0.0);
const double zoom = 4.0;
const double exp_zoom = exp(zoom * log(1.1));

const int max_calc = 500;

BMP mandelbrot(bmp_width, bmp_height);
BMP buddahbrot(bmp_width, bmp_height);


vector <complex<double> > z_samples(max_calc);
vector <long double> c_samples; // contributions


complex<double> current_orbit[max_calc];
int current_orbit_length = 0;

long double complex_distance(complex<double> a, complex<double> b) {
	return ((a.real() - b.real())*(a.real() - b.real()) + (a.imag() - b.imag())*(a.imag() - b.imag()));
}

long double random_range(double low, double high) {
	long double s = fabs(low - high);
	long double f = rand() / (double)RAND_MAX;
	return low + (s * f);
}

void random_complex_range(complex<double> &c, long double r) {
	while (1) {
		c = *(new complex<double>(random_range(-r, r), random_range(-r, r)));
		
		if ((c.real()*c.real() + (c.imag()*c.imag())) < (r*r)) {
			return;	
		}
	}
}

void complex_to_screen(complex<double> value, int &x, int &y) {
	x = (((value.real() - center.real()) / exp_zoom) + 1.0) / bmp_width_2;
	y = ((((value.imag() * aspect) - center.imag()) / exp_zoom) + 1.0) / bmp_height_2;
} 

bool on_screen(int x, int y) {
	return ((x >= 0 && x < bmp_width) && (y >= 0 && y < bmp_height));
}

bool evaluate(complex<double> z, int max_calc_local) {
	current_orbit_length = 0;
	
	complex<double> zc(z);
	for (int i = 0; i < max_calc_local; i++) {
		printf("evaluate z=(%0.2f, %0.2f) %d\n", z.real(), z.imag(), i);
		
		z = (z * z) + zc;
		if ((z.real()*z.real() + z.imag()*z.imag()) > 4.0) {
			return true;
		}
		
		current_orbit[current_orbit_length++] = complex<double>(z);
		
		if (current_orbit_length >= max_calc_local) {
			return false;
		}
	}
	return false;
}

long double get_contribution() {
	int x, y;
	
	long double contribution = 0;

	for (int i = 0; i < current_orbit_length; i++) {
		complex_to_screen(current_orbit[i], x, y);
		
		if (on_screen(x, y)) {
			contribution++;
		}
	}

	if (current_orbit_length != 0) {
		return contribution / (long double)(current_orbit_length);
	} else {
		return 0;
	}
}

bool find_initial_sample(complex<double> &c, long double x, long double y, long double rad, int f) {
	printf("find_initial_sample(c=(%0.2f, %0.2f), x=%0.2f, y=%0.2f, rad=%0.2f, f=%d)\n", c.real(), c.imag(), x, y, rad, f);
	if (f > 500)
		return false;
	
	complex<double> tmp, seed;
	
	long double closest = 1e20;
	
	for (int i = 0; i < 200; i++) {
		random_complex_range(tmp, rad);
		tmp += *(new complex<double>(x, y));
		printf("tmp=(%0.2f, %0.2f) ", tmp.real(), tmp.imag());
		
		if (!evaluate(tmp, 500)) {
			printf("skipped\n");
			continue;
		}
		printf("evaluated ");
		
		if (get_contribution() > 0.0) {
			c = tmp;
			printf("contributed\n");
			return true;
		}
		
		printf("NOT contributed\n");
		
		for (int q = 0; q < current_orbit_length; q++) {
			long double d = complex_distance(current_orbit[q], center);
			if (d < closest) {
				closest = d;
				seed = tmp;
			}
		}
	}
	return find_initial_sample(c, seed.real(), seed.imag(), rad / 2.0, f+1);
}

void build_initial_sample_points() {
	for (int i = 0; i < 30; i++) {
		complex<double> m(0.0, 0.0);
		
		printf("current sample %d\n", i);
		if (!find_initial_sample(m, 0, 0, 2.0, 0)) {
			printf("couldn't find seed %d\n", i);
			continue;
		}
		printf("m=(%0.2f, %0.2f)\n", m.real(), m.imag());
		
		evaluate(m, 500);
		z_samples.push_back(m);
		c_samples.push_back(get_contribution());
	}
}

int main() {
	build_initial_sample_points();
	
	mandelbrot.save("test.bmp");
	
	return 0;
}